/*********************************************************************************
* Copyright (C) 2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file NPI.cpp
  \brief This file implements the norm-preserving integrator of Meek and Levine
    
*/

#include "NPI.h"
#include "../math_specialfunctions/libspecialfunctions.h"
#include "../math_linalg/liblinalg.h"
#include "../math_meigen/libmeigen.h"

#include "../Units.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;
using namespace libmeigen;


/// libcalculators namespace
namespace libcalculators{

namespace {

const double NPI_ORTHOGONALITY_TOL = 1.0e-6;
const double NPI_SINGULARITY_TOL = 1.0e-12;

double clamp_unit_interval(double value){
  return std::max(-1.0, std::min(1.0, value));
}

double determinant(MATRIX matrix){
  double result = 1.0;
  int sign = 1;

  for(int col=0; col<matrix.n_cols; col++){
    int pivot = col;
    double pivot_abs = std::fabs(matrix.get(col, col));
    for(int row=col+1; row<matrix.n_rows; row++){
      double candidate = std::fabs(matrix.get(row, col));
      if(candidate > pivot_abs){
        pivot = row;
        pivot_abs = candidate;
      }
    }

    if(pivot_abs < NPI_SINGULARITY_TOL){ return 0.0; }
    if(pivot != col){
      for(int k=col; k<matrix.n_cols; k++){
        double tmp = matrix.get(col, k);
        matrix.set(col, k, matrix.get(pivot, k));
        matrix.set(pivot, k, tmp);
      }
      sign = -sign;
    }

    double pivot_value = matrix.get(col, col);
    result *= pivot_value;
    for(int row=col+1; row<matrix.n_rows; row++){
      double factor = matrix.get(row, col) / pivot_value;
      for(int k=col+1; k<matrix.n_cols; k++){
        matrix.set(row, k, matrix.get(row, k) - factor * matrix.get(col, k));
      }
    }
  }
  return sign * result;
}

void validate_npi_input(const MATRIX& St, double dt){
  if(St.n_rows != St.n_cols || St.n_rows == 0){
    std::ostringstream msg;
    msg << "nac_npi: the time-overlap matrix must be a non-empty square matrix; "
        << "received " << St.n_rows << " x " << St.n_cols
        << ". Supply overlaps between the same complete set of states at the two time steps.";
    throw std::invalid_argument(msg.str());
  }
  if(!std::isfinite(dt) || dt <= 0.0){
    std::ostringstream msg;
    msg << "nac_npi: dt must be finite and positive; received " << dt
        << ". Pass the positive time interval separating the two overlap matrices.";
    throw std::invalid_argument(msg.str());
  }

  double max_orthogonality_error = 0.0;
  for(int i=0; i<St.n_rows; i++){
    for(int j=0; j<St.n_cols; j++){
      double value = St.get(i, j);
      if(!std::isfinite(value)){
        std::ostringstream msg;
        msg << "nac_npi: overlap element (" << i << ", " << j << ") is not finite. "
            << "Check the electronic-structure calculation and state-overlap construction.";
        throw std::invalid_argument(msg.str());
      }
    }
    if(St.get(i, i) < -NPI_ORTHOGONALITY_TOL){
      std::ostringstream msg;
      msg << "nac_npi: diagonal overlap S(" << i << ", " << i << ") = " << St.get(i, i)
          << " is negative, indicating a discontinuous electronic-state phase. "
          << "Phase-match and reorder the states so corresponding-state overlaps are non-negative.";
      throw std::invalid_argument(msg.str());
    }
  }

  for(int i=0; i<St.n_cols; i++){
    for(int j=0; j<St.n_cols; j++){
      double dot = 0.0;
      for(int k=0; k<St.n_rows; k++){
        dot += St.get(k, i) * St.get(k, j);
      }
      double expected = (i == j) ? 1.0 : 0.0;
      max_orthogonality_error = std::max(max_orthogonality_error, std::fabs(dot - expected));
    }
  }
  if(max_orthogonality_error > NPI_ORTHOGONALITY_TOL){
    std::ostringstream msg;
    msg << "nac_npi: the time-overlap matrix is not orthogonal; max|S^T S - I| = "
        << max_orthogonality_error << " exceeds " << NPI_ORTHOGONALITY_TOL
        << ". Use the same complete state space at both steps and orthogonalize the overlap "
        << "matrix (for example, by polar/Lowdin orthogonalization) after state matching.";
    throw std::invalid_argument(msg.str());
  }

  double det = determinant(St);
  if(det <= 0.0){
    std::ostringstream msg;
    msg << "nac_npi: the phase-matched overlap matrix must represent a proper rotation, but det(S) = "
        << det << ". Correct state phases/permutations so det(S) is positive before applying NPI.";
    throw std::invalid_argument(msg.str());
  }
}

} // namespace

MATRIX nac_npi(MATRIX& St, double dt){
/** 
 Computes derivative coupling matrix elements using NPI of Meek and Levine

 Meek, G.A; Levine, B. G. "Evaluation of the Time-Derivative Coupling for Accurate
 Electronic State Transition Probabilities from Numerical Simulations" J. Phys. Chem. Lett.
 5, 2351-2356, 2014

*/

  validate_npi_input(St, dt);

  int nstates = St.n_cols;
  MATRIX nac(nstates, nstates);

  for(int i=0; i<nstates; i++){
    for(int j=i+1; j<nstates; j++){

// W_jk = <j|d/dt|k>      j->i;   k->j
// d_kj = ...

      double W00 = clamp_unit_interval(St.get(i,i));
      double W01 = clamp_unit_interval(St.get(i,j));
      double W10 = clamp_unit_interval(St.get(j,i));
      double W11 = clamp_unit_interval(St.get(j,j));

      double A = acos(W00) - asin(W01);
      double B = acos(W00) + asin(W01);
      double C = acos(W11) - asin(W10);
      double D = acos(W11) + asin(W10);

      // Not sure what is this
      //  if (Wlj != Wlj){
      //      Wlj = 0.0;}

      if(A == 0.0){ A = -1.0; }
      else{  A = -1.0 * sin(A) / A; }

      if(B == 0.0){  B = 1.0; }
      else{   B = sin(B) / B; }

      if(C == 0.0){  C = 1.0; }
      else{   C = sin(C) / C; }

      if(D == 0.0){  D = 1.0; }
      else{   D = sin(D) / D; }

      //cout << "Flag A:" << A << endl;
      //cout << "Flag B:" << B << endl;
      //cout << "Flag C:" << C << endl;
      //cout << "Flag D:" << D << endl;
      //cout << "Flag W00:" << W00 << endl;
      //cout << "Flag W10:" << W10 << endl;

      double Wlj = 1.0 - (W00 * W00) - (W10 * W10);
      if(Wlj < 0.0){ Wlj = 0.0; }
      Wlj = sqrt(Wlj);

      double E; 
      if(Wlj == 0.0){     E = 0.0;      }
      else{
         double Wlk = clamp_unit_interval(-1.0 * (W01 * W00 + W11 * W10) / Wlj);
         E = (1.0 - Wlj * Wlj) * (1.0 - Wlk * Wlk);
         if(E < 0.0 ){   E = 0.0; }
         else{

           double sWlj = asin(Wlj);
           double sWlk = asin(Wlk);
           //cout << "Flag Wlj:" << Wlj << " , sWlj: " << sWlj << endl;
           //cout << "Flag Wlk:" << Wlk << " , sWlk: " << sWlk << endl;
           E = sqrt(E); 
           double denom = sWlj * sWlj - sWlk * sWlk;
           //cout << "Flag NPI, denom:" << denom << endl;
           if(std::fabs(denom) <= NPI_SINGULARITY_TOL){
             E = (Wlk < 0.0) ? -Wlj * Wlj : Wlj * Wlj;
           }
           else{
             E = 2.0 * sWlj * (Wlj * Wlk * sWlj + (E - 1.0) * sWlk) / denom;
           }
         } // else Wlk > 1.0       
      }// else: Wlj != 0.0
      
      //cout << "Flag E:" << E << endl;
      //cout << "Flag W00:" << W00 << endl;
      //cout << "Flag W10:" << W10 << endl;
      //cout << "----------------" << endl;

      double tdc = (0.5 / dt) * ( acos(W00) * (A + B)  + asin(W10) * (C + D) + E);
      //cout << "Flag NPI:" << tdc << endl;
      nac.set(j, i, tdc);
      nac.set(i, j,-tdc);
    
    }// for j  
  }// for i

  return nac;

}


}// namespace libcalculators

}// liblibra


