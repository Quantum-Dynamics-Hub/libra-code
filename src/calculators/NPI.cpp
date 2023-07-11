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

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;
using namespace libmeigen;


/// libcalculators namespace
namespace libcalculators{


MATRIX nac_npi(MATRIX& St, double dt){
/** 
 Computes derivative coupling matrix elements using NPI of Meek and Levine

 Meek, G.A; Levine, B. G. "Evaluation of the Time-Derivative Coupling for Accurate
 Electronic State Transition Probabilities from Numerical Simulations" J. Phys. Chem. Lett.
 5, 2351-2356, 2014

*/

  int nstates = St.n_cols;
  MATRIX nac(nstates, nstates);

  double tol = 0.0;

  for(int i=0; i<nstates; i++){
    for(int j=i+1; j<nstates; j++){

// W_jk = <j|d/dt|k>      j->i;   k->j
// d_kj = ...

      double W00 = St.get(i,i); 
      double W01 = St.get(i,j);
      double W10 = St.get(j,i);
      double W11 = St.get(j,j);

      //if( (W00>1.0) && (W00<1.01)) { W00 = 0.99999; }
      //if( (W00>-1.01) && (W00<-1.0)) { W00 = -0.99999; }
      //if( (W11>1.0) && (W11<1.01)) { W11 = 0.99999; }
      //if( (W11>-1.01) && (W11<-1.0)) { W11 = -0.99999; }

      if( (W00>1.0) && (W00<1.01)) { W00 = 1.0; }
      if( (W00>-1.01) && (W00<-1.0)) { W00 = -1.0; }
      if( (W11>1.0) && (W11<1.01)) { W11 = 1.0; }
      if( (W11>-1.01) && (W11<-1.0)) { W11 = -1.0; }

      double A = acos(W00) - asin(W01);
      double B = acos(W00) + asin(W01);
      double C = acos(W11) - asin(W10);
      double D = acos(W11) + asin(W10);

      // Not sure what is this
      //  if (Wlj != Wlj){
      //      Wlj = 0.0;}

      if( fabs(A) <= tol){ A = -1.0; }
      else{  A = -1.0 * sin(A) / A; }

      if(fabs(B) <= tol){  B = 1.0; }
      else{   B = sin(B) / B; }

      if(fabs(C) <= tol){  C = 1.0; }
      else{   C = sin(C) / C; }

      if(fabs(D) <= tol){  D = 1.0; }
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
         double Wlk = -1.0 * (W01 * W00 + W11 * W10) / Wlj;
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
           E = 2.0 * asin(Wlj) * (Wlj * Wlk * sWlj + (E - 1.0) * sWlk) / denom;
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


