/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include "mEigen.h"
#include <cmath>

/// liblibra namespace
namespace liblibra{


using namespace Eigen;
using namespace std;
using namespace liblinalg;

/// libmeigen namespace
namespace libmeigen{





void sqrt_matrix(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half, double thresh, int do_phase_correction){
/**
  This function computes S^{1/2} and S^{-1/2} for given matrix S
  \param[in] S Input matrix
  \param[out] S_half Computed S^{1/2} matrix
  \param[out] S_i_half Computed S^{-1/2} matrix
  \param[in] threshold - if an absolute value of any eigenvalue of S is below this level, we stop,
   throwing an error message


*/

  if(S.n_cols != S.n_rows){
    cout<<"Error in libmeigen::sqrt_matrix : the input matrix is not square\n"; exit(0); 
  }
  if(S_half.n_cols != S_half.n_rows){
    cout<<"Error in libmeigen::sqrt_matrix : the output S^{1/2} matrix is not square\n"; exit(0); 
  }
  if(S_i_half.n_cols != S_i_half.n_rows){
    cout<<"Error in libmeigen::sqrt_matrix : the output S^{-1/2} matrix is not square\n"; exit(0); 
  }
  if(S.n_cols != S_half.n_cols){
    cout<<"Error in libmeigen::sqrt_matrix : size of matrix S is not the same as that of matrix S^{1/2}\n"; exit(0); 
  }
  if(S.n_cols != S_i_half.n_cols){
    cout<<"Error in libmeigen::sqrt_matrix : size of matrix S is not the same as that of matrix S^{-1/2}\n"; exit(0); 
  }


  int i,j;
 
  // Let us first diagonalize the overlap matrix S
  int sz = S.n_cols;  
  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0, 0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);


  // Find the eigenvalues of the the S matrix
  solve_eigen(S, *Seig, *C, 0);  // S * C = C * Seig  ==>  S = C * Seig * C.H()

  //if(do_phase_correction){   correct_phase(C);  }

  // Diagonal form of the S^{-1/2} and S^{1/2} matrices
  S_i_half = complex<double>(0.0,0.0);  // S^{-1/2}
  S_half = complex<double>(0.0,0.0);    // S^{1/2}

  for(i=0;i<sz;i++){
    complex<double> val = std::sqrt(Seig->get(i,i));

    double nrm  = std::abs(val);
    if(nrm<thresh){  
      std::cout<<"\n Error in sqrt_matrix: One of the eigenvalues of the matrix S is "<< val
               <<"\n this is below the used threshold of "<<thresh
               <<"\n So... the matrix is likely singular or your threshold is too large"
               <<"\n Exiting now...\n";
      exit(0);
    }
    else{
      S_i_half.M[i*sz+i] = 1.0/val;
      S_half.M[i*sz+i] = val;
    }
  }

  // Convert to the original basis
  S_i_half = (*C) * S_i_half * ((*C).H());
  S_half = (*C) * S_half * ((*C).H());

  delete C;
  delete Seig;


}// sqrt_matrix


void sqrt_matrix(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half, double thresh){

  sqrt_matrix(S, S_half, S_i_half, thresh, 0);

}

void sqrt_matrix(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half){

  sqrt_matrix(S, S_half, S_i_half, -1.0, 0);

}



void exp_matrix(CMATRIX& res, CMATRIX& S, complex<double> dt, int do_phase_correction){
/**
  This function computes exp(S*dt) for a given matrix S
  \param[in] S input matrix
  \param[in] dt scaling factor

*/

  if(S.n_cols != S.n_rows){
    cout<<"Error in libmeigen::exp_matrix : the input matrix is not square\n"; exit(0); 
  }


  int i,j;
 
  // Let us first diagonalize the input matrix x
  int sz = S.n_cols;  
  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0, 0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);


  // Find the eigenvalues of the the S matrix
  solve_eigen(S, *Seig, *C, 0);  // S * C = C * Seig  ==>  S = C * Seig * C.H()

  //if(do_phase_correction){   correct_phase(C);  }

  
  for(i=0;i<sz;i++){ res.M[i*sz+i]= std::exp(dt * Seig->get(i,i)); }

  // Convert to the original basis
  res = (*C) * res * ((*C).H());

  delete C;
  delete Seig;

}// exp_matrix


void exp_matrix(CMATRIX& res, CMATRIX& S, complex<double> dt){
  exp_matrix(res,S, dt, 0);
}



}// namespace libmeigen
}// namespace liblibra
