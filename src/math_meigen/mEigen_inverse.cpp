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


/// liblibra namespace
namespace liblibra{

using namespace Eigen;
using namespace std;
using namespace liblinalg;

/// libmeigen namespace
namespace libmeigen{



void FullPivLU_rank_invertible(MATRIX& A, int& rank, int& is_inver){
/** Find whether the matrix is invertible. This result is stored in the
  passed variable is_inver. The function will also let you know the
  rank of the provided matrix (will be stored in the "rank" variable)   
  
*/

  int N = A.n_cols;
  int i,j;

  MatrixXd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i

  Eigen::FullPivLU<MatrixXd> lu(a);

  rank = lu.rank();

  is_inver = (int)lu.isInvertible();

  
}


void FullPivLU_rank_invertible(CMATRIX& A, int& rank, int& is_inver){
/** Find whether the matrix is invertible. This result is stored in the
  passed variable is_inver. The function will also let you know the
  rank of the provided matrix (will be stored in the "rank" variable)   
  
*/


  int N = A.n_cols;
  int i,j;

  MatrixXcd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i

  Eigen::FullPivLU<MatrixXcd> lu(a);

  rank = lu.rank();

  is_inver = (int)lu.isInvertible();

  
}

boost::python::list FullPivLU_rank_invertible(MATRIX& A){

  int rank = -1;
  int is_inver = -1;

  FullPivLU_rank_invertible(A, rank, is_inver);

  boost::python::list res;
  res.append(rank);
  res.append(is_inver);

  return res;
}

boost::python::list FullPivLU_rank_invertible(CMATRIX& A){

  int rank = -1;
  int is_inver = -1;

  FullPivLU_rank_invertible(A, rank, is_inver);

  boost::python::list res;
  res.append(rank);
  res.append(is_inver);

  return res;
}





void FullPivLU_inverse(MATRIX& A, MATRIX& invA){
/** A wrapper of the Eigen::FullPivLU<MatrixXcd>.inverse() method */

  int N = A.n_cols;
  int i,j;

  MatrixXd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i

  Eigen::FullPivLU<MatrixXd> lu(a);

  MatrixXd inva; inva = lu.inverse();

  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      invA.M[i*N+j] = inva(i,j);
    }// for j
  }// for i


}


void FullPivLU_inverse(CMATRIX& A, CMATRIX& invA){
/** A wrapper of the Eigen::FullPivLU<MatrixXcd>.inverse() method */

  int N = A.n_cols;
  int i,j;

  MatrixXcd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i

  Eigen::FullPivLU<MatrixXcd> lu(a);

  MatrixXcd inva; inva = lu.inverse();

  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      invA.M[i*N+j] = inva(i,j);
    }// for j
  }// for i


}





void inv_matrix(MATRIX& S, MATRIX& S_inv, double thresh, int do_phase_correction){
/**
  This function computes S^{-1} of a given matrix S
  \param[in] S Input matrix
  \param[out] S_inv Computed S^{-1} matrix
  \param[in] threshold - if an absolute value of any eigenvalue of S is below this level, we stop,
   throwing an error message

*/

  if(S.n_cols != S.n_rows){
    cout<<"Error in libmeigen::inv_matrix : the input matrix is not square\n"; exit(0); 
  }
  if(S_inv.n_cols != S_inv.n_rows){
    cout<<"Error in libmeigen::inv_matrix : the output S^{-1} matrix is not square\n"; exit(0); 
  }
  if(S.n_cols != S_inv.n_cols){
    cout<<"Error in libmeigen::inv_matrix : size of matrix S is not the same as that of matrix S^{-1}\n"; exit(0); 
  }


  int i,j;
 
  // Let us first diagonalize the overlap matrix S
  int sz = S.n_cols;  

  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0,0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);
  CMATRIX* tmp; tmp = new CMATRIX(sz, sz);  *tmp = complex<double>(0.0,0.0);

  // Find the eigenvalues of the the S matrix
  libmeigen::solve_eigen(S, *Seig, *C, 0);  // S * C = C * Seig  ==>  S = C * Seig * C.H()

  //if(do_phase_correction){   correct_phase(C);  }


  for(i=0;i<sz;i++){
    complex<double> val = Seig->get(i,i);
    double nrm = std::abs(val);
    if(nrm<thresh){  
      std::cout<<"\n Error in inv_matrix: One of the eigenvalues of the matrix S is "<< val
               <<"\n this is below the used threshold of "<<thresh
               <<"\n So... the matrix is likely singular or your threshold is too large"
               <<"\n Exiting now...\n";
      exit(0);
    }
    else{    tmp->M[i*sz+i] = 1.0/val; }
  }

  // Convert to the original basis
  S_inv = ((*C) * (*tmp) * ((*C).H())).real();

  delete C;
  delete Seig;
  delete tmp;


}// inv_matrix


void inv_matrix(MATRIX& S, MATRIX& S_inv, double thresh){
/**
  This function computes S^{-1} of a given matrix S
  \param[in] S Input matrix
  \param[out] S_inv Computed S^{-1} matrix
  \param[in] threshold - if an absolute value of any eigenvalue of S is below this level, we stop,
   throwing an error message

*/
  inv_matrix(S, S_inv, thresh, 1);

}

void inv_matrix(MATRIX& S, MATRIX& S_inv){

  inv_matrix(S, S_inv, -1.0, 1);

}



void inv_matrix(CMATRIX& S, CMATRIX& S_inv, double thresh, int do_phase_correction){
/**
  This function computes S^{-1} of a given matrix S
  \param[in] S Input matrix
  \param[out] S_inv Computed S^{-1} matrix
  \param[in] threshold - if an absolute value of any eigenvalue of S is below this level, we stop,
   throwing an error message

*/

  if(S.n_cols != S.n_rows){
    cout<<"Error in libmeigen::inv_matrix : the input matrix is not square\n"; exit(0); 
  }
  if(S_inv.n_cols != S_inv.n_rows){
    cout<<"Error in libmeigen::inv_matrix : the output S^{-1} matrix is not square\n"; exit(0); 
  }
  if(S.n_cols != S_inv.n_cols){
    cout<<"Error in libmeigen::inv_matrix : size of matrix S is not the same as that of matrix S^{-1}\n"; exit(0); 
  }


  int i,j;
 
  // Let us first diagonalize the overlap matrix S
  int sz = S.n_cols;  

  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0,0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);

  // Find the eigenvalues of the the S matrix
  libmeigen::solve_eigen(S, *Seig, *C, 0);  // S * C = C * Seig  ==>  S = C * Seig * C.H()

  //if(do_phase_correction){   correct_phase(C);  }

  // Diagonal form of the S^{-1} matrix
  S_inv = complex<double>(0.0,0.0);  // S^{-1}


  for(i=0;i<sz;i++){
    complex<double> val = Seig->get(i,i);
    double nrm = std::abs(val);
    if(nrm<thresh){  
      std::cout<<"\n Error in inv_matrix: One of the eigenvalues of the matrix S is "<< val
               <<"\n this is below the used threshold of "<<thresh
               <<"\n So... the matrix is likely singular or your threshold is too large"
               <<"\n Exiting now...\n";
      exit(0);
    }
    else{    S_inv.M[i*sz+i] = 1.0/val; }
  }

  // Convert to the original basis
  S_inv = (*C) * S_inv * ((*C).H());

  delete C;
  delete Seig;


}// inv_matrix

void inv_matrix(CMATRIX& S, CMATRIX& S_inv, double thresh){

  inv_matrix(S, S_inv, thresh, 1);

}

void inv_matrix(CMATRIX& S, CMATRIX& S_inv){

  inv_matrix(S, S_inv, -1.0, 1);

}



}// namespace libmeigen
}// namespace liblibra
