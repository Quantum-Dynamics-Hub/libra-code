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

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <Eigen/SVD>
#endif 

#include "mEigen.h"

/// liblibra namespace
namespace liblibra{

using namespace Eigen;
using namespace std;
using namespace liblinalg;

/// libmeigen namespace
namespace libmeigen{



void FullPivLU_decomposition(MATRIX& A, MATRIX& P, MATRIX& L, MATRIX& U, MATRIX& Q){
/** A wrapper of Eigen::FullPivLU<MatrixXd>.matrixLU() 
   A - the source matrix
   P - the permutation matrix
   L - lower triangular
   U - upper triangular
   Q - another permutation matrix

  Really, what happens is:

  A = P * L * U * Q,  
  P, Q - are the permutation matrices

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


  MatrixXd l(N,N); 
  l = lu.matrixLU().triangularView<StrictlyLower>();

  MatrixXd u(N,N);
  u = lu.matrixLU().triangularView<Upper>();

  MatrixXd p(N,N);
  p = lu.permutationP(); 

  MatrixXd q(N,N);
  q = lu.permutationQ(); 


  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      P.M[i*N+j] = p(i,j);

      L.M[i*N+j] = l(i,j);
      if(i==j){   L.M[i*N+j] = 1.0;  }

      U.M[i*N+j] = u(i,j);

      Q.M[i*N+j] = q(i,j);

    }// for j
  }// for i
  
}


void FullPivLU_decomposition(CMATRIX& A, CMATRIX& P, CMATRIX& L, CMATRIX& U, CMATRIX& Q){
/** A wrapper of Eigen::FullPivLU<MatrixXcd>.matrixLU() 
   A - the source matrix
   P - the permutation matrix
   L - lower triangular
   U - upper triangular
   Q - another permutation matrix

  Really, what happens is:

  A = P * L * U * Q,  
  P, Q - are the permutation matrices

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


  MatrixXcd l(N,N); 
  l = lu.matrixLU().triangularView<StrictlyLower>();

  MatrixXcd u(N,N);
  u = lu.matrixLU().triangularView<Upper>();

  MatrixXcd p(N,N);
  p = lu.permutationP(); 

  MatrixXcd q(N,N);
  q = lu.permutationQ(); 


  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      P.M[i*N+j] = p(i,j);

      L.M[i*N+j] = l(i,j);
      if(i==j){   L.M[i*N+j] = 1.0;  }

      U.M[i*N+j] = u(i,j);

      Q.M[i*N+j] = q(i,j);

    }// for j
  }// for i
  
}




void JacobiSVD_decomposition(CMATRIX& A, CMATRIX& U, CMATRIX& S, CMATRIX& V){
/** 

   A wrapper of Eigen::JacobiSVD<MatrixXd> class

   A - the source matrix:  A = U * S * V^{*}
   U - left matrix
   S - singular values matrix
   V - right matrix

*/

  

  int N = A.n_cols;
  int i,j;

  MatrixXcd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i


  Eigen::JacobiSVD<Eigen::MatrixXcd> svd(a, Eigen::ComputeFullU | Eigen::ComputeFullV);

  MatrixXcd u(N,N); 
  MatrixXcd s(N,1); 
  MatrixXcd v(N,N); 

  u = svd.matrixU();
  s = svd.singularValues();
  v = svd.matrixV();

  for(i=0;i<N;i++){

    for(j=0;j<N;j++){

      U.M[i*N+j] = u(i,j);
      V.M[i*N+j] = v(i,j);

      if(j==i){    S.M[i*N+i] = s(i,0); }
      else{ S.M[i*N+j] = complex<double>(0.0, 0.0); }


    }// for j
  }// for i
  
}


void BDCSVD_decomposition(CMATRIX& A, CMATRIX& U, CMATRIX& S, CMATRIX& V){
/** 

   A wrapper of Eigen::BDCSVD<MatrixXd> class

   A - the source matrix:  A = U * S * V^{*}
   U - left matrix
   S - singular values matrix
   V - right matrix

*/

  

  int N = A.n_cols;
  int i,j;

  MatrixXcd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i


  Eigen::BDCSVD<Eigen::MatrixXcd> svd(a, Eigen::ComputeFullU | Eigen::ComputeFullV);

  MatrixXcd u(N,N); 
  MatrixXcd s(N,1); 
  MatrixXcd v(N,N); 

  u = svd.matrixU();
  s = svd.singularValues();
  v = svd.matrixV();

  for(i=0;i<N;i++){

    for(j=0;j<N;j++){

      U.M[i*N+j] = u(i,j);
      V.M[i*N+j] = v(i,j);

      if(j==i){    S.M[i*N+i] = s(i,0); }
      else{ S.M[i*N+j] = complex<double>(0.0, 0.0); }

    }// for j
  }// for i
  
}




}// namespace libmeigen
}// namespace liblibra
