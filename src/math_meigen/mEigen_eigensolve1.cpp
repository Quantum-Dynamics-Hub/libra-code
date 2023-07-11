/*********************************************************************************
* Copyright (C) 2015-2019 Alexey V. Akimov
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
#endif 

#include "mEigen.h"

/// liblibra namespace
namespace liblibra{

using namespace Eigen;
using namespace std;
using namespace liblinalg;

/// libmeigen namespace
namespace libmeigen{



void solve_eigen(CMATRIX* H, CMATRIX* S, CMATRIX* E, CMATRIX* C, int symm){
/** Solve the generalized eigenvalue problem: H * C = S * C * E
  i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
  C[i][j] - is the weight of j-th AO in i-th MO
  More, generally C.col(i) is the eigenvector corresponding the eigenstate E_i

  int symm - is the flag for H and S matrices symmetrization:
   symm = 0 - don't symmetrize (use the input matrices as they are)  
   symm = 1 - do symmetrize them

  H * C = (N_bas x N_bas)   x  ( N_bas x N_mo ) = N_mo x N_mo

  S * C * E = (N_bas x N_bas) x (N_bas x N_mo ) x (N_mo x N_mo) = N_bas x N_mo
 

*/


  int i,j;

  if(H->n_rows!=S->n_rows){  // N_bas
    std::cout<<"Error in solve_eigen: The # of rows of the H and S matrices must be equal\n";
    exit(0);
  }
  if(H->n_cols!=C->n_rows){  // N_bas
    std::cout<<"Error in solve_eigen: The # of cols of H and the # of rows of C matrices must be equal\n";
    std::cout<<"Make sure C is allocated \n";
    exit(0);
  }
  if(C->n_cols!=E->n_rows){  // N_mo
    std::cout<<"Error in solve_eigen: The # of cols of C and the # of rows of E matrices must be equal\n";
    std::cout<<"Make sure C and E are allocated \n";
    exit(0);
  }
  if(C->n_cols!=E->n_cols){  // N_mo
    std::cout<<"Error in solve_eigen: The # of cols of C and the # of cols of E matrices must be equal\n";
    std::cout<<"Make sure C and E are allocated \n";
    exit(0);
  }
  // The above two conditions imply E is square matrix
  int N_bas = H->n_cols;
  int N_mo  = C->n_cols;


  // Wrapper matrices for Eigen3
  MatrixXcd A(N_bas,N_bas), B(N_bas,N_bas);

  if(symm==0){
    for(i=0;i<N_bas;i++){
      for(j=0;j<N_bas;j++){
        // Symmetrize, to reduce numerical errors:
        A(i,j) = H->M[i*N_bas+j];     
        B(i,j) = S->M[i*N_bas+j];     
      }// for j
    }// for i
  }// no symmetrization

  else if(symm==1){
    for(i=0;i<N_bas;i++){
      for(j=0;j<N_bas;j++){
        // Symmetrize, to reduce numerical errors:
        A(i,j) = 0.5*(H->M[i*N_bas+j] + std::conj(H->M[j*N_bas+i]));     
        B(i,j) = 0.5*(S->M[i*N_bas+j] + std::conj(S->M[j*N_bas+i]));     
      }// for j
    }// for i
  } // symmetrize


  // Solve eigenvalue problem
  GeneralizedSelfAdjointEigenSolver<MatrixXcd> solution(A,B);
  if(solution.info()!=Success){ cout<<"Eigen fails\n";   exit(0); }


  // Copy results into output matrices
  for(i=0;i<N_mo;i++){

     E->M[i*N_mo+i] = solution.eigenvalues()[i]; 

     for(j=0;j<N_bas;j++){  C->M[j*N_mo+i] = solution.eigenvectors().col(i)[j];    }// for j

  }// for i

  if(0){
    // Checking solution - for debug purposes, inactive by default
    cout<<"in solve_eigen (with MATRIX):\n";
    cout<<"A = "<<A<<endl;
    cout<<"B = "<<B<<endl;
    cout<<"E = "<<(*E)<<endl;
    cout<<"C = "<<(*C)<<endl;
    cout<<"S = "<<(*S)<<endl;
    cout<<"H = "<<(*H)<<endl;
    cout<<"A*C = "<<(*H) * (*C)<<endl;
    cout<<"S*C*E = "<<(*S) * (*C) * (*E)<<endl;
    cout<<"A*C - S*C*E = "<<*H * (*C) - (*S) * (*C) * (*E)<<endl;
    cout<<"C.H()*S*C= "<<((*C).H()) * (*S) * (*C)<<endl; // This gives unity, indeed.
  }


}//void solve_eigen(CMATRIX* H, CMATRIX* S, CMATRIX* E, CMATRIX* C, int symm)


void solve_eigen(CMATRIX& H, CMATRIX& S, CMATRIX& E, CMATRIX& C, int symm){
/** Solve the generalized eigenvalue problem: H * C = S * C * E
  i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
  C[i][j] - is the weight of j-th AO in i-th MO
  More, generally C.col(i) is the eigenvector corresponding the eigenstate E_i

  int symm - is the flag for H and S matrices symmetrization:
   symm = 0 - don't symmetrize (use the input matrices as they are)  
   symm = 1 - do symmetrize them

  H * C = (N_bas x N_bas)   x  ( N_bas x N_mo ) = N_mo x N_mo

  S * C * E = (N_bas x N_bas) x (N_bas x N_mo ) x (N_mo x N_mo) = N_bas x N_mo
*/ 

  solve_eigen(&H, &S, &E, &C, symm);

}



void solve_eigen(MATRIX* H, MATRIX* S, CMATRIX* E, CMATRIX* C, int symm){
/** Solve the generalized eigenvalue problem: H * C = S * C * E
  i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
  C[i][j] - is the weight of j-th AO in i-th MO
  More, generally C.col(i) is the eigenvector corresponding the eigenstate E_i

  int symm - is the flag for H and S matrices symmetrization:
   symm = 0 - don't symmetrize (use the input matrices as they are)  
   symm = 1 - do symmetrize them

  H * C = (N_bas x N_bas)   x  ( N_bas x N_mo ) = N_mo x N_mo

  S * C * E = (N_bas x N_bas) x (N_bas x N_mo ) x (N_mo x N_mo) = N_bas x N_mo
 

*/


  int i,j;

  if(H->n_rows!=S->n_rows){  // N_bas
    std::cout<<"Error in solve_eigen: The # of rows of the H and S matrices must be equal\n";
    exit(0);
  }
  if(H->n_cols!=C->n_rows){  // N_bas
    std::cout<<"Error in solve_eigen: The # of cols of H and the # of rows of C matrices must be equal\n";
    std::cout<<"Make sure C is allocated \n";
    exit(0);
  }
  if(C->n_cols!=E->n_rows){  // N_mo
    std::cout<<"Error in solve_eigen: The # of cols of C and the # of rows of E matrices must be equal\n";
    std::cout<<"Make sure C and E are allocated \n";
    exit(0);
  }
  if(C->n_cols!=E->n_cols){  // N_mo
    std::cout<<"Error in solve_eigen: The # of cols of C and the # of cols of E matrices must be equal\n";
    std::cout<<"Make sure C and E are allocated \n";
    exit(0);
  }
  // The above two conditions imply E is square matrix
  int N_bas = H->n_cols;
  int N_mo  = C->n_cols;


  *E = 0.0;
  *C = 0.0;

  // Wrapper matrices for Eigen3
  MatrixXcd A(N_bas,N_bas), B(N_bas,N_bas);

  if(symm==0){
    for(i=0;i<N_bas;i++){
      for(j=0;j<N_bas;j++){
        // Symmetrize, to reduce numerical errors:
        A(i,j) = complex<double>(H->M[i*N_bas+j], 0.0);     
        B(i,j) = complex<double>(S->M[i*N_bas+j], 0.0);     
      }// for j
    }// for i
  }// no symmetrization

  else if(symm==1){
    for(i=0;i<N_bas;i++){
      for(j=0;j<N_bas;j++){
        // Symmetrize, to reduce numerical errors:
        A(i,j) = complex<double>( 0.5*(H->M[i*N_bas+j] + H->M[j*N_bas+i]), 0.0);     
        B(i,j) = complex<double>( 0.5*(S->M[i*N_bas+j] + S->M[j*N_bas+i]), 0.0);     
      }// for j
    }// for i
  } // symmetrize


  // Solve eigenvalue problem
  GeneralizedSelfAdjointEigenSolver<MatrixXcd> solution(A,B);
  if(solution.info()!=Success){ cout<<"Eigen fails\n";   exit(0); }


  // Copy results into output matrices
  for(i=0;i<N_mo;i++){

     E->M[i*N_mo+i] = solution.eigenvalues()[i]; 

     for(j=0;j<N_bas;j++){  C->M[j*N_mo+i] = solution.eigenvectors().col(i)[j];    }// for j

  }// for i

  if(0){
    // Checking solution - for debug purposes, inactive by default
    cout<<"in solve_eigen (with MATRIX):\n";
    cout<<"A = "<<A<<endl;
    cout<<"B = "<<B<<endl;
    cout<<"E = "<<(*E)<<endl;
    cout<<"C = "<<(*C)<<endl;
    cout<<"S = "<<(*S)<<endl;
    cout<<"H = "<<(*H)<<endl;
//    cout<<"A*C = "<<(*H) * (*C)<<endl;
//    cout<<"S*C*E = "<<(*S) * (*C) * (*E)<<endl;
//    cout<<"A*C - S*C*E = "<<*H * (*C) - (*S) * (*C) * (*E)<<endl;
//    cout<<"C.H()*S*C= "<<((*C).H()) * (*S) * (*C)<<endl; // This gives unity, indeed.

  }


}//void solve_eigen(MATRIX* H, MATRIX* S, CMATRIX* E, CMATRIX* C, int symm)



void solve_eigen(MATRIX& H, MATRIX& S, CMATRIX& E, CMATRIX& C, int symm){
/** Solve the generalized eigenvalue problem: H * C = S * C * E
  i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
  C[i][j] - is the weight of j-th AO in i-th MO
  More, generally C.col(i) is the eigenvector corresponding the eigenstate E_i

  int symm - is the flag for H and S matrices symmetrization:
   symm = 0 - don't symmetrize (use the input matrices as they are)  
   symm = 1 - do symmetrize them

  H * C = (N_bas x N_bas)   x  ( N_bas x N_mo ) = N_mo x N_mo

  S * C * E = (N_bas x N_bas) x (N_bas x N_mo ) x (N_mo x N_mo) = N_bas x N_mo
 

*/

  solve_eigen(&H, &S, &E, &C, symm);

}



void solve_eigen(MATRIX* H, MATRIX* S, MATRIX* E, MATRIX* C, int symm){
/** Solve the generalized eigenvalue problem: H * C = S * C * E
  i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
  C[i][j] - is the weight of j-th AO in i-th MO
  More, generally C.col(i) is the eigenvector corresponding the eigenstate E_i

  int symm - is the flag for H and S matrices symmetrization:
   symm = 0 - don't symmetrize (use the input matrices as they are)  
   symm = 1 - do symmetrize them

  H * C = (N_bas x N_bas)   x  ( N_bas x N_mo ) = N_mo x N_mo

  S * C * E = (N_bas x N_bas) x (N_bas x N_mo ) x (N_mo x N_mo) = N_bas x N_mo
 

*/


  int i,j;

  if(H->n_rows!=S->n_rows){  // N_bas
    std::cout<<"Error in solve_eigen: The # of rows of the H and S matrices must be equal\n";
    exit(0);
  }
  if(H->n_cols!=C->n_rows){  // N_bas
    std::cout<<"Error in solve_eigen: The # of cols of H and the # of rows of C matrices must be equal\n";
    std::cout<<"Make sure C is allocated \n";
    exit(0);
  }
  if(C->n_cols!=E->n_rows){  // N_mo
    std::cout<<"Error in solve_eigen: The # of cols of C and the # of rows of E matrices must be equal\n";
    std::cout<<"Make sure C and E are allocated \n";
    exit(0);
  }
  if(C->n_cols!=E->n_cols){  // N_mo
    std::cout<<"Error in solve_eigen: The # of cols of C and the # of cols of E matrices must be equal\n";
    std::cout<<"Make sure C and E are allocated \n";
    exit(0);
  }
  // The above two conditions imply E is square matrix
  int N_bas = H->n_cols;
  int N_mo  = C->n_cols;


  *E = 0.0;
  *C = 0.0;

  // Wrapper matrices for Eigen3
  MatrixXd A(N_bas,N_bas), B(N_bas,N_bas);

  if(symm==0){
    for(i=0;i<N_bas;i++){
      for(j=0;j<N_bas;j++){
        // Symmetrize, to reduce numerical errors:
        A(i,j) = H->M[i*N_bas+j];     
        B(i,j) = S->M[i*N_bas+j];     
      }// for j
    }// for i
  }// no symmetrization

  else if(symm==1){
    for(i=0;i<N_bas;i++){
      for(j=0;j<N_bas;j++){
        // Symmetrize, to reduce numerical errors:
        A(i,j) = 0.5*(H->M[i*N_bas+j] + H->M[j*N_bas+i]);     
        B(i,j) = 0.5*(S->M[i*N_bas+j] + S->M[j*N_bas+i]);     
      }// for j
    }// for i
  } // symmetrize


  // Solve eigenvalue problem
  GeneralizedSelfAdjointEigenSolver<MatrixXd> solution(A,B);
  if(solution.info()!=Success){ cout<<"Eigen fails\n";   exit(0); }


  // Copy results into output matrices
  for(i=0;i<N_mo;i++){

     E->M[i*N_mo+i] = solution.eigenvalues()[i]; 

     for(j=0;j<N_bas;j++){  C->M[j*N_mo+i] = solution.eigenvectors().col(i)[j];    }// for j

  }// for i

  if(0){
    // Checking solution - for debug purposes, inactive by default
    cout<<"in solve_eigen (with MATRIX):\n";
    cout<<"A = "<<A<<endl;
    cout<<"B = "<<B<<endl;
    cout<<"E = "<<(*E)<<endl;
    cout<<"C = "<<(*C)<<endl;
    cout<<"S = "<<(*S)<<endl;
    cout<<"H = "<<(*H)<<endl;
    cout<<"A*C = "<<(*H) * (*C)<<endl;
    cout<<"S*C*E = "<<(*S) * (*C) * (*E)<<endl;
    cout<<"A*C - S*C*E = "<<*H * (*C) - (*S) * (*C) * (*E)<<endl;
    cout<<"C.T()*S*C= "<<((*C).T()) * (*S) * (*C)<<endl; // This gives unity, indeed.
  }


}//void solve_eigen(MATRIX* H, MATRIX* S, MATRIX* E, MATRIX* C, int symm)



void solve_eigen(MATRIX& H, MATRIX& S, MATRIX& E, MATRIX& C, int symm){
/** Solve the generalized eigenvalue problem: H * C = S * C * E
  i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
  C[i][j] - is the weight of j-th AO in i-th MO
  More, generally C.col(i) is the eigenvector corresponding the eigenstate E_i

  int symm - is the flag for H and S matrices symmetrization:
   symm = 0 - don't symmetrize (use the input matrices as they are)  
   symm = 1 - do symmetrize them

  H * C = (N_bas x N_bas)   x  ( N_bas x N_mo ) = N_mo x N_mo

  S * C * E = (N_bas x N_bas) x (N_bas x N_mo ) x (N_mo x N_mo) = N_bas x N_mo
 

*/

  solve_eigen(&H, &S, &E, &C, symm);

}







}// namespace libmeigen
}// namespace liblibra
