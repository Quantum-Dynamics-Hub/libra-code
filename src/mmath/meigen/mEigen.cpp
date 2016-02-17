/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
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

using namespace Eigen;
using namespace std;

/// libmmath namespace
namespace libmmath{

using namespace liblinalg;

/// libmeigen namespace
namespace libmeigen{

/****************************************************************************
  This file contains following functions:

  void solve_eigen(int Norb, MATRIX* H, MATRIX* S, MATRIX* E, MATRIX* C)
  void solve_eigen(int Norb, MATRIX* H, MATRIX* E, MATRIX* C)
  int merge_sort(vector< pair<int,double> >& in, vector< pair<int,double> >& out)

****************************************************************************/

void solve_eigen(int Norb, MATRIX* H, MATRIX* S, MATRIX* E, MATRIX* C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  *E = 0.0;
  *C = 0.0;

  int i,j;

  // Wrapper matrices for Eigen3
  MatrixXd A(Norb,Norb), B(Norb,Norb);
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      // Symmetrize, to reduce numerical errors:
      A(i,j) = 0.5*(H->M[i*Norb+j] + H->M[j*Norb+i]);     
      B(i,j) = 0.5*(S->M[i*Norb+j] + S->M[j*Norb+i]);     
    }// for j
  }// for i


  // Solve eigenvalue problem
  GeneralizedSelfAdjointEigenSolver<MatrixXd> solution(A,B);
  if(solution.info()!=Success){ cout<<"Eigen fails\n";   exit(0); }


  // Copy results into output matrices
  for(i=0;i<Norb;i++){

     E->M[i*Norb+i] = solution.eigenvalues()[i]; 

    for(j=0;j<Norb;j++){  C->M[j*Norb+i] = solution.eigenvectors().col(i)[j];    }// for j

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


}//void solve_eigen(int Norb, MATRIX* H, MATRIX* S, MATRIX* E, MATRIX* C)




void solve_eigen(int Norb, MATRIX& H, MATRIX& S, MATRIX& E, MATRIX& C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  solve_eigen(Norb, &H, &S, &E, &C);

}


void solve_eigen_gen(int Norb, MATRIX* H, MATRIX* S, MATRIX* E, MATRIX* C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  *E = 0.0;
  *C = 0.0;

  int i,j;

  // Wrapper matrices for Eigen3
  MatrixXd A(Norb,Norb), B(Norb,Norb);
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      // Symmetrize, to reduce numerical errors:
      A(i,j) = H->M[i*Norb+j];     
      B(i,j) = S->M[i*Norb+j];     
    }// for j
  }// for i


  // Solve eigenvalue problem
  GeneralizedSelfAdjointEigenSolver<MatrixXd> solution(A,B);
  if(solution.info()!=Success){ cout<<"Eigen fails\n";   exit(0); }


  // Copy results into output matrices
  for(i=0;i<Norb;i++){

     E->M[i*Norb+i] = solution.eigenvalues()[i]; 

    for(j=0;j<Norb;j++){  C->M[j*Norb+i] = solution.eigenvectors().col(i)[j];    }// for j

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


}//void solve_eigen(int Norb, MATRIX* H, MATRIX* S, MATRIX* E, MATRIX* C)




void solve_eigen_gen(int Norb, MATRIX& H, MATRIX& S, MATRIX& E, MATRIX& C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  solve_eigen_gen(Norb, &H, &S, &E, &C);

}










void solve_eigen(int Norb, MATRIX* H, MATRIX* S, CMATRIX* E, CMATRIX* C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  *E = 0.0;
  *C = 0.0;

  int i,j;

  // Wrapper matrices for Eigen3
  MatrixXcd A(Norb,Norb), B(Norb,Norb);
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      // Symmetrize, to reduce numerical errors:
      A(i,j) = complex<double>(0.5*(H->M[i*Norb+j] + H->M[j*Norb+i]), 0.0);     
      B(i,j) = complex<double>(0.5*(S->M[i*Norb+j] + S->M[j*Norb+i]), 0.0);     
    }// for j
  }// for i


  // Solve eigenvalue problem
  GeneralizedSelfAdjointEigenSolver<MatrixXcd> solution(A,B);
  if(solution.info()!=Success){ cout<<"Eigen fails\n";   exit(0); }


  // Copy results into output matrices
  for(i=0;i<Norb;i++){

     E->M[i*Norb+i] = solution.eigenvalues()[i]; 

    for(j=0;j<Norb;j++){  C->M[j*Norb+i] = solution.eigenvectors().col(i)[j];    }// for j

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


}//void solve_eigen(int Norb, MATRIX* H, MATRIX* S, CMATRIX* E, CMATRIX* C)


void solve_eigen(int Norb, MATRIX& H, MATRIX& S, CMATRIX& E, CMATRIX& C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  solve_eigen(Norb, &H, &S, &E, &C);

}


void solve_eigen_gen(int Norb, MATRIX* H, MATRIX* S, CMATRIX* E, CMATRIX* C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  *E = 0.0;
  *C = 0.0;

  int i,j;

  // Wrapper matrices for Eigen3
  MatrixXcd A(Norb,Norb), B(Norb,Norb);
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      // Symmetrize, to reduce numerical errors:
      A(i,j) = complex<double>(H->M[i*Norb+j], 0.0);     
      B(i,j) = complex<double>(S->M[i*Norb+j], 0.0);     
    }// for j
  }// for i


  // Solve eigenvalue problem
  GeneralizedSelfAdjointEigenSolver<MatrixXcd> solution(A,B);
  if(solution.info()!=Success){ cout<<"Eigen fails\n";   exit(0); }


  // Copy results into output matrices
  for(i=0;i<Norb;i++){

     E->M[i*Norb+i] = solution.eigenvalues()[i]; 

    for(j=0;j<Norb;j++){  C->M[j*Norb+i] = solution.eigenvectors().col(i)[j];    }// for j

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


}//void solve_eigen(int Norb, MATRIX* H, MATRIX* S, CMATRIX* E, CMATRIX* C)


void solve_eigen_gen(int Norb, MATRIX& H, MATRIX& S, CMATRIX& E, CMATRIX& C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  solve_eigen_gen(Norb, &H, &S, &E, &C);

}






void solve_eigen(int Norb, CMATRIX* H, CMATRIX* S, CMATRIX* E, CMATRIX* C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  *E = 0.0;
  *C = 0.0;

  int i,j;

  // Wrapper matrices for Eigen3
  MatrixXcd A(Norb,Norb), B(Norb,Norb);
  
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      // Symmetrize, to reduce numerical errors:
      A(i,j) = 0.5*(H->M[i*Norb+j] + std::conj(H->M[j*Norb+i]));     
      B(i,j) = 0.5*(S->M[i*Norb+j] + std::conj(S->M[j*Norb+i]));     
    }// for j
  }// for i


  // Solve eigenvalue problem
  GeneralizedSelfAdjointEigenSolver<MatrixXcd> solution(A,B);
  if(solution.info()!=Success){ cout<<"Eigen fails\n";   exit(0); }


  // Copy results into output matrices
  for(i=0;i<Norb;i++){

     E->M[i*Norb+i] = solution.eigenvalues()[i]; 

    for(j=0;j<Norb;j++){  C->M[j*Norb+i] = solution.eigenvectors().col(i)[j];    }// for j

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


}//void solve_eigen(int Norb, MATRIX* H, MATRIX* S, CMATRIX* E, CMATRIX* C)


void solve_eigen(int Norb, CMATRIX& H, CMATRIX& S, CMATRIX& E, CMATRIX& C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  solve_eigen(Norb, &H, &S, &E, &C);

}


void solve_eigen_gen(int Norb, CMATRIX* H, CMATRIX* S, CMATRIX* E, CMATRIX* C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO
// No additional symmetrzation!!! - good for non-Hermitian Hamiltonians

  *E = 0.0;
  *C = 0.0;

  int i,j;

  // Wrapper matrices for Eigen3
  MatrixXcd A(Norb,Norb), B(Norb,Norb);
  
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      // Symmetrize, to reduce numerical errors:
      A(i,j) = H->M[i*Norb+j];     
      B(i,j) = S->M[i*Norb+j];     
    }// for j
  }// for i


  // Solve eigenvalue problem
  GeneralizedSelfAdjointEigenSolver<MatrixXcd> solution(A,B);
  //GeneralizedEigenSolver<MatrixXcd> solution(A,B);
  if(solution.info()!=Success){ cout<<"Eigen fails\n";   exit(0); }


  // Copy results into output matrices
  for(i=0;i<Norb;i++){

     E->M[i*Norb+i] = solution.eigenvalues()[i]; 

    for(j=0;j<Norb;j++){  C->M[j*Norb+i] = solution.eigenvectors().col(i)[j];    }// for j

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


}//void solve_eigen(int Norb, MATRIX* H, MATRIX* S, CMATRIX* E, CMATRIX* C)


void solve_eigen_gen(int Norb, CMATRIX& H, CMATRIX& S, CMATRIX& E, CMATRIX& C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  solve_eigen_gen(Norb, &H, &S, &E, &C);

}







void solve_eigen(int Norb, MATRIX* H, MATRIX* E, MATRIX* C){
// Solve H * C = C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// ACHTUNG:    EigenSolver<MatrixXd> solution(A); - messes up the ordering of eigenvalues
// GeneralizedSelfAdjointEigenSolver<MatrixXd> solution(A,B); does not!, so 
// we eventually would need to abandon the former (special case) solver

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  *E = 0.0;
  *C = 0.0;

  int i,j;

  // Wrapper matrices for Eigen3
  MatrixXd A(Norb,Norb);
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      // Symmetrize, to reduce numerical errors:
      A(i,j) = 0.5*(H->M[i*Norb+j] + H->M[j*Norb+i]);     
    }// for j
  }// for i


  // Solve eigenvalue problem
  EigenSolver<MatrixXd> solution(A);


  // Copy results into output matrices
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){

      if(i==j){   E->M[i*Norb+i] = solution.eigenvalues()[i].real();  }// i == j    

      C->M[j*Norb+i] = solution.eigenvectors().col(i)[j].real(); 

    }// for j
  }// for i

  if(0){
    // Checking solution - for debug purposes, inactive by default
    cout<<"in solve_eigen (with MATRIX):\n";
    cout<<"A*C = "<<*H * (*C)<<endl;
    cout<<"C*E = "<<(*C) * (*E)<<endl;
    cout<<"C.T()*C= "<<((*C).T()) * (*C)<<endl; // This gives unity, indeed.
  }

}// void solve_eigen(int Norb, MATRIX* H, MATRIX* E, MATRIX* C)


}// namespace libmeigen
}// namespace libmmath
