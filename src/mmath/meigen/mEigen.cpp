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


double det(MATRIX& A){

  // Wrapper matrices for Eigen3

  if(A.num_of_cols!=A.num_of_rows){
    std::cout<<"Error in det(MATRIX): Can not compute a determinant of non-square matrix\n";
    exit(0);
  }

  int N = A.num_of_cols;
  int i,j;

  MatrixXd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i

  return a.determinant();


}

double FullPivLU_det(MATRIX& A){


  // Wrapper matrices for Eigen3

  if(A.num_of_cols!=A.num_of_rows){
    std::cout<<"Error in det(MATRIX): Can not compute a determinant of non-square matrix\n";
    exit(0);
  }

  int N = A.num_of_cols;
  int i,j;

  MatrixXd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i

  Eigen::FullPivLU<MatrixXd> lu(a);

  return lu.determinant();


}



complex<double> det(CMATRIX& A){

  // Wrapper matrices for Eigen3

  if(A.n_cols!=A.n_rows){
    std::cout<<"Error in det(MATRIX): Can not compute a determinant of non-square matrix\n";
    exit(0);
  }

  int N = A.n_cols;
  int i,j;

  MatrixXcd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i

  return a.determinant();

}

complex<double> FullPivLU_det(CMATRIX& A){

  // Wrapper matrices for Eigen3

  if(A.n_cols!=A.n_rows){
    std::cout<<"Error in det(MATRIX): Can not compute a determinant of non-square matrix\n";
    exit(0);
  }

  int N = A.n_cols;
  int i,j;

  MatrixXcd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i


  Eigen::FullPivLU<MatrixXcd> lu(a);

  return lu.determinant();


}




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



void sqrt_matrix(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half, double thresh){
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
  CMATRIX* I; I = new CMATRIX(sz, sz);  I->load_identity(); //->Init_Unit_Matrix(complex<double>(1.0,0.0));
  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0, 0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);


  // Find the eigenvalues of the the S matrix
  libmmath::libmeigen::solve_eigen(sz, S, *I, *Seig, *C);  // S * C = I * C * Seig  ==>  S = C * Seig * C.H()

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
  delete I;
  delete Seig;


}// sqrt_matrix


void sqrt_matrix(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half){

  sqrt_matrix(S, S_half, S_i_half, -1.0);

}


void inv_matrix(CMATRIX& S, CMATRIX& S_inv, double thresh){
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
  CMATRIX* I; I = new CMATRIX(sz, sz);  I->load_identity(); //->Init_Unit_Matrix(complex<double>(1.0,0.0));
  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0, 0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);


  // Find the eigenvalues of the the S matrix
  libmmath::libmeigen::solve_eigen(sz, S, *I, *Seig, *C);  // S * C = I * C * Seig  ==>  S = C * Seig * C.H()

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
  delete I;
  delete Seig;


}// inv_matrix


void inv_matrix(CMATRIX& S, CMATRIX& S_inv){

  inv_matrix(S, S_inv, -1.0);

}


void FullPivLU_rank_invertible(MATRIX& A, int& rank, int& is_inver){

  int N = A.num_of_cols;
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





void FullPivLU_decomposition(MATRIX& A, MATRIX& P, MATRIX& L, MATRIX& U, MATRIX& Q){

  int N = A.num_of_cols;
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
/**
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


void FullPivLU_inverse(MATRIX& A, MATRIX& invA){

  int N = A.num_of_cols;
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

/**
  This function is adapted from the OpenBabel Qeq method:

  Solving a linear system of equations: A * x = b

  The result goes into the x argument

  //! Wrapper around the Eigen linear solver routines
  // First attempts to solve using Gaussian elimination
  // if that fails, tries again using singular value decomposition (SVD)
*/
bool linsys_solver(const MATRIX& A, MATRIX& X, const MATRIX& B, const double NormThreshold){

    int i,j;   // counters

    if(A.num_of_rows!=A.num_of_cols){
        std::cout<<"Error: The matrix A in equation A * x = b must be square\n"; exit(1);  // N
    }
    if(A.num_of_cols!=X.num_of_rows){
        std::cout<<"Error: The number of cols of matrix A and the number of rows in matrix x in equation A * x = b must be equal\n"; exit(1);  // N
    }
    if( (X.num_of_cols!=B.num_of_cols) || (X.num_of_rows!=B.num_of_rows)){
        std::cout<<"Error: The dimensions of matrices x and b in equation A * x = b must be equal\n"; exit(1); 
    }
    if(X.num_of_cols!=1){
        std::cout<<"Error: The matrices x and b in equation A * x = b should be column-vectors (num_of_cols ==1) \n"; exit(1); 
    }

    // Set dimentions
    int N = A.num_of_rows; // the dimension of the matrix A (square)
/*

    // Convert out matrices to the Eigen types:
    MatrixXd a(N,N); for(i=0;i<N;i++){  for(j=0;j<N;j++){  a(i,j) = A.M[i*N+j];  }    }// for i
    VectorXd x(N);   for(i=0;i<N;i++){  x(i) = 0.0;   }// for i
    VectorXd b(N);   for(i=0;i<N;i++){  b(i) = B.M[i];   }// for i


    // using a LU factorization
    bool SolverOK = true;
    x = A.partialPivLu().solve(b);   //  bool SolverOK = A.lu().solve(b, &x);

    VectorXd resid = A*x - b;
    double resnorm = resid.norm();
    if(IsNan(resnorm) || resnorm > NormThreshold || !SolverOK){
        std::cout << "Warning, LU solver failed." << endl;
        if(!SolverOK){ std::cout << "Solver returned error." << endl; }
        if(IsNan(resnorm)){ std::cout << "NaNs were returned" << endl; }
        if(resnorm > NormThreshold){ std::cout << "Residual has norm " << resnorm
                                         << " which exceeds the recommended threshold of " << NormThreshold
                                         << endl;
        }
        std::cout << "Proceeding with singular value decomposition.";
               
        x = A.jacobiSvd().solve(b);   //     SolverOK = A.svd().solve(b, &x);
        resid = A*x - b;
        resnorm = resid.norm();

        if(IsNan(resnorm) || !SolverOK){
            std::cout << "SVD solver returned an error. Solution may be not reliable!\n";
            return false;
        }

    }// if go into Jacobi SVD

    std::cout << "The residual of the solution has norm " << resnorm << endl;

    if (resnorm > NormThreshold) {
      std::cout  << "Warning, the norm of the residual is " << resnorm
                 << "which exceeds the recommended threshold of " << NormThreshold << endl;
    }

    // Copy results from the intermediate variable back into the argument
    for(i=0;i<N;i++){  X.M[i] = x(i);   }// for i

*/

    return true;
}



}// namespace libmeigen
}// namespace libmmath
