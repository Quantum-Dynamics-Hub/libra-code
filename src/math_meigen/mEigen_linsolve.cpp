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



bool linsys_solver(const MATRIX& A, MATRIX& X, const MATRIX& B, const double NormThreshold){
/**
  This function is adapted from the OpenBabel Qeq method:

  Solving a linear system of equations: A * x = b

  The result goes into the x argument

  //! Wrapper around the Eigen linear solver routines
  // First attempts to solve using Gaussian elimination
  // if that fails, tries again using singular value decomposition (SVD)
*/

    std::cout<<"Not implemented. Exiting...\n";
    exit(0);

    int i,j;   // counters

    if(A.n_rows!=A.n_cols){
        std::cout<<"Error: The matrix A in equation A * x = b must be square\n"; exit(1);  // N
    }
    if(A.n_cols!=X.n_rows){
        std::cout<<"Error: The number of cols of matrix A and the number of rows in matrix x in equation A * x = b must be equal\n"; exit(1);  // N
    }
    if( (X.n_cols!=B.n_cols) || (X.n_rows!=B.n_rows)){
        std::cout<<"Error: The dimensions of matrices x and b in equation A * x = b must be equal\n"; exit(1); 
    }
    if(X.n_cols!=1){
        std::cout<<"Error: The matrices x and b in equation A * x = b should be column-vectors (num_of_cols ==1) \n"; exit(1); 
    }

    // Set dimentions
    int N = A.n_rows; // the dimension of the matrix A (square)

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


void solve_linsys(MATRIX& C,MATRIX& D, MATRIX& X,double eps,int maxiter){
/*********************************************
 Here we solve the system of linear equations
              AX = D
 using Gauss-Seidel iterative procedure

 Inputs: A, D - matrices
         eps  - precision criterion
 Output: X

 Some preliminary transformations are made in order
 to be able to use Gauss-Seidel method for any matrix A

 More details:
 80.47 An iterative Algorithm for Matrix Inversion
 Which is Always Convergent
 Authors: S. Simons
 Source: The Mathematical Gazette, Vol. 80, No. 489
 (Nov., 1996), pp. 567-569

 url: http://www.jstor.org/stable/pdfplus/3618529.pdf

**********************************************/

// Do the transformations A = C^T * C and b = C^T * d
// If matrices d and c have more then 1 columns we do the
// procedure for each column


    int i,j,k;   // counters
    int n,m,p;   // dimetions
    double s;    // sums
    double error;// error
    int iter;    // number of iterations

    if(C.n_rows!=D.n_rows)
        {std::cout<<"Error: The number of rows of matrices C and D in equation CX = D must be equal\n"; exit(35); } // n
    if(C.n_cols!=X.n_rows)
        {std::cout<<"Error: The number of cols of matrix C and num of rows in matrix D in equation CX = D must be equal\n"; exit(35); } // m
    if(X.n_cols!=D.n_cols)
        {std::cout<<"Error: The number of cols of matrices X and D in equation CX = D must be equal\n"; exit(35); } // p

    // Set dimentions
    n = C.n_rows;
    m = C.n_cols;
    p = D.n_cols;

    MATRIX A(m,m); A = C.T() * C;
    X = 0.0;
    error = 2.0*eps;
    iter = 0;

    while((error>eps)&&(iter<maxiter)){

    error = 0.0;

    for( k = 0; k < p; k++ ){

        //------- Matrix preparation step -----------

        MATRIX d(n,1);  for(i = 0;i<n;i++){ d.M[i] = D.M[i*p+k]; }
        MATRIX x(m,1);  for(i = 0;i<m;i++){ x.M[i] = X.M[i*p+k]; }
        MATRIX xprev(m,1); xprev = x;
        MATRIX b(m,1);  b = C.T() * d;

        //------- Gauss-Seidel step -----------------

        for( i = 0; i < m; i++ ){

            s = 0.0;
            for( j = 0; j < i; j++ ){ s += A.M[i*m + j]*x.M[j];   }// for j
            for( j = i+1; j < m; j++ ){ s += A.M[i*m + j]*xprev.M[j]; }// for j

            x.M[i] = (b.M[i] - s)/A.M[i*m + i];

        }// for i - all elements of vector x


        //-------- Now calculate the error and update X ---------
        for( i = 0; i < m; i++ ){
            error += (x.M[i] - xprev.M[i])*(x.M[i] - xprev.M[i]);
            X.M[i*p + k] = x.M[i];
        }// for i


    }// for k - all columns of D

    error = sqrt(error/double(m));
    iter++;

    }// loop over convergence

    if(error>eps){
      cout<<"Error in solve_linsys: convergence to eps= "<<eps<<" is not achieved for "<<iter<<" iterations\n";
      exit(0);
    }

}





}// namespace libmeigen
}// namespace liblibra
