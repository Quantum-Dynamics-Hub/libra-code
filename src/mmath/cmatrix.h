/*********************************************************************************
* Copyright (C) 2012 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/


#ifndef cmatrix_h
#define cmatrix_h

#include <complex>
#include <vector>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>

using namespace std;


class matrix{

  void max_nondiagonal(int& row,int& col);

public:
  int n_rows,n_cols,n_elts;
  complex<double>* M;

  // Constructors
  matrix(){ n_rows = n_cols = n_elts = 0; M = NULL;}
  matrix(int n_rows_,int n_cols_){
    n_rows = n_rows_; n_cols = n_cols_; n_elts = n_rows * n_cols;
    M = new complex<double>[n_elts];
    for(int i=0;i<n_elts;i++){  M[i] = std::complex<double>(0.0,0.0); }
  }
  matrix(vector<vector<double> >& re_part,vector<vector<double> >& im_part);
 
  // Copy constructor
  matrix(const matrix& ob);  

  // Destructor
  ~matrix(){ delete [] M; n_rows = n_cols = n_elts = 0;}

  // Operatiions
  // Important! For memory efficiency it is crucial to use
  // (const matrix& ob) instead of (matrix ob) in some of the below
  // operators  - this avoid creation of the copy of the matrix object
  matrix operator-();                 // Negation;
  matrix operator*(const matrix& ob);
  matrix operator+(const matrix& ob);
  matrix operator-(const matrix& ob);
  void operator*=(const double&);
  void operator*=(const complex<double>&);
  void operator*=(const matrix& ob);
  void operator+=(const matrix& ob);
  void operator-=(const matrix& ob);

  matrix operator/(double num);
  matrix operator/(complex<double> num);
  matrix& operator=(const matrix& ob);
  matrix operator=(double num);
  matrix operator=(complex<double> num);


  friend matrix operator*(const double& f,  const matrix& m1);  // Multiplication of matrix and double;
  friend matrix operator*(const matrix& m1, const double  &f);  // Multiplication of matrix and double;
  friend matrix operator*(const complex<double>& f,  const matrix& m1);  // Multiplication of matrix and double;
  friend matrix operator*(const matrix& m1, const complex<double>  &f);  // Multiplication of matrix and double;
  friend ostream &operator<<(ostream &strm,matrix ob);
  friend istream& operator>>(istream& strm,matrix &ob);



  // Some basic functions
  matrix conj();
  matrix T();   // transpose
  matrix H();   // Hermitian conj
  void load_identity(); 
  void dot(const matrix& ob1,const matrix& ob2);

  matrix col(int); // takes given column and makes it n x 1 matrix
  matrix row(int); // takes given row and makes it 1 x n matrix

  // More advanced functions
  void QR(matrix& w,matrix& R);  // QR for general (Hermitian or symmetric) matrices
  void QR1(matrix& w,matrix& R); // QR for tridiagonal matrices
  void tridiagonalize(matrix& T);// for Hermitian or symmetric matrix - only resulting tridiagonal matrix
  void tridiagonalize(matrix& T,matrix& H); // ---//---  also keep track of Householder transformation matrices

  // Eigenvalues
  void eigen(double EPS,matrix& EVAL,matrix& EVECT,int opt); // interface
  void eigen0(matrix& EVAL, matrix& EVECT,double EPS,int max_num_iter,int is_cycle,int alg); // Schur decomposition or Jacobi rotation
  void eigen1(double EPS,vector<double>& Eval);  // only eigenvalues - fast
  void eigen2(double EPS,vector<double>& Eval,matrix& EVECT); // also eigenvectors, slower - keep track of transformation matrixes
  void eigen3(double EPS,vector<double>& Eval,matrix& EVECT); // also eigenvectors - solve for each eigenvector independently

  // Matrix inverse
  void inverse(double EPS,matrix& INV,int opt); // interface
  void direct_inverse(double EPS,matrix& INV);

  // Functions of matrix
  friend matrix exp(matrix& m1,complex<double> scl,double eps);
  friend matrix sin(matrix& m1,complex<double> scl,double eps);
  friend matrix cos(matrix& m1,complex<double> scl,double eps);
  friend matrix pow(matrix& m1,double scl,double eps);

};

void qr(double EPS,int n,matrix& M,vector<double>& Eval); // Compute eigenvalues of M
void qr(double EPS,int n,matrix& M,vector<double>& Eval,matrix& Evec);  // Computes eigenvalues and eigenvectors of M

void solve_linsys(matrix& C,matrix& D, matrix& X,double eps,int maxiter,double omega);
void solve_linsys1(matrix& C,matrix& X,double eps,int maxiter,double omega);


//---------- Fourier transforms ----------------
void dft(matrix& in,matrix& out);
void inv_dft(matrix& in,matrix& out);

void cft(matrix& in,matrix& out,double xmin,double dx);
void inv_cft(matrix& in,matrix& out,double xmin,double dx);

void cft1(matrix& in,matrix& out,double xmin,double kmin,double dx);
void inv_cft1(matrix& in,matrix& out,double xmin,double kmin,double dx);

void cft2(matrix& in,matrix& out,double xmin,double kmin,double dx,double dk);
void inv_cft2(matrix& in,matrix& out,double xmin,double kmin,double dx,double dk);

void cft1_2D(matrix& in, matrix& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);
void inv_cft1_2D(matrix& in, matrix& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);

void convolve(matrix& f,matrix& g, matrix& conv,double dx);
void convolve_2D(matrix& f,matrix& g, matrix& conv,double dx,double dy);

//-------- Fast Fourier Transforms -------------
void cfft1(matrix& in,matrix& out,double xmin,double kmin,double dx);
void inv_cfft1(matrix& in,matrix& out,double xmin,double kmin,double dx);

void cfft1_2D(matrix& in, matrix& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);
void inv_cfft1_2D(matrix& in, matrix& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);


#endif // cmatrix_h
