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


#ifndef CMATRIX_H
#define CMATRIX_H

#include <complex>
#include <vector>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>

using namespace std;

namespace libmmath{
namespace liblinalg{

class CMATRIX{

  void max_nondiagonal(int& row,int& col);

public:
  int n_rows,n_cols,n_elts;
  complex<double>* M;

  void set(int,double,double);
  complex<double> get(int);
  void set(int,int,double,double);
  complex<double> get(int,int);


  // Constructors
  CMATRIX(){ n_rows = n_cols = n_elts = 0; M = NULL;}
  CMATRIX(int n_rows_,int n_cols_){
    n_rows = n_rows_; n_cols = n_cols_; n_elts = n_rows * n_cols;
    M = new complex<double>[n_elts];
    for(int i=0;i<n_elts;i++){  M[i] = std::complex<double>(0.0,0.0); }
  }
  CMATRIX(vector<vector<double> >& re_part,vector<vector<double> >& im_part);
 
  // Copy constructor
  CMATRIX(const CMATRIX& ob);  

  // Destructor
  ~CMATRIX(){ delete [] M; n_rows = n_cols = n_elts = 0;}

  // Operatiions
  // Important! For memory efficiency it is crucial to use
  // (const CMATRIX& ob) instead of (CMATRIX ob) in some of the below
  // operators  - this avoid creation of the copy of the CMATRIX object
  CMATRIX operator-();                 // Negation;
  CMATRIX operator*(const CMATRIX& ob);
  CMATRIX operator+(const CMATRIX& ob);
  CMATRIX operator-(const CMATRIX& ob);
  void operator*=(const double&);
  void operator*=(const complex<double>&);
  void operator*=(const CMATRIX& ob);
  void operator+=(const CMATRIX& ob);
  void operator-=(const CMATRIX& ob);

  CMATRIX operator/(double num);
  CMATRIX operator/(complex<double> num);
  CMATRIX& operator=(const CMATRIX& ob);
  CMATRIX operator=(double num);
  CMATRIX operator=(complex<double> num);


  friend int operator == (const CMATRIX& m1, const CMATRIX& m2);  // Are matrices equal;

  friend CMATRIX operator*(const double& f,  const CMATRIX& m1);  // Multiplication of CMATRIX and double;
  friend CMATRIX operator*(const CMATRIX& m1, const double  &f);  // Multiplication of CMATRIX and double;
  friend CMATRIX operator*(const complex<double>& f,  const CMATRIX& m1);  // Multiplication of CMATRIX and double;
  friend CMATRIX operator*(const CMATRIX& m1, const complex<double>  &f);  // Multiplication of CMATRIX and double;
  friend ostream &operator<<(ostream &strm,CMATRIX ob);
  friend istream& operator>>(istream& strm,CMATRIX &ob);

  // Basic printing
  void show();

  // Some basic functions
  CMATRIX conj();
  CMATRIX T();   // transpose
  CMATRIX H();   // Hermitian conj
  void load_identity(); 
  void dot(const CMATRIX& ob1,const CMATRIX& ob2);

  CMATRIX col(int); // takes given column and makes it n x 1 CMATRIX
  CMATRIX row(int); // takes given row and makes it 1 x n CMATRIX

  // More advanced functions
  void QR(CMATRIX& w,CMATRIX& R);  // QR for general (Hermitian or symmetric) matrices
  void QR1(CMATRIX& w,CMATRIX& R); // QR for tridiagonal matrices
  void tridiagonalize(CMATRIX& T);// for Hermitian or symmetric CMATRIX - only resulting tridiagonal CMATRIX
  void tridiagonalize(CMATRIX& T,CMATRIX& H); // ---//---  also keep track of Householder transformation matrices

  // Eigenvalues
  void eigen(double EPS,CMATRIX& EVAL,CMATRIX& EVECT,int opt); // interface
  void eigen0(CMATRIX& EVAL, CMATRIX& EVECT,double EPS,int max_num_iter,int is_cycle,int alg); // Schur decomposition or Jacobi rotation
  void eigen1(double EPS,vector<double>& Eval);  // only eigenvalues - fast
  void eigen2(double EPS,vector<double>& Eval,CMATRIX& EVECT); // also eigenvectors, slower - keep track of transformation matrixes
  void eigen3(double EPS,vector<double>& Eval,CMATRIX& EVECT); // also eigenvectors - solve for each eigenvector independently

  // Matrix inverse
  void inverse(double EPS,CMATRIX& INV,int opt); // interface
  void direct_inverse(double EPS,CMATRIX& INV);

  // Functions of CMATRIX
  friend CMATRIX exp(CMATRIX& m1,complex<double> scl,double eps);
  friend CMATRIX sin(CMATRIX& m1,complex<double> scl,double eps);
  friend CMATRIX cos(CMATRIX& m1,complex<double> scl,double eps);
  friend CMATRIX pow(CMATRIX& m1,double scl,double eps);

};

void qr(double EPS,int n,CMATRIX& M,vector<double>& Eval); // Compute eigenvalues of M
void qr(double EPS,int n,CMATRIX& M,vector<double>& Eval,CMATRIX& Evec);  // Computes eigenvalues and eigenvectors of M

void solve_linsys(CMATRIX& C,CMATRIX& D, CMATRIX& X,double eps,int maxiter,double omega);
void solve_linsys1(CMATRIX& C,CMATRIX& X,double eps,int maxiter,double omega);


//---------- Fourier transforms ----------------
void dft(CMATRIX& in,CMATRIX& out);
void inv_dft(CMATRIX& in,CMATRIX& out);

void cft(CMATRIX& in,CMATRIX& out,double xmin,double dx);
void inv_cft(CMATRIX& in,CMATRIX& out,double xmin,double dx);

void cft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx);
void inv_cft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx);

void cft2(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx,double dk);
void inv_cft2(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx,double dk);

void cft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);
void inv_cft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);

void convolve(CMATRIX& f,CMATRIX& g, CMATRIX& conv,double dx);
void convolve_2D(CMATRIX& f,CMATRIX& g, CMATRIX& conv,double dx,double dy);

//-------- Fast Fourier Transforms -------------
void cfft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx);
void inv_cfft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx);

void cfft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);
void inv_cfft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);



  typedef std::vector<CMATRIX> CMATRIXList;

}//namespace liblinalg
}// namespace libmmath

#endif // CMATRIX_H
