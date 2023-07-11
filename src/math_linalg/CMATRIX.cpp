/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file CMATRIX.cpp
  \brief The file implements the CMATRIX class for representing arbitrary size complex-valued matrices as well as the
  set of functions of complex Fourier transforms and convolution
    
*/


#if defined(USING_PCH)
#include "../pch.h"
#else
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <iostream>
#endif

#include "CMATRIX.h"


/// liblibra 
namespace liblibra{

using namespace std;

/// liblinalg namespace
namespace liblinalg{


CMATRIX::CMATRIX(vector<vector<double> >& re_part,vector<vector<double> >& im_part){
/**
  Constructor: creates a CMATRIX from 2 2D-arrays - real and imaginary parts
*/

  if(re_part.size()!=im_part.size()){ 
    cout<<"Error in CMATRIX constructor: y-dimensions(num of rows) of the real and imaginary arrays are not equal\n"; exit(0); 
  }
  if(re_part[0].size()!=im_part[0].size()){
    cout<<"Error in CMATRIX constructor: x-dimensions(num of cols) of the real and imaginary arrays are not equal\n"; exit(0);
  }
  n_rows = re_part.size();
  n_cols = re_part[0].size();
  n_elts = n_rows * n_cols;

  M = new complex<double>[n_elts];
  int n = 0;
  for(int i=0;i<n_rows;i++){ 
    for(int j=0;j<n_cols;j++){
      M[n] = complex<double>(re_part[i][j],im_part[i][j]); n++;
    }
  }
 
}

/*
CMATRIX::CMATRIX(IMATRIX& re_part){

  n_rows = re_part.n_rows;
  n_cols = re_part.n_cols;
  n_elts = n_rows * n_cols;

  M = new complex<double>[n_elts];
  int n = 0;
  for(int i=0;i<n_rows;i++){ 
    for(int j=0;j<n_cols;j++){
      M[n] = complex<double>((double)re_part.get(i,j), 0.0); n++;
    }
  }

}

CMATRIX::CMATRIX(const IMATRIX& re_part){

  n_rows = re_part.n_rows;
  n_cols = re_part.n_cols;
  n_elts = n_rows * n_cols;

  M = new complex<double>[n_elts];
  int n = 0;
  for(int i=0;i<n_rows;i++){ 
    for(int j=0;j<n_cols;j++){
      M[n] = complex<double>((double)re_part.get(i,j), 0.0); n++;
    }
  }

}
*/

CMATRIX::CMATRIX(MATRIX& re_part){

  n_rows = re_part.n_rows;
  n_cols = re_part.n_cols;
  n_elts = n_rows * n_cols;

  M = new complex<double>[n_elts];
  int n = 0;
  for(int i=0;i<n_rows;i++){ 
    for(int j=0;j<n_cols;j++){
      M[n] = complex<double>(re_part.get(i,j), 0.0); n++;
    }
  }

}


CMATRIX::CMATRIX(const MATRIX& re_part){

  n_rows = re_part.n_rows;
  n_cols = re_part.n_cols;
  n_elts = n_rows * n_cols;

  M = new complex<double>[n_elts];
  int n = 0;
  for(int i=0;i<n_rows;i++){ 
    for(int j=0;j<n_cols;j++){
      M[n] = complex<double>(re_part.get(i,j), 0.0); n++;
    }
  }

}




CMATRIX::CMATRIX(MATRIX& re_part,MATRIX& im_part){

  if(re_part.n_rows!=im_part.n_rows){ 
    cout<<"Error in CMATRIX constructor: num_of_rows of the real and imaginary matrices are not equal\n"; exit(0); 
  }
  if(re_part.n_cols!=im_part.n_cols){
    cout<<"Error in CMATRIX constructor: num_of_cols of the real and imaginary arrays are not equal\n"; exit(0);
  }
  n_rows = re_part.n_rows;
  n_cols = re_part.n_cols;
  n_elts = n_rows * n_cols;

  M = new complex<double>[n_elts];
  int n = 0;
  for(int i=0;i<n_rows;i++){ 
    for(int j=0;j<n_cols;j++){
      M[n] = complex<double>(re_part.get(i,j),im_part.get(i,j)); n++;
    }
  }

}

CMATRIX CMATRIX::T(){   
// Returns the matrix which is transposed w.r.t. the caller matrix 

  CMATRIX res(*this); res.Transpose();
  return res;    
}

CMATRIX CMATRIX::H(){
/** Returns the matrix which is Hermitian-conjugate to the caller matrix */

  CMATRIX m(n_cols,n_rows);
  for(int i=0;i<n_rows;i++){
    for(int j=0;j<n_cols;j++){
      m.M[j*n_rows+i] = std::conj(M[i*n_cols+j]);
    }
  }
  return m;
}


CMATRIX CMATRIX::conj(){
/** Returns the matrix which is complex-conjugate to the caller matrix */

  CMATRIX m(n_rows,n_cols);
  for(int i=0;i<m.n_elts;i++){ m.M[i] = std::conj(M[i]); }
  return m;
}


MATRIX CMATRIX::real(){
  MATRIX res(n_rows, n_cols);
  for(int i=0;i<n_elts;i++){ res.M[i] = M[i].real(); }
  return res;
}

MATRIX CMATRIX::imag(){
  MATRIX res(n_rows, n_cols);
  for(int i=0;i<n_elts;i++){ res.M[i] = M[i].imag(); }
  return res;
}

void CMATRIX::get_components(MATRIX& re_part,MATRIX& im_part){
  if(re_part.n_cols != n_cols){
    std::cout<<"Error in CMATRIX::get_components : The number of columns of the target real component matrix ("
             <<re_part.n_cols<<") is not equal to the number of columns of the source complex matrix ("<<n_cols<<")\n";
    exit(0);
  }
  if(im_part.n_cols != n_cols){
    std::cout<<"Error in CMATRIX::get_components : The number of columns of the target imaginary component matrix ("
             <<im_part.n_cols<<") is not equal to the number of columns of the source complex matrix ("<<n_cols<<")\n";
    exit(0);
  }
  if(re_part.n_rows != n_rows){
    std::cout<<"Error in CMATRIX::get_components : The number of rows of the target real component matrix ("
             <<re_part.n_rows<<") is not equal to the number of rows of the source complex matrix ("<<n_rows<<")\n";
    exit(0);
  }
  if(im_part.n_rows != n_rows){
    std::cout<<"Error in CMATRIX::get_components : The number of rows of the target imaginary component matrix ("
             <<im_part.n_rows<<") is not equal to the number of rows of the source complex matrix ("<<n_rows<<")\n";
    exit(0);
  }

  for(int i=0;i<n_elts;i++){ re_part.M[i] = M[i].real();  im_part.M[i] = M[i].imag(); }

}



CMATRIX CMATRIX::col(int i){ 
// takes given column and makes it n x 1 CMATRIX 

  CMATRIX tmp(n_rows,1);
  for(int j=0;j<n_rows;j++){ tmp.M[j] = M[j*n_cols+i]; }
  return tmp;
}

CMATRIX CMATRIX::row(int i){ 
// takes given row and makes it 1 x n CMATRIX 

  CMATRIX tmp(1,n_cols);
  for(int j=0;j<n_cols;j++){ tmp.M[j] = M[i*n_cols+j]; }
  return tmp;
}




double CMATRIX::NonOrtogonality_Measure(){
/** The sum of the scalar products of all pairs of distinct columns
*/
  double sum;   sum = 0.0;   
  complex<double> aa;

  for(int i=0;i<n_cols-1;i++){
    for(int j=i+1;j<n_cols;j++){

      aa = 0.0;
      for(int k=0;k<n_rows;k++){ aa += std::conj(M[k*n_cols+i]) * M[k*n_cols+j];  }
      sum += (std::conj(aa) * aa).real();
    }
  }
  return sum;
}

complex<double> CMATRIX::max_elt(){
/** Finds the maximal (in absolute value) element and its position  */

  int max_indx = 0;
  double x = abs(M[0]);
  double y;
  for(int i=0;i<n_elts;i++){  
    y = abs(M[i]); 
    if(y>=x){ x = y; max_indx = i; } 
  }

  return M[max_indx];

}

void CMATRIX::FindMaxNondiagonalElement(int& row,int& col,complex<double>& value){
/** Finds the maximal (in absolute value) off-diagonal element and its position  */

  int k=0;
  double elem, max_elem;
  value = M[1]; max_elem = abs(value); row = 0 ; col = 1;

  for(int rw=0;rw<n_rows;rw++){
    for(int cl=rw+1;cl<n_cols;cl++){
      k = rw*n_cols + cl;
      elem = abs(M[k]);
      if(elem > max_elem) {max_elem = elem; value = M[k]; col = cl; row = rw;}
    }
  }
}



void CMATRIX::max_nondiagonal(int& row,int& col){
  double maxeps = norm(M[1]); row = 0; col = 1;
  double eps;
  for(int r=0;r<n_rows;r++){
    for(int c=r+1;c<n_cols;c++){
      eps = norm(M[r*n_cols+c]);
      if(eps>=maxeps){ row = r; col = c; maxeps = eps; }
    }
  }  
}





void CMATRIX::max_col_elt(int I, complex<double>& val, int& max_elt_indx){
 ///< Finds the maximal element (in abs. value) and its index in a given column

  val = M[0*n_cols+I];
  max_elt_indx = 0;
  double max_norm = norm(val);
  
  for(int row=0; row<n_rows; row++){

    double nrm = norm(M[row*n_cols+I]);
    if(nrm>max_norm){  
      max_norm = nrm;   
      val = M[row*n_cols+I];
      max_elt_indx = row;
    }

  }// for row

}

void CMATRIX::min_col_elt(int I, complex<double>& val, int& min_elt_indx){
 ///< Finds the minimal element (in abs. value) and its index in a given column

  val = M[0*n_cols+I];
  min_elt_indx = 0;
  double min_norm = norm(val);
  
  for(int row=0; row<n_rows; row++){

    double nrm = norm(M[row*n_cols+I]);
    if(nrm < min_norm){  
      min_norm = nrm;   
      val = M[row*n_cols+I];
      min_elt_indx = row;
    }

  }// for row

}



void CMATRIX::max_row_elt(int I, complex<double>& val, int& max_elt_indx){
 ///< Finds the maximal element (in abs. value) and its index in a given row

  val = M[I*n_cols+0];
  max_elt_indx = 0;
  double max_norm = norm(val);
  
  for(int col=0; col<n_cols; col++){

    double nrm = norm(M[I*n_cols+col]);
    if(nrm>max_norm){  
      max_norm = nrm;   
      val = M[I*n_cols+col];
      max_elt_indx = col;
    }

  }// for col

}

void CMATRIX::min_row_elt(int I, complex<double>& val, int& min_elt_indx){
 ///< Finds the minimal element (in abs. value) and its index in a given row

  val = M[I*n_cols+0];
  min_elt_indx = 0;
  double min_norm = norm(val);
  
  for(int col=0; col<n_cols; col++){

    double nrm = norm(M[I*n_cols+col]);
    if(nrm<min_norm){  
      min_norm = nrm;   
      val = M[I*n_cols+col];
      min_elt_indx = col;
    }

  }// for col

}


boost::python::list CMATRIX::max_col_elt(int I){
 ///< Finds the maximal element (in abs. value) and its index in a given column
  complex<double> val;
  int indx;

  this->max_col_elt(I, val, indx);
 
  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;
}

boost::python::list CMATRIX::min_col_elt(int I){
  ///< Finds the maximal element (in abs. value) and its index in a given column

  complex<double> val;
  int indx;

  this->min_col_elt(I, val, indx);
 
  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;

}

boost::python::list CMATRIX::max_row_elt(int I){
 ///< Finds the maximal element (in abs. value) and its index in a given row

  complex<double> val;
  int indx;

  this->max_row_elt(I, val, indx);
 
  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;

}

boost::python::list CMATRIX::min_row_elt(int I){
 ///< Finds the maximal element (in abs. value) and its index in a given row

  complex<double> val;
  int indx;

  this->min_row_elt(I, val, indx);
 
  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;

}


CMATRIX CMATRIX::operator+(const CMATRIX& rhs){ 
  CMATRIX res(*this);  res += rhs;
  return res;
}

///< Addition operator
CMATRIX CMATRIX::operator+(int rhs){ 
  CMATRIX res(*this);  res += rhs;
  return res;
}

CMATRIX CMATRIX::operator+(double rhs){ 
  CMATRIX res(*this);  res += rhs;
  return res;
}

CMATRIX CMATRIX::operator+(complex<double> rhs){ 
  CMATRIX res(*this);  res += rhs;
  return res;
}


CMATRIX CMATRIX::operator-(const CMATRIX& rhs){ 
  CMATRIX res(*this);  res -= rhs;
  return res;
}

///< Addition operator
CMATRIX CMATRIX::operator-(int rhs){ 
  CMATRIX res(*this);  res -= rhs;
  return res;
}

CMATRIX CMATRIX::operator-(double rhs){ 
  CMATRIX res(*this);  res -= rhs;
  return res;
}

CMATRIX CMATRIX::operator-(complex<double> rhs){ 
  CMATRIX res(*this);  res -= rhs;
  return res;
}


void CMATRIX::operator+=(complex<double> f){  
  for(int i=0;i<n_elts;i++) { M[i] += f; }
}
void CMATRIX::operator-=(complex<double> f){  
  for(int i=0;i<n_elts;i++) { M[i] -= f; }
}


CMATRIX CMATRIX::operator*(const CMATRIX& ob){
  CMATRIX res(n_rows, ob.n_cols);  res.product(*this, ob);
  return res;
}

CMATRIX operator*(const CMATRIX& ob, int f){  
  CMATRIX res(ob);  res *= f;
  return res;
}

CMATRIX operator*(const CMATRIX& ob, double f){  
  CMATRIX res(ob);  res *= f;
  return res;
}

CMATRIX operator*(const CMATRIX& ob, complex<double> f){  
  CMATRIX res(ob);  res *= f;
  return res;
}



CMATRIX operator*(int f, const CMATRIX& ob){  
  CMATRIX res(ob);   res *= f;
  return res;
}

CMATRIX operator*(double f, const CMATRIX& ob){  
  CMATRIX res(ob);   res *= f;
  return res;
}

CMATRIX operator*(complex<double> f, const CMATRIX& ob){  
  CMATRIX res(ob);   res *= f;
  return res;
}


/*
CMATRIX operator*(const MATRIX& mtx1, const CMATRIX& mtx2){
  CMATRIX res(mtx1); res = res * mtx2;
  return res;
}

CMATRIX operator*(const CMATRIX& mtx1, const MATRIX& mtx2){  
  CMATRIX res(mtx1); res = res * CMATRIX(mtx2);  
  return res;
}


CMATRIX operator*(const IMATRIX& mtx1, const CMATRIX& mtx2){
  CMATRIX res(mtx1); res = res * mtx2;
  return res;
}

CMATRIX operator*(const CMATRIX& mtx1, const IMATRIX& mtx2){
  CMATRIX res(mtx1); res = res * CMATRIX(mtx2);
  return res;
}
*/

void CMATRIX::operator*=(complex<double> f){ 
  for(int i=0; i<n_elts; i++){  M[i] *= f;  }
}

void CMATRIX::operator/=(complex<double> f){ 
  for(int i=0; i<n_elts; i++){  M[i] /= f;  }
}


CMATRIX CMATRIX::operator/(int f){ 
  CMATRIX res(*this); res /= f; 
  return res;
}

CMATRIX CMATRIX::operator/(double f){ 
  CMATRIX res(*this); res /= f; 
  return res;
}


CMATRIX CMATRIX::operator/(complex<double> f){ 
  CMATRIX res(*this); res /= f; 
  return res;
}








}// namespace liblinalg
}// liblibra

