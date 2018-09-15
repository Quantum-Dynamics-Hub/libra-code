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
/**
  \file CMATRIX.cpp
  \brief The file implements the CMATRIX class for representing arbitrary size complex-valued matrices as well as the
  set of functions of complex Fourier transforms and convolution
    
*/

#include "CMATRIX.h"
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <iostream>



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
/** Returns the matrix which is transposed w.r.t. the caller matrix */

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
/** takes given column and makes it n x 1 CMATRIX */

  CMATRIX tmp(n_rows,1);
  for(int j=0;j<n_rows;j++){ tmp.M[j] = M[j*n_cols+i]; }
  return tmp;
}

CMATRIX CMATRIX::row(int i){ 
/** takes given row and makes it 1 x n CMATRIX */

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

  double x = abs(M[0]);
  double y;
  for(int i=0;i<n_elts;i++){  y = abs(M[i]); if(y>=x){ x = y; } }
  return x;

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






vector<int> get_reordering(CMATRIX& time_overlap){
    /**
    """ This function identifies which states have changed their identities during some
    calculations (usually the eigenvalue problem) in comparison to what they might have been.

    In the dynamics, this situation occurs when the system passes the conical intersection region
    and the identity of the states may change. The states' energies become non-descriptive for this
    purpose and one needs to look at the changes of the orbitals (e.g. eigenvectors). In particular,
    we can look at the overlap of the sates at adjacent times: <phi_i(t)|phi_i(t+dt)> 
    If no spurious state changes happens, the diagonal elements should be close to 1.0. 
    If they are not - we locate to which state the transitions might have happened.

    In general context, the "time_overlap" matrix is compused of the overlaps of the eigenvectors
    for two problems - the original one and a perturbed one.

    \param[in] time_overlap ( CMATRIX ) the time overlap matrix, <phi_i(t)|phi_j(t+dt)>.

    Returns:
    perm - list of integers that describe the permutation. That is:
    perm[i] - is the index identifying the "older" state "i". Now, it may be labeled
    by some other index, j = perm[i]. 

    """
    */

    // extract the indices where <phi_i(t)|phi_i(t+dt)> is not close to 1. 
    CMATRIX S(time_overlap);  // just a temporary working object
    int sz = time_overlap.n_rows;
    int i;

    // Original permutation
    vector<int> perm(sz, 0);  
    vector<int> perm_cum(sz, 0);  
    for(i=0;i<sz;i++){ perm_cum[i] = i; } 
    
    for(int col=0; col<sz; col++){

      int indx = -1;
      complex<double> val(0.0, 0.0);
      
      while(indx!=col){

        // Find the max element in the given column "col"
        S.max_col_elt(col, val, indx);
            
        // Apply the permutation (col, indx) to the present "perm" list
        for(i=0;i<sz;i++){ perm[i] = i; } 

        int tmp = perm[col];
        perm[col] = perm[indx];
        perm[indx] = tmp;

        // Do the corresponding swap of the columns in the S matrix
        S.swap_cols(col,indx);

        update_permutation(perm, perm_cum);


      }// while indx!=col
    }// for col

    return perm_cum;
}


vector<int> compute_signature(CMATRIX& Ref, CMATRIX& X){
/**
  Compute signature of one matrix with respect to a given reference matrix
  This is essentially determining whether the system of vectors given by
  the columns of the matrix X forms a "right" or "left" system in multi-
  dimensional space of vectors. What is "right" and what is "left" is
  given by the matrix Ref.

  Essentially, we compute a projection of each column of X onto each of the columns 
  of Ref and define the sign on that projection. Positive sign is a signature +1 for that
  dimension, whereas negative sign gives -1 signature. A (hopefully) rare case of the projection
  equal 0 is not yet handled, but if this will happen, we will need to simply rotate the 
  reference system of vectors by a random (small) angle in a random direction to make a non-zero 
  projection onto each of the directions.
*/

  // Check that the dimensions of the two matrices are equal
  if(Ref.n_cols!=X.n_cols){
    cout<<"ERROR in vector<int> compute_signature(CMATRIX& Ref, CMATRIX& X):\
          The number of columns of the reference matrix ("<<Ref.n_cols<<") \
          should be equal to the number of columns of the matrix of interest \
          ("<<X.n_cols<<")\nExiting...\n";
    exit(0);
  }

  if(Ref.n_rows!=X.n_rows){
    cout<<"ERROR in vector<int> compute_signature(CMATRIX& Ref, CMATRIX& X):\
          The number of rows of the reference matrix ("<<Ref.n_rows<<") \
          should be equal to the number of rows of the matrix of interest \
          ("<<X.n_rows<<")\nExiting...\n";
    exit(0);
  }

  vector<int> res(X.n_cols, 1);

  CMATRIX tmp(X.n_cols, X.n_cols);
  tmp = Ref.H() * X;

  for(int n=0; n<X.n_cols; n++){

    int val = tmp.get(n,n).real();


    if(fabs(val)<1e-50){   
/*
      cout<<"ERROR in vector<int> compute_signature(CMATRIX& Ref, CMATRIX& X):\
            One of the projections of the X matrix onto one of the \
            reference vectors is zero (less thant 1e-50)\nExiting...\n";
      exit(0);
*/
    }
    else{
      if(val<0.0){ res[n] = -1; }
    }  
  }// for n

  return res;

}

vector<int> compute_signature(CMATRIX& X){

  if(X.n_cols != X.n_rows){
    cout<<"ERROR in vector<int> compute_signature(CMATRIX& X): The number of columns (currently = \
    "<<X.n_cols<<") of matrix X should be equal to the number of rows (currently = "<<X.n_rows<<")\
    Exiting...\n";
    exit(0);
  }
  CMATRIX ref(X.n_cols, X.n_cols); 
  ref.identity();

  return compute_signature(ref, X);

}




void correct_phase(CMATRIX& Ref, CMATRIX& X){
/**
  This function checks that the phases of all vectors in matrix X are consistent with the
  phases of the same vectors in the reference set of vectors Ref.
  Correct the phases of vectors in X if needed.
*/

  vector<int> res(1, X.n_cols);
  res = compute_signature(Ref, X);

  for(int c=0; c<X.n_cols; c++){
    for(int r=0; r<X.n_rows; r++){
        X.scale(r,c, res[c]); 
    }
  }
}


void correct_phase(CMATRIX& X){
/**
  This function checks that the phases of all vectors in matrix X are consistent with the
  phases of the same vectors in the reference set of vectors Ref.
  Correct the phases of vectors in X if needed.
*/

  vector<int> res(1, X.n_cols);
  res = compute_signature(X);

  for(int c=0; c<X.n_cols; c++){
    for(int r=0; r<X.n_rows; r++){
        X.scale(r,c, res[c]); 
    }
  }
}



void correct_phase(CMATRIX& Ref, CMATRIX* X){
/**
  This function checks that the phases of all vectors in matrix X are consistent with the
  phases of the same vectors in the reference set of vectors Ref.
  Correct the phases of vectors in X if needed.
*/

  vector<int> res(1, X->n_cols);
  res = compute_signature(Ref, *X);

  for(int c=0; c<X->n_cols; c++){
    for(int r=0; r<X->n_rows; r++){
        X->scale(r,c, res[c]); 
    }
  }
}

void correct_phase(CMATRIX* X){
/**
  This function checks that the phases of all vectors in matrix X are consistent with the
  phases of the same vectors in the reference set of vectors Ref.
  Correct the phases of vectors in X if needed.
*/

  vector<int> res(1, X->n_cols);
  res = compute_signature(*X);

  for(int c=0; c<X->n_cols; c++){
    for(int r=0; r<X->n_rows; r++){
        X->scale(r,c, res[c]); 
    }
  }
}








}// namespace liblinalg
}// liblibra

