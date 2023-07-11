/*********************************************************************************
* Copyright (C) 2020-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file IMATRIX.cpp
  \brief The file implements the IMATRIX class for representing arbitrary size integer-valued matrices     
*/


#if defined(USING_PCH)
#include "../pch.h"
#else
#include <cstdio>
#include <cstdlib>
#include <iostream>
#endif

#include "IMATRIX.h"

/// liblibra 
namespace liblibra{

using namespace std;

/// liblinalg namespace
namespace liblinalg{


IMATRIX::IMATRIX(vector<vector<int> >& mtx){
/**
  Constructor: creates a IMATRIX from a 2D-array
*/

  n_rows = mtx.size();
  n_cols = mtx[0].size();
  n_elts = n_rows * n_cols;

  M = new int[n_elts];
  int n = 0;
  for(int i=0;i<n_rows;i++){ 
    for(int j=0;j<n_cols;j++){
      M[n] = mtx[i][j]; n++;
    }
  }
 
}

/*
IMATRIX::IMATRIX(MATRIX& re_part){

  n_rows = re_part.n_rows;
  n_cols = re_part.n_cols;
  n_elts = n_rows * n_cols;

  M = new int[n_elts];
  int n = 0;
  for(int i=0;i<n_rows;i++){ 
    for(int j=0;j<n_cols;j++){
      M[n] = (int)re_part.get(i,j); n++;
    }
  }

}
*/


IMATRIX IMATRIX::T(){   
// Returns the matrix which is transposed w.r.t. the caller matrix 

  IMATRIX res(*this); res.Transpose();
  return res;    
}

IMATRIX IMATRIX::col(int i){ 
// takes given column and makes it n x 1 IMATRIX /

  IMATRIX tmp(n_rows,1);
  for(int j=0;j<n_rows;j++){ tmp.M[j] = M[j*n_cols+i]; }
  return tmp;
}

IMATRIX IMATRIX::row(int i){ 
// takes given row and makes it 1 x n IMATRIX /

  IMATRIX tmp(1,n_cols);
  for(int j=0;j<n_cols;j++){ tmp.M[j] = M[i*n_cols+j]; }
  return tmp;
}


int IMATRIX::max_elt(){
/** Finds the maximal (in absolute value) element and its position  */

  int max_indx = 0;
  int x = abs(M[0]);
  int y;
  for(int i=0;i<n_elts;i++){  
    y = abs(M[i]); 
    if(y>=x){ x = y; max_indx = i; } 
  }

  return M[max_indx];

}

void IMATRIX::FindMaxNondiagonalElement(int& row, int& col, int& value){
/** Finds the maximal (in absolute value) off-diagonal element and its position  */

  int k=0;
  int elem, max_elem;
  value = M[1]; max_elem = abs(value); row = 0 ; col = 1;

  for(int rw=0;rw<n_rows;rw++){
    for(int cl=rw+1;cl<n_cols;cl++){
      k = rw*n_cols + cl;
      elem = abs(M[k]);
      if(elem > max_elem) {max_elem = elem; value = M[k]; col = cl; row = rw;}
    }
  }
}



void IMATRIX::max_nondiagonal(int& row, int& col){
  int maxeps = abs(M[1]); row = 0; col = 1;
  int eps;
  for(int r=0;r<n_rows;r++){
    for(int c=r+1;c<n_cols;c++){
      eps = abs(M[r*n_cols+c]);
      if(eps>=maxeps){ row = r; col = c; maxeps = eps; }
    }
  }  
}





void IMATRIX::max_col_elt(int I, int& val, int& max_elt_indx){
 ///< Finds the maximal element (in abs. value) and its index in a given column

  val = M[0*n_cols+I];
  max_elt_indx = 0;
  int max_norm = abs(val);
  
  for(int row=0; row<n_rows; row++){

    int nrm = abs(M[row*n_cols+I]);
    if(nrm>max_norm){  
      max_norm = nrm;   
      val = M[row*n_cols+I];
      max_elt_indx = row;
    }

  }// for row

}

void IMATRIX::min_col_elt(int I, int& val, int& min_elt_indx){
 ///< Finds the minimal element (in abs. value) and its index in a given column

  val = M[0*n_cols+I];
  min_elt_indx = 0;
  int min_norm = abs(val);
  
  for(int row=0; row<n_rows; row++){

    int nrm = abs(M[row*n_cols+I]);
    if(nrm < min_norm){  
      min_norm = nrm;   
      val = M[row*n_cols+I];
      min_elt_indx = row;
    }

  }// for row

}



void IMATRIX::max_row_elt(int I, int& val, int& max_elt_indx){
 ///< Finds the maximal element (in abs. value) and its index in a given row

  val = M[I*n_cols+0];
  max_elt_indx = 0;
  int max_norm = abs(val);
  
  for(int col=0; col<n_cols; col++){

    int nrm = abs(M[I*n_cols+col]);
    if(nrm>max_norm){  
      max_norm = nrm;   
      val = M[I*n_cols+col];
      max_elt_indx = col;
    }

  }// for col

}

void IMATRIX::min_row_elt(int I, int& val, int& min_elt_indx){
 ///< Finds the minimal element (in abs. value) and its index in a given row

  val = M[I*n_cols+0];
  min_elt_indx = 0;
  int min_norm = abs(val);
  
  for(int col=0; col<n_cols; col++){

    int nrm = abs(M[I*n_cols+col]);
    if(nrm<min_norm){  
      min_norm = nrm;   
      val = M[I*n_cols+col];
      min_elt_indx = col;
    }

  }// for col

}


boost::python::list IMATRIX::max_col_elt(int I){
 ///< Finds the maximal element (in abs. value) and its index in a given column
  int val;
  int indx;

  this->max_col_elt(I, val, indx);
 
  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;
}

boost::python::list IMATRIX::min_col_elt(int I){
  ///< Finds the maximal element (in abs. value) and its index in a given column

  int val;
  int indx;

  this->min_col_elt(I, val, indx);
 
  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;

}

boost::python::list IMATRIX::max_row_elt(int I){
 ///< Finds the maximal element (in abs. value) and its index in a given row

  int val;
  int indx;

  this->max_row_elt(I, val, indx);
 
  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;

}

boost::python::list IMATRIX::min_row_elt(int I){
 ///< Finds the maximal element (in abs. value) and its index in a given row

  int val;
  int indx;

  this->min_row_elt(I, val, indx);
 
  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;

}


IMATRIX IMATRIX::operator+(const IMATRIX& rhs){ 
  IMATRIX res(*this);  res += rhs;
  return res;
}

///< Addition operator
IMATRIX IMATRIX::operator+(int rhs){ 
  IMATRIX res(*this);  res += rhs;
  return res;
}


IMATRIX IMATRIX::operator-(const IMATRIX& rhs){ 
  IMATRIX res(*this);  res -= rhs;
  return res;
}

///< Addition operator
IMATRIX IMATRIX::operator-(int rhs){ 
  IMATRIX res(*this);  res -= rhs;
  return res;
}


/*
void IMATRIX::operator+=(int f){  
  for(int i=0;i<n_elts;i++) { M[i] += f; }
}
*/

IMATRIX IMATRIX::operator*(const IMATRIX& ob){
  IMATRIX res(n_rows, ob.n_cols);  res.product(*this, ob);
  return res;
}

IMATRIX operator*(const IMATRIX& ob, int f){  
  IMATRIX res(ob);  res *= f;
  return res;
}



IMATRIX operator*(int f, const IMATRIX& ob){  
  IMATRIX res(ob);   res *= f;
  return res;
}


/*
void IMATRIX::operator*=(int f){ 
  for(int i=0; i<n_elts; i++){  M[i] *= f;  }
}

void IMATRIX::operator/=(int f){ 
  for(int i=0; i<n_elts; i++){  M[i] /= f;  }
}
*/

IMATRIX IMATRIX::operator/(int f){ 
  IMATRIX res(*this); res /= f; 
  return res;
}








}// namespace liblinalg
}// liblibra

