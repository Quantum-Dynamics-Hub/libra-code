/*********************************************************************************
* Copyright (C) 2015-2020 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file CMATRIX.h
  \brief The file describes the CMATRIX class for representing arbitrary size complex-valued matrices as well as the
  set of functions of complex Fourier transforms and convolution

  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  NOTE:  the functions that return CMATRIX should not be relpaced by the base template functions

  for some reason, the base_matrix< complex<double> > is not recognized as CMATRIX so 

  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
*/


#ifndef CMATRIX_H
#define CMATRIX_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <complex>
#endif 

#include "base_matrix.h"
#include "MATRIX.h"
//#include "IMATRIX.h"

/// liblibra 
namespace liblibra{

using namespace std;


/// liblinalg namespace
namespace liblinalg{



class CMATRIX : public base_matrix< complex<double> >{
/**
  The class representing an arbitrary-sized complex valued matrices
*/

public:

  ///========= Constructors and destructors ===============
  CMATRIX() : base_matrix< complex<double> >() { }
  CMATRIX(int i, int j) : base_matrix< complex<double> >(i,j) { }
/*
  CMATRIX(const CMATRIX& ob) : base_matrix< complex<double> >(ob) {   }

*/
  CMATRIX(const CMATRIX& ob) {
    n_rows = ob.n_rows;  ///< The number of rows
    n_cols = ob.n_cols;  ///< The number of colomns
    n_elts = ob.n_elts;  ///< The number of elements
    M = new complex<double>[n_elts];
    memcpy(M, ob.M, sizeof(complex<double>)*n_elts);    
  }


  /// Type-specific Constructors
  ///< Create the complex-valued matrix from two tables: one for real, one for imaginary components
  CMATRIX(vector<vector<double> >& re_part,vector<vector<double> >& im_part); 

  ///< Create the complex-valued matrix (with zero imaginary components) from an int-valued matrix
//  CMATRIX(IMATRIX& re_part);  
//  CMATRIX(const IMATRIX& re_part);  

  ///< Create the complex-valued matrix (with zero imaginary components) from a real-valued matrix
  CMATRIX(MATRIX& re_part);  
  CMATRIX(const MATRIX& re_part);  

  ///< Create the complex-valued matrix from two real-valued matrices: one for real, one for imaginary components
  CMATRIX(MATRIX& re_part,MATRIX& im_part); 


 ~CMATRIX(){}


  ///========== Getters and setters ====================
  /// Inherit the base methods

  using base_matrix<complex<double> >::set;
  using base_matrix<complex<double> >::get;


  ///< Sets the indx's emelent of the M array to the input value (real and imaginary components)  
  void set(int indx,double value1, double value2){
    set(indx, complex<double>(value1,value2));
  }

  void set(CMATRIX& value1){ set(value1); }

  ///< Sets the "row","col" matrix emelent of the M array to the input value (real and imaginary components) 
  void set(int row,int col,double value1,double value2){
    set(row, col, complex<double>(value1,value2));
  }


  ///=========== Extractions ======================
//  using base_matrix<complex<double> >::col;
//  using base_matrix<complex<double> >::row;

  ///========= Initialization =====================
  using base_matrix<complex<double> >::diag;
  using base_matrix<complex<double> >::identity;
  using base_matrix<complex<double> >::Init;
  using base_matrix<complex<double> >::InitSquareMatrix;
  using base_matrix<complex<double> >::Init_Unit_Matrix;

  /// For Backward-compatibility
  void load_identity(){ identity(); } ///< Reset the caller matrix to the identity matrix (all real)


  ///========= Operations =====================
  using base_matrix<complex<double> >::add;
  using base_matrix<complex<double> >::scale;
  using base_matrix<complex<double> >::product;
  using base_matrix<complex<double> >::dot_product;


  ///========= Transformation =====================
  using base_matrix<complex<double> >::Transpose;
//  using base_matrix<complex<double> >::T;
  using base_matrix<complex<double> >::swap_cols;
  using base_matrix<complex<double> >::swap_rows;
  using base_matrix<complex<double> >::permute_cols;
  using base_matrix<complex<double> >::permute_rows;
  using base_matrix<complex<double> >::RightRotation;
  using base_matrix<complex<double> >::LeftRotation;
 

  ///========== Return derivative matrices ===========
  CMATRIX T();   ///< Returns the matrix which is transposed w.r.t. the caller matrix
  CMATRIX H();   ///< Returns the matrix which is Hermitian conjugate w.r.t. the caller matrix
  CMATRIX conj();///< Returns the matrix which is complex conjugate w.r.t. the caller matrix

  MATRIX real(); ///< Returns the real component of the matrix
  MATRIX imag(); ///< Returns the imaginary component of the matrix
  void get_components(MATRIX& re_part,MATRIX& im_part); ///< Split the matrix into real and imaginary components 

  CMATRIX col(int i); ///< takes given column and makes it n x 1 CMATRIX
  CMATRIX row(int i); ///< takes given column and makes it n x 1 CMATRIX

  ///================ Matrix properties =====================
  /// Inherited properties

  using base_matrix<complex<double> >::tr;
  using base_matrix<complex<double> >::sum;
  using base_matrix<complex<double> >::sum_col;
  using base_matrix<complex<double> >::sum_row;
  using base_matrix<complex<double> >::prod_col;
  using base_matrix<complex<double> >::prod_row;



  double NonOrtogonality_Measure();
  complex<double> max_elt();
  void FindMaxNondiagonalElement(int& row,int& col,complex<double>& value);
  void max_nondiagonal(int& row,int& col); // Backward-compatibility

  void max_col_elt(int, complex<double>&, int&); ///< Finds the maximal element (in abs. value) and its index in a given column
  void min_col_elt(int, complex<double>&, int&); ///< Finds the maximal element (in abs. value) and its index in a given column
  void max_row_elt(int, complex<double>&, int&); ///< Finds the maximal element (in abs. value) and its index in a given row
  void min_row_elt(int, complex<double>&, int&); ///< Finds the maximal element (in abs. value) and its index in a given row
  boost::python::list max_col_elt(int); ///< Finds the maximal element (in abs. value) and its index in a given column
  boost::python::list min_col_elt(int); ///< Finds the maximal element (in abs. value) and its index in a given column
  boost::python::list max_row_elt(int); ///< Finds the maximal element (in abs. value) and its index in a given row
  boost::python::list min_row_elt(int); ///< Finds the maximal element (in abs. value) and its index in a given row



  ///=================== Matrix IO and preparation ===================  
  using base_matrix<complex<double> >::bin_dump;
  using base_matrix<complex<double> >::bin_load;
  using base_matrix<complex<double> >::show_matrix;
  using base_matrix<complex<double> >::Load_Matrix_From_File;
  using base_matrix<complex<double> >::show_matrix_address;


  ///================ Operator overloads =====================

  using base_matrix<complex<double> >::operator=;
  using base_matrix<complex<double> >::operator+=;
  using base_matrix<complex<double> >::operator-=;
  using base_matrix<complex<double> >::operator*=;
  using base_matrix<complex<double> >::operator/=;


  CMATRIX operator+(const CMATRIX& rhs);
  CMATRIX operator+(int rhs);
  CMATRIX operator+(double rhs);
  CMATRIX operator+(complex<double> rhs);

  CMATRIX operator-(const CMATRIX& rhs);
  CMATRIX operator-(int rhs);
  CMATRIX operator-(double rhs);
  CMATRIX operator-(complex<double> rhs);

  void operator+=(complex<double> f);
  void operator-=(complex<double> f);

  CMATRIX operator*(const CMATRIX& ob);
  friend CMATRIX operator*(const CMATRIX& ob, int f);
  friend CMATRIX operator*(const CMATRIX& ob, double f);
  friend CMATRIX operator*(const CMATRIX& ob, complex<double> f);
  friend CMATRIX operator*(int f, const CMATRIX& ob);
  friend CMATRIX operator*(double f, const CMATRIX& ob);
  friend CMATRIX operator*(complex<double> f, const CMATRIX& ob);

//  friend CMATRIX operator*(const MATRIX& mtx1, const CMATRIX& mtx2);
//  friend CMATRIX operator*(const CMATRIX& mtx1, const MATRIX& mtx2);
//  friend CMATRIX operator*(const IMATRIX& mtx1, const CMATRIX& mtx2);
//  friend CMATRIX operator*(const CMATRIX& mtx1, const IMATRIX& mtx2);


  void operator*=(complex<double> f);
  void operator/=(complex<double> f);
  CMATRIX operator/(int f);
  CMATRIX operator/(double f);
  CMATRIX operator/(complex<double> f);


                 
};


typedef std::vector<CMATRIX> CMATRIXList;  ///< Data type holding a list of arbitrary-size complex-valued matrices
typedef std::vector<vector<CMATRIX> > CMATRIXMap; ///< Data type for storing the table (grid) of the arbitrary-size complex-valued matrices

}//namespace liblinalg
}// liblibra

#endif // CMATRIX_H
