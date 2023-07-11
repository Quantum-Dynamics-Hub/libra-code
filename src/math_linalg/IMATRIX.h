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
  \file IMATRIX.h
  \brief The file describes the IMATRIX class for representing arbitrary size integer-valued matrices    

  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  NOTE:  the functions that return IMATRIX should not be relpaced by the base template functions

  for some reason, the base_matrix<int> is not recognized as IMATRIX so 

  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*/


#ifndef IMATRIX_H
#define IMATRIX_H

#include "base_matrix.h"
//#include "MATRIX.h"

/// liblibra 
namespace liblibra{

using namespace std;


/// liblinalg namespace
namespace liblinalg{



class IMATRIX : public base_matrix< int >{
/**
  The class representing an arbitrary-sized integer-valued matrices
*/

public:

  ///========= Constructors and destructors ===============
  IMATRIX() : base_matrix< int >() { }
  IMATRIX(int i, int j) : base_matrix< int >(i,j) { }
/*
  CMATRIX(const CMATRIX& ob) : base_matrix< complex<double> >(ob) {   }

*/
  IMATRIX(const IMATRIX& ob) {
    n_rows = ob.n_rows;  ///< The number of rows
    n_cols = ob.n_cols;  ///< The number of colomns
    n_elts = ob.n_elts;  ///< The number of elements
    M = new int[n_elts];
    memcpy(M, ob.M, sizeof(int)*n_elts);    
  }


  /// Type-specific Constructors
  IMATRIX(vector< vector<int> >& mtx);  

  ///< Create the int-valued matrix from a real-valued matrix by converting the numbers to ints
//  IMATRIX(MATRIX& mtx);  

 ~IMATRIX(){}


  ///========== Getters and setters ====================
  /// Inherit the base methods

  using base_matrix< int >::set;
  using base_matrix< int >::get;


  ///< Sets the indx's emelent of the M array to the input value (real and imaginary components)  
/*
  void set(int indx, int value){
    set(indx, value);
  }

  ///< Sets the "row","col" matrix emelent of the M array to the input value (real and imaginary components) 
  void set(int row, int col, int value){
    set(row, col, value);
  }


  void set(int indx, double value){
    set(indx, (int)value);
  }

  void set(int row, int col, double value){
    set(row, col, (int)value);
  }
*/

  ///=========== Extractions ======================
//  using base_matrix<int>::col;
//  using base_matrix<int>::row;


  ///========= Initialization =====================
  using base_matrix< int >::diag;
  using base_matrix< int >::identity;
  using base_matrix< int >::Init;
  using base_matrix< int >::InitSquareMatrix;
  using base_matrix< int >::Init_Unit_Matrix;

  /// For Backward-compatibility
  void load_identity(){ identity(); } ///< Reset the caller matrix to the identity matrix (all real)


  ///========= Operations =====================
  using base_matrix< int >::add;
  using base_matrix< int >::scale;
  using base_matrix< int >::product;
  using base_matrix< int >::dot_product;


  ///========= Transformation =====================
  using base_matrix< int >::Transpose;
//  using base_matrix< int >::T;
  using base_matrix< int >::swap_cols;
  using base_matrix< int >::swap_rows;
  using base_matrix< int >::permute_cols;
  using base_matrix< int >::permute_rows;
 

  ///========== Return derivative matrices ===========
  IMATRIX T();   ///< Returns the matrix which is transposed w.r.t. the caller matrix
  IMATRIX col(int i); ///< takes given column and makes it n x 1 IMATRIX
  IMATRIX row(int i); ///< takes given column and makes it n x 1 IMATRIX


  ///================ Matrix properties =====================
  /// Inherited properties

  using base_matrix< int >::tr;
  using base_matrix< int >::sum;
  using base_matrix< int >::sum_col;
  using base_matrix< int >::sum_row;
  using base_matrix< int >::prod_col;
  using base_matrix< int >::prod_row;


  int max_elt();
  void FindMaxNondiagonalElement(int& row, int& col, int& value);
  void max_nondiagonal(int& row,int& col); // Backward-compatibility

  void max_col_elt(int, int&, int&); ///< Finds the maximal element (in abs. value) and its index in a given column
  void min_col_elt(int, int&, int&); ///< Finds the maximal element (in abs. value) and its index in a given column
  void max_row_elt(int, int&, int&); ///< Finds the maximal element (in abs. value) and its index in a given row
  void min_row_elt(int, int&, int&); ///< Finds the maximal element (in abs. value) and its index in a given row
  boost::python::list max_col_elt(int); ///< Finds the maximal element (in abs. value) and its index in a given column
  boost::python::list min_col_elt(int); ///< Finds the maximal element (in abs. value) and its index in a given column
  boost::python::list max_row_elt(int); ///< Finds the maximal element (in abs. value) and its index in a given row
  boost::python::list min_row_elt(int); ///< Finds the maximal element (in abs. value) and its index in a given row



  ///=================== Matrix IO and preparation ===================  
  using base_matrix< int >::bin_dump;
  using base_matrix< int >::bin_load;
  using base_matrix< int >::show_matrix;
  using base_matrix< int >::Load_Matrix_From_File;
  using base_matrix< int >::show_matrix_address;


  ///================ Operator overloads =====================

  using base_matrix< int >::operator=;
  using base_matrix< int >::operator+=;
  using base_matrix< int >::operator-=;
  using base_matrix< int >::operator*=;
  using base_matrix< int >::operator/=;


  IMATRIX operator+(const IMATRIX& rhs);
  IMATRIX operator+(int rhs);

  IMATRIX operator-(const IMATRIX& rhs);
  IMATRIX operator-(int rhs);

/*
  void operator+=(int f);
  void operator-=(int f);
*/
  IMATRIX operator*(const IMATRIX& ob);
  friend IMATRIX operator*(const IMATRIX& ob, int f);
  friend IMATRIX operator*(int f, const IMATRIX& ob);

/*
  void operator*=(int f);
  void operator/=(int f);
*/
  IMATRIX operator/(int f);

                 
};



typedef std::vector<IMATRIX> IMATRIXList;  ///< Data type holding a list of arbitrary-size complex-valued matrices
typedef std::vector<vector<IMATRIX> > IMATRIXMap; ///< Data type for storing the table (grid) of the arbitrary-size complex-valued matrices

}//namespace liblinalg
}// liblibra

#endif // IMATRIX_H
