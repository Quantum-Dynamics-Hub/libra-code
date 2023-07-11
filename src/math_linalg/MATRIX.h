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
  \file MATRIX.h
  \brief The file describes the MATRIX class for representing arbitrary real-valued matrices

  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  NOTE:  the functions that return MATRIX should not be relpaced by the base template functions

  for some reason, the base_matrix<double> is not recognized as MATRIX so 

  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
    
*/

#ifndef MATRIX_H
#define MATRIX_H

#include "base_matrix.h"
#include "VECTOR.h"

/// liblibra namespace
namespace liblibra{

/// liblinalg namespace
namespace liblinalg{

//============================================================
// Forward declared dependencies
//class VECTOR;
//class MATRIX;



class MATRIX : public base_matrix<double>{

  public:

//  using base_matrix<double>::M;

  ///========= Constructors and destructors ===============
  MATRIX() : base_matrix<double>(){
//    cout<<"In MATRIX constructor 1\n";
  }
  MATRIX(int i, int j) : base_matrix<double>(i,j){
//    cout<<"In MATRIX constructor 2\n";
  }

/*
  MATRIX(const MATRIX& ob){ 
    n_rows = ob.n_rows;  ///< The number of rows
    n_cols = ob.n_cols;  ///< The number of colomns
    n_elts = ob.n_elts;  ///< The number of elements
    M = new double[n_elts];
    memcpy(M, ob.M, sizeof(double)*n_elts);    

    cout<<"In MATRIX cctor\n";
  }
*/

//  MATRIX(boost::python::tuple ob);

  MATRIX(const MATRIX& ob){
    n_rows = ob.n_rows;  ///< The number of rows
    n_cols = ob.n_cols;  ///< The number of colomns
    n_elts = ob.n_elts;  ///< The number of elements
    M = new double[n_elts];
    memcpy(M, ob.M, sizeof(double)*n_elts);    
  }


  MATRIX(const VECTOR& u1, const VECTOR& u2, const VECTOR& u3);

  ~MATRIX(){  /* cout<<"In MATRIX destructor\n";*/ }

  ///========== Getters and setters ====================
  /// Inherit the base methods - this is enough for this class

  using base_matrix<double>::set;
  using base_matrix<double>::get;


  ///=========== Extractions ======================
//  using base_matrix<double >::col;
//  using base_matrix<double >::row;


  ///========= Initialization =====================
  using base_matrix<double>::diag;
  using base_matrix<double>::identity;
  using base_matrix<double>::Init;
  using base_matrix<double>::InitSquareMatrix;
  using base_matrix<double>::Init_Unit_Matrix;

  /// Type-specific methods
  void init(VECTOR& u1, VECTOR& u2, VECTOR& u3);
  void init(const VECTOR& u1, const VECTOR& u2, const VECTOR& u3);


  ///========= Operations =====================
  using base_matrix<double>::add;
  using base_matrix<double>::scale;
  using base_matrix<double>::product;
  using base_matrix<double>::dot_product;

  ///========= Transformation =====================
  using base_matrix<double>::Transpose;
//  using base_matrix<double>::T;
  using base_matrix<double>::swap_cols;
  using base_matrix<double>::swap_rows;
  using base_matrix<double>::permute_cols;
  using base_matrix<double>::permute_rows;
  using base_matrix<double>::RightRotation;
  using base_matrix<double>::LeftRotation;
 

  ///========== Return derivative matrices ===========

  MATRIX T();   ///< Returns the matrix which is transposed w.r.t. the caller matrix
  MATRIX col(int i); ///< takes given column and makes it n x 1 CMATRIX
  MATRIX row(int i); ///< takes given column and makes it n x 1 CMATRIX

  ///================ Matrix properties =====================
  /// Inherited properties

  using base_matrix<double >::tr;
  using base_matrix<double >::sum;
  using base_matrix<double >::sum_col;
  using base_matrix<double >::sum_row;
  using base_matrix<double >::prod_col;
  using base_matrix<double >::prod_row;


  double NonOrtogonality_Measure();
  double max_elt();
  void FindMaxNondiagonalElement(int& row,int& col,double& value);
  void max_col_elt(int I, double& val, int& max_elt_indx);
  void min_col_elt(int I, double& val, int& min_elt_indx);
  void max_row_elt(int I, double& val, int& max_elt_indx);
  void min_row_elt(int I, double& val, int& min_elt_indx);
  boost::python::list max_col_elt(int I);
  boost::python::list min_col_elt(int I);
  boost::python::list max_row_elt(int I);
  boost::python::list min_row_elt(int I);



  ///=================== Matrix IO and preparation ===================  
  using base_matrix<double >::bin_dump;
  using base_matrix<double >::bin_load;
  using base_matrix<double >::show_matrix;
  using base_matrix<double >::Load_Matrix_From_File;
  using base_matrix<double >::show_matrix_address;




  ///================ Operator overloads =====================

/*
void operator=(const base_matrix<T1>& ob)
void operator=(int f)
void operator=(double f)
void operator=(complex<double> f)
void operator+=(const base_matrix<T1>& ob)
void operator+=(int f)
void operator+=(double f)
void operator-=(const base_matrix<T1>& ob)
void operator-=(int f)
void operator-=(double f)
void operator*=(int f)
void operator*=(double f)
void operator/=(int f)
void operator/=(double f)
*/



  using base_matrix<double>::operator=;  
  using base_matrix<double>::operator+=;
  using base_matrix<double>::operator-=;
  using base_matrix<double>::operator*=;
  using base_matrix<double>::operator/=;


  MATRIX operator-() const;   ///< Negation  

  MATRIX operator+(const MATRIX& rhs);  ///< Addition operator
  MATRIX operator+(int rhs);     ///< Addition operator
  MATRIX operator+(double rhs);  ///< Addition operator

  MATRIX operator-(const MATRIX& rhs); ///< Subtraction operator by a matrix type
  MATRIX operator-(int rhs);     ///< Subtraction operator 
  MATRIX operator-(double rhs);  ///< Subtraction operator 

  MATRIX operator*(const MATRIX& ob); ///< Multiplication operator
  friend MATRIX operator*(const MATRIX& ob, int f);
  friend MATRIX operator*(int f, const MATRIX& ob);

  friend MATRIX operator*(const MATRIX& ob, double f);
  friend MATRIX operator*(double f, const MATRIX& ob);
 
  ///< Multiplication of vector and matrix:  matrix * vector
  friend VECTOR operator*(const MATRIX& m,  const VECTOR& v);   

 
  MATRIX operator/(int f); ///< Division by a scalar
  MATRIX operator/(double f); ///< Division by a scalar

  ///< Tensor product of two 3D vectors to yield a 3x3 matrix
  friend MATRIX operator^(const VECTOR& v1, const VECTOR& v2); 
  

  ///=================== Misc matrix methods  ===================  

  void tensor_product(const VECTOR& v1, const VECTOR& v2);
  void get_vectors(VECTOR& u1,VECTOR& u2,VECTOR& u3);

  void skew(const VECTOR& v);
  void skew1(const VECTOR& v);
  void Rotation(const VECTOR& u);
  void Rx(double phi);
  void Ry(double phi);
  void Rz(double phi);



};


typedef std::vector<MATRIX> MATRIXList;  ///< Data type for storing a vector of MATRIX objects
typedef std::vector<vector<MATRIX> > MATRIXMap;  ///< Data type for storing a table (grid) of MATRIX objects


//-------- IO functions --------
void set_value(int& defined, MATRIX& value, boost::python::object obj, std::string attrName);
void save(boost::property_tree::ptree& pt,std::string path,MATRIX& vt);
void save(boost::property_tree::ptree& pt,std::string path, char path_separator, MATRIX& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt);
void save(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<MATRIX>& vt);

void load(boost::property_tree::ptree& pt,std::string path,MATRIX& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path, char path_separator, MATRIX& vt, int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path, char path_separator,vector<MATRIX>& vt,int& status);


}// namespace liblinalg
}// namespace liblibra

#endif // MATRIX.h


