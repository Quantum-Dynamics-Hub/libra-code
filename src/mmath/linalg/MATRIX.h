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
  \file MATRIX.h
  \brief The file describes the MATRIX class for representing arbitrary real-valued matrices
    
*/

#ifndef MATRIX_H
#define MATRIX_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vector>
using namespace std;

#include "../../io/libio.h"
using namespace libio;

/// libmmath namespace
namespace libmmath{

/// liblinalg namespace
namespace liblinalg{

//============================================================
// Forward declared dependencies
class VECTOR;
class MATRIX;

//===============================================================================

class MATRIX{
/**
  The class for representing an arbitrary real-valued matrices
*/
    
  double *Orthogonalization_Matrix;
    
public:
  // Data; 
  int num_of_rows;  ///< The number of the matrix rows
  int num_of_cols;  ///< The number of the matrix coloumns
  int num_of_elems; ///< The number of the matrix elements = num_of_rows x num_of_cols
  double *M;        ///< The array storing matrix elements
    
  inline void set(int indx,double value){    M[indx] = value; }  ///< Sets the indx's emelent of the M array to the input value 
  inline double get(int indx){    return M[indx];  }   ///< Return the matrix element with the 1-D index "indx"
  inline void set(int row,int col,double value){   M[row*num_of_cols+col] = value;  }  ///< Set the "row", "col" matrix element to the value "value"
  inline double get(int row,int col){  return M[row*num_of_cols+col]; }  ///< Return the "row", "col" matrix element

          
  double NonOrtMeasure;  ///< Non-orthogonality measure - used internally in Jacobi rotation algorithm
  // Constructors;
  MATRIX(){
  /** Default constructor */

  MATRIX_PRECISION=8;
     MATRIX_WIDTH=15;
     num_of_rows=3;
     num_of_cols=3;
     num_of_elems=9;    
     M = new double[9];
     for(int i=0;i<9;i++) { M[i]=0.0;}
  } 
  
  MATRIX(int m,int n);
  MATRIX(const MATRIX& ob);   // Copy constructor;
  MATRIX(const VECTOR&, const VECTOR&, const VECTOR&);
 ~MATRIX();

  // Initializations;
  void Init(double x);            ///< Set all matrix elements to x
  void InitSquareMatrix(int dim,double x);  ///< initialize dim x dim matrix to be x * I, where I - identity matrix
  bool Init_Unit_Matrix(double x);///< Set the matrix (already allocated) to be x * I, where I - identity matrix
  bool Load_Matrix_From_File(char * FileName);  ///< Initialize matrix by reading the data from file (text format)
  void init(VECTOR&, VECTOR&, VECTOR&); ///< Initialize 3x3 matrix with 3 input elements being its columns
  void init(const VECTOR&, const VECTOR&, const VECTOR&);  ///< Initialize 3x3 matrix with 3 input elements being its columns
  void init(MATRIX&); ///< Initialize the matrix to be a copy of an imput matrix
  void init(const MATRIX&);  ///< Initialize the matrix to be a copy of an imput matrix

  // Rotation matrices
  void Rotation(const VECTOR&); ///< Create a rotation (3x3) matrix for rotation around direction given by the direction
                                ///< of the argument vector, and on the amount given by the norm of the vector

  // Rotation matrixes around corresponding axes
  void Rx(double);  ///< The 3x3 matrix for rotation around x axis
  void Ry(double);  ///< The 3x3 matrix for rotation around y axis
  void Rz(double);  ///< The 3x3 matrix for rotation around z axis

  // Manipulations with colomns, rows and elements;               
  int show_num_of_rows() { return num_of_rows;}  ///< Return the number of rows
  int show_num_of_cols() { return num_of_cols;}  ///< Return the number of coloumns
  int show_num_of_elems(){ return num_of_elems;}  ///< Return the number of elements (rows x coloumns)

  void Add_To_Element(int row,int col,double x) {int n=row*num_of_cols+col; *(M+n)+=x;}  ///< Increment the "row", "col" matrix element by "x"
  void FindMaxNondiagonalElement(int& rw,int& cl,double& value); ///< An auxiliary function - finds the maximal off-diagonal element and 
                                                                 ///< returns its position and magnitude into the input references

  // Tehniques for ortogonalization and decomposition;
  double NonOrtogonality_Measure();    ///< Compute a certain non-orthogonality measure 
  void Show_Orthogonalization_Matrix();  ///< Print the orthogonalization matrix (available after Jacobi rotation)
  void Show_Orthogonalization_Matrix(char *Output_File);  ///< Print the orthogonalization matrix (available after Jacobi rotation) to the output file
  void Delete_Orthogonalization_Matrix();  ///< De-allocate memory used by the orthogonalizaion matrix
  void RightRotation(int i,int j,double sine,double cosine); 
  void LeftRotation(int i,int j,double sine,double cosine);
  void Ortogonalization(double Eps);
  void RL_Decomposition(MATRIX *L,MATRIX *R);

  // Basic matrix operations;
  void Transpose();   ///< Transposes this matrix. That is the caller object is changed
  MATRIX T();         ///< Returns the matrix which is transposed w.r.t. the caller object
  void Inverse(MATRIX* INV);  ///< Compute the matrix which is inverse w.r.t. to the caller matrix and stores it in the matrix pointed to be the input argument
  void Inverse(MATRIX& INV);  ///< Compute the matrix which is inverse w.r.t. to the caller matrix and stores it in the matrix referenced by the input argument
//  void tensor_product(VECTOR v1,VECTOR v2);
  void tensor_product(VECTOR v1,VECTOR v2); ///< Compute a tensor (outer) product of 2 3D vectors and store the result in the 3x3 matrix

  double dot_product(MATRIX&);  ///< Compute the dot product of the caller matrix and the input argument
  double dot_product(MATRIX*);  ///< Compute the dot product of the caller matrix and the input argument

  MATRIX col(int); ///< takes given column and makes it n x 1 MATRIX
  MATRIX row(int); ///< takes given row and makes it 1 x n MATRIX

/*
  inline void tensor_product(VECTOR& v1,VECTOR& v2){
    M[0] = v1.x*v2.x;   M[1] = v1.x*v2.y;  M[2] = v1.x*v2.z;
    M[3] = v1.y*v2.x;   M[4] = v1.y*v2.y;  M[5] = v1.y*v2.z;
    M[6] = v1.z*v2.x;   M[7] = v1.z*v2.y;  M[8] = v1.z*v2.z;
  }
*/
  double& operator[](const int indx); ///< This returns M[indx], where indx - the input argument
  MATRIX operator-();                 ///< Negation of the caller matrix. This will change the caller matrix.
  MATRIX operator*(const MATRIX&);    ///< Matrix multiplication
  MATRIX operator+(const MATRIX&);    ///< Matrix addition
  MATRIX operator-(const MATRIX&);    ///< Matrix subtraction
  MATRIX& operator+=(const MATRIX&);  ///< Matrix increment by a matrix
  MATRIX& operator-=(const MATRIX&);  ///< Matrix decrement by a matrix
  MATRIX& operator*=(double);         ///< Matrix scaling by a scalar
  MATRIX operator/(double num);       ///< Matrix division by a scalar
  MATRIX operator=(const MATRIX&);    ///< Copying one matrix into the other one
  MATRIX operator=(double num);       ///< Setting all matrix elements to a given scal

  friend int operator == (const MATRIX& m1, const MATRIX& m2);  ///< Are matrices equal?;
  friend int operator != (const MATRIX& m1, const MATRIX& m2);  ///< Are matrices not equal?;

  friend double operator%(MATRIX& m1, MATRIX& m2);  ///< Scalar product of matrices;
  friend MATRIX operator^(VECTOR& v1,VECTOR& v2);   ///< Tensor product of two 3D vectors to yield a 3x3 matrix
  friend VECTOR operator*(const MATRIX& m,  const VECTOR& v);   ///< Multiplication of vector and matrix:  matrix * vector
//  friend VECTOR operator*(const VECTOR& v,  const MATRIX& m);   // Multiplication of vector and matrix; 
  friend MATRIX operator*(const double& f,  const MATRIX& m1);  ///< Multiplication of matrix and double:  scalar * matrix
  friend MATRIX operator*(const MATRIX &m1, const double  &f);  ///< Multiplication of matrix and double:  matrix * scalar


  // Input & output functions;
  void show_matrix();  ///< Print out the matrix
  void show_matrix(char *Output_File);  ///< Print the matrix to the file
  friend ostream &operator<<(ostream &strm,MATRIX ob);  ///< Output the matrix to a stream
  friend istream& operator>>(istream& strm,MATRIX &ob); ///< Input a matrix from a stream
  void Delete_Matrix();  ///< Free memory allocated to store matrix elements (almost destructor)

  // Additional methods for matrixes.
  void get_vectors(VECTOR&,VECTOR&,VECTOR&);  ///< For 3x3 matrix - take the columns of the matrix as the vectors
  void skew(VECTOR);   /// Generate a 3x3 skew-symmetric matrix from a 3D vector
  void skew1(VECTOR);  /// Generate a 4x4 skew-symmetric matrix from a 3D vector - for quaternionic representation
  void exp(MATRIX&);   /// Comute the exponential of the input matrix - result is stored in the caller object. Use Taylor series
//  void exp(const MATRIX&);
  void JACOBY_EIGEN(MATRIX&, MATRIX&);  ///< Jacoby-rotations based matrix diagonalization procedure for real symmetric matrices
                                        ///< - returns eigenvalues and eigenvectors
  void JACOBY_EIGEN(MATRIX&, MATRIX&,double);///< Jacoby-rotations based matrix diagonalization procedure for real symmetric matrices
                                             ///< - returns eigenvalues and eigenvectors. This version allows precision control.

  // Properties;
  int MATRIX_PRECISION;  ///< The precision of the matrix elements printing (how many positions)
  int MATRIX_WIDTH;      ///< Width for each printed matrix element (includes number and whitespaces around)
  double Determinant();  ///< Compute matrix determinant: so far is only available for 2x2 and 3x3 matrices
  double tr();          ///< Compute the trace of the matrix
  double max_elt();     ///< Compute the maximal in magnitude element

  // Binary output/input
  void bin_dump(std::string filename);  ///< Store the matrix into a file in a binary format
  void bin_load(std::string filename);  ///< Load the matrix from a binary file (created with bin_dump)
};


void pop_submatrix(MATRIX*,MATRIX*,vector<int>&);  ///< Take a smaller sub-matrix from a bigger matrix according to the specified pattern
void push_submatrix(MATRIX*,MATRIX*,vector<int>&); ///< Push a smaller matrix into a bigger one according to the specified pattern

typedef std::vector<MATRIX> MATRIXList;  ///< Data type for storing a vector of MATRIX objects
typedef std::vector<vector<MATRIX> > MATRIXMap;  ///< Data type for storing a table (grid) of MATRIX objects

//-------- IO functions --------
void set_value(int& defined, MATRIX& value, boost::python::object obj, std::string attrName);
void save(boost::property_tree::ptree& pt,std::string path,MATRIX& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt);
void load(boost::property_tree::ptree& pt,std::string path,MATRIX& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt,int& status);



}// namespace liblinalg
}// namespace libmmath

#endif // MATRIX.h


