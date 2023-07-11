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
  \file MATRIX.cpp
  \brief The file implements the MATRIX class for representing arbitrary real-valued matrices
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <stdio.h>
#include <string.h>
#endif

#include "MATRIX.h"
#include "VECTOR.h"


// ========================= Matrices ================================
// ------------------------- Constructors ----------------------------

/// liblibra namespace
namespace liblibra{

/// liblinalg namespace
namespace liblinalg{

/*
MATRIX::MATRIX(boost::python::tuple ob){

    n_rows = extract<int>( ob[0] );
    n_cols = extract<int>( ob[1] );
    n_elts = n_rows * n_cols;

    M = new double[n_elts];

    for(int i=0;i<n_elts;i++){  
      M[i] = extract<double>(ob[i+2]);
    }

}
*/


MATRIX::MATRIX(const VECTOR& u1, const VECTOR& u2, const VECTOR& u3){
/** The constructor of a 3 x 3 matrix from 3 input vectors  */

  n_rows=3;   n_cols=3;   n_elts=9;
  M=new double[n_elts];

  M[0] = u1.x;  M[1] = u2.x;  M[2] = u3.x;
  M[3] = u1.y;  M[4] = u2.y;  M[5] = u3.y;
  M[6] = u1.z;  M[7] = u2.z;  M[8] = u3.z;

} 



void MATRIX::init(VECTOR& u1, VECTOR& u2, VECTOR& u3){
/** Initialize the 3 x 3 matrix by setting each of its columns to 
the corresponding given input vectors
*/
  if(n_cols==3 && n_rows==3){
    M[0] = u1.x;  M[1] = u2.x;  M[2] = u3.x;
    M[3] = u1.y;  M[4] = u2.y;  M[5] = u3.y;
    M[6] = u1.z;  M[7] = u2.z;  M[8] = u3.z;
  }
  else{ 
    cout<<"Error in MATRIX::init function.\n";
    cout<<"This function works only for 3 x 3 matrices\n";
    exit(0);
  }

}

void MATRIX::init(const VECTOR& u1, const VECTOR& u2, const VECTOR& u3){
/** Initialize the 3 x 3 matrix by setting each of its columns to 
the corresponding given input vectors
*/
  if(n_cols==3 && n_rows==3){
    M[0] = u1.x;  M[1] = u2.x;  M[2] = u3.x;
    M[3] = u1.y;  M[4] = u2.y;  M[5] = u3.y;
    M[6] = u1.z;  M[7] = u2.z;  M[8] = u3.z;
  }
  else{ 
    cout<<"Error in MATRIX::init function.\n";
    cout<<"This function works only for 3 x 3 matrices\n";
    exit(0);
  }

}




MATRIX MATRIX::T(){   
// Returns the matrix which is transposed w.r.t. the caller matrix  

  MATRIX res(*this); res.Transpose();
  return res;    
}


MATRIX MATRIX::col(int i){
// takes given column and makes it n x 1 CMATRIX  

  MATRIX tmp(n_rows,1);
  for(int j=0;j<n_rows;j++){ tmp.M[j] = M[j*n_cols+i]; }
  return tmp;
}

MATRIX MATRIX::row(int i){ 
// takes given row and makes it 1 x n CMATRIX  

  MATRIX tmp(1,n_cols);
  for(int j=0;j<n_cols;j++){ tmp.M[j] = M[i*n_cols+j]; }
  return tmp;
}


double MATRIX::NonOrtogonality_Measure(){
/** The sum of the scalar products of all pairs of distinct columns
*/
  double sum, aa; sum = 0.0; 

  for(int i=0;i<n_cols-1;i++){
    for(int j=i+1;j<n_cols;j++){

      aa = 0.0;
      for(int k=0;k<n_rows;k++){ aa += M[k*n_cols+i] * M[k*n_cols+j];  }
      sum += aa*aa;
    }
  }
  return sum;
}

double MATRIX::max_elt(){
/** Finds the maximal (in absolute value) element and its position  */

  double x = fabs(M[0]);
  double y;
  for(int i=0;i<n_elts;i++){  y = fabs(M[i]); if(y>=x){ x = y; } }
  return x;

}

void MATRIX::FindMaxNondiagonalElement(int& row,int& col,double& value){
/** Finds the maximal (in absolute value) off-diagonal element and its position  */

  int k=0;
  double elem, max_elem;
  value = M[1]; max_elem = fabs(value); row = 0 ; col = 1;

  for(int rw=0;rw<n_rows;rw++){
    for(int cl=rw+1;cl<n_cols;cl++){
      k = rw*n_cols + cl;
      elem = fabs(M[k]);
      if(elem>max_elem) {max_elem = elem; value = M[k]; col = cl; row = rw;}
    }
  }
}


void MATRIX::max_col_elt(int I, double& val, int& max_elt_indx){
/** Finds the maximal element (in abs. value) and its index in a given column */

  val = M[0*n_cols+I];
  max_elt_indx = 0;
  double max_norm = fabs(val);

  for(int row=0; row<n_rows; row++){

    double nrm = fabs(M[row*n_cols+I]);
    if(nrm>max_norm){  
      max_norm = nrm;   
      val = M[row*n_cols+I];
      max_elt_indx = row;
    }
  }

}


void MATRIX::min_col_elt(int I, double& val, int& min_elt_indx){
/** Finds the minimal element (in abs. value) and its index in a given column  */

  val = M[0*n_cols+I];
  min_elt_indx = 0;
  double min_norm = fabs(val);

  for(int row=0; row<n_rows; row++){

    double nrm = fabs(M[row*n_cols+I]);
    if(nrm < min_norm){  
      min_norm = nrm;   
      val = M[row*n_cols+I];
      min_elt_indx = row;
    }
  }// for row

}


void MATRIX::max_row_elt(int I, double& val, int& max_elt_indx){
/** Finds the maximal element (in abs. value) and its index in a given row  */

  val = M[I*n_cols+0];
  max_elt_indx = 0;
  double max_norm = fabs(val);

  for(int col=0; col<n_cols; col++){

    double nrm = fabs(M[I*n_cols+col]);
    if(nrm>max_norm){  
      max_norm = nrm;   
      val = M[I*n_cols+col];
      max_elt_indx = col;
    }

  }// for col

}


void MATRIX::min_row_elt(int I, double& val, int& min_elt_indx){
/** Finds the minimal element (in abs. value) and its index in a given row  */

  val = M[I*n_cols+0];
  min_elt_indx = 0;
  double min_norm = fabs(val);

  for(int col=0; col<n_cols; col++){

    double nrm = fabs(M[I*n_cols+col]);
    if(nrm<min_norm){  
      min_norm = nrm;   
      val = M[I*n_cols+col];
      min_elt_indx = col;
    }

  }// for col
}



boost::python::list MATRIX::max_col_elt(int I){
/** Finds the maximal element (in abs. value) and its index in a given column  */

  double val;
  int indx;

  this->max_col_elt(I, val, indx);

  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;
}


boost::python::list MATRIX::min_col_elt(int I){
/** Finds the maximal element (in abs. value) and its index in a given column  */

  double val;
  int indx;

  this->min_col_elt(I, val, indx);

  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;

}


boost::python::list MATRIX::max_row_elt(int I){
/** Finds the maximal element (in abs. value) and its index in a given row  */

  double val;
  int indx;

  this->max_row_elt(I, val, indx);

  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;

}


boost::python::list MATRIX::min_row_elt(int I){
/** Finds the maximal element (in abs. value) and its index in a given row  */

  double val;
  int indx;

  this->min_row_elt(I, val, indx);

  boost::python::list res;

  res.append(indx);
  res.append(val);

  return res;

}




MATRIX MATRIX::operator-() const{
  MATRIX tmp(n_rows,n_cols);

  for(int i=0;i<n_elts;i++){ tmp.M[i] = -M[i]; }
  return tmp;
}


MATRIX MATRIX::operator+(const MATRIX& rhs){ 
  MATRIX res(*this);  res += rhs;
  return res;
}


MATRIX MATRIX::operator+(int rhs){ 
  MATRIX res(*this);  res += rhs;
  return res;
}


MATRIX MATRIX::operator+(double rhs){ 
  MATRIX res(*this);  res += rhs;
  return res;
}


///< Subtraction operator by a matrix type
MATRIX MATRIX::operator-(const MATRIX& rhs){ 
  MATRIX res(*this);  res -= rhs;
  return res;
}

///< Subtraction by a numeric type
MATRIX MATRIX::operator-(int rhs){ 
  MATRIX res(*this);  res -= rhs;
  return res;
}

MATRIX MATRIX::operator-(double rhs){ 
  MATRIX res(*this);  res -= rhs;
  return res;
}



MATRIX MATRIX::operator*(const MATRIX& ob){
  MATRIX res(n_rows, ob.n_cols);  res.product(*this, ob);
  return res;
}


MATRIX operator*(const MATRIX& ob, int f){  
  MATRIX res(ob);  res *= f;
  return res;
}

MATRIX operator*(const MATRIX& ob, double f){  
  MATRIX res(ob);  res *= f;
  return res;
}

MATRIX operator*(int f, const MATRIX& ob){  
  MATRIX res(ob);   res *= f;
  return res;
}

MATRIX operator*(double f, const MATRIX& ob){  
  MATRIX res(ob);   res *= f;
  return res;
}


MATRIX MATRIX::operator/(int f){  
  MATRIX ob(*this); ob /= f;
  return ob;
}

MATRIX MATRIX::operator/(double f){  
  MATRIX ob(*this); ob /= f;
  return ob;
}




void MATRIX::tensor_product(const VECTOR& v1, const VECTOR& v2){

  M[0] = v1.x*v2.x;   M[1] = v1.x*v2.y;  M[2] = v1.x*v2.z;
  M[3] = v1.y*v2.x;   M[4] = v1.y*v2.y;  M[5] = v1.y*v2.z;
  M[6] = v1.z*v2.x;   M[7] = v1.z*v2.y;  M[8] = v1.z*v2.z;
}

void MATRIX::get_vectors(VECTOR& u1,VECTOR& u2,VECTOR& u3){
/** Extract vectors from the 3 x 3 matrix   */

  if(n_cols==3 && n_rows==3){
    u1.x = M[0];  u2.x = M[1];  u3.x = M[2];
    u1.y = M[3];  u2.y = M[4];  u3.y = M[5];
    u1.z = M[6];  u2.z = M[7];  u3.z = M[8];
  }
  else{ 
    cout<<"Error in MATRIX::get_vectors function.\n";
    cout<<"This function works only for 3 x 3 matrices\n";
    exit(0);
  }

}



void MATRIX::skew(const VECTOR& v){
/**
                | 0  -v.z  v.y |
  W = skew(n) = | v.z  0  -v.x |
                |-v.y v.x  0   |
*/

  M[0] = 0.0; M[1] =-v.z; M[2] = v.y;
  M[3] = v.z; M[4] = 0.0; M[5] =-v.x;
  M[6] =-v.y; M[7] = v.x; M[8] = 0.0;

}

void MATRIX::skew1(const VECTOR& v){
/**
                 | 0    -v.x  -v.y  -v.z |
  W = skew(n1) = | v.x    0    v.z  -v.y |
                 | v.y  -v.z    0    v.x |
                 | v.z   v.y  -v.x   0   |
*/

  M[0] = 0.0; M[1] =-v.x; M[2] =-v.y; M[3] =-v.z;
  M[4] = v.x; M[5] = 0.0; M[6] = v.z; M[7] =-v.y;
  M[8] = v.y; M[9] =-v.z; M[10]= 0.0; M[11]= v.x;
  M[12]= v.z; M[13]= v.y; M[14]=-v.x; M[15]= 0.0;

}



void MATRIX::Rotation(const VECTOR& u){
/**
This function initializes the matrix to a rotation matrix
for rotation around the axis given by direction of the vector: (u/|u|)
on amount given by norm of vector u: |u|.
If the vector has a zero length - this will be an identity matrix
*/

  double umod,umod2;

  umod2 = u.x*u.x + u.y*u.y + u.z*u.z;

  if(umod2==0.0){
    M[0] = 1.0;  M[1] = 0.0;  M[2] = 0.0;
    M[3] = 0.0;  M[4] = 1.0;  M[5] = 0.0;
    M[6] = 0.0;  M[7] = 0.0;  M[8] = 1.0;
  }
  else{
    umod = sqrt(umod2);
    VECTOR n(u.x/umod, u.y/umod, u.z/umod);

    double cs,sn;
    cs = (1.0 - cos(umod));
    sn = sin(umod);

    /**
    Here is efficient implementation of Rodrigues formula: 
    M = I + W * sin_psi + W*W*(1.0 - cos_psi)
    where I - identity matrix
                        | 0  -n.z  n.y |
          W = skew(n) = | n.z  0  -n.x |
                        |-n.y n.x  0   |

    */

    double x,y,z,xy,xz,yz,x2,y2,z2;
    x  = n.x * sn;
    y  = n.y * sn;
    z  = n.z * sn;
    x2 = (n.x * n.x - 1.0) * cs;
    y2 = (n.y * n.y - 1.0) * cs;
    z2 = (n.z * n.z - 1.0) * cs;
    xy = n.x * n.y * cs;
    xz = n.x * n.z * cs;
    yz = n.y * n.z * cs;

    M[0] = 1.0 + x2;  M[1] =-z   + xy;  M[2] = y   + xz;
    M[3] = z   + xy;  M[4] = 1.0 + y2;  M[5] =-x   + yz;
    M[6] =-y   + xz;  M[7] = x   + yz;  M[8] = 1.0 + z2;

  }

}

void MATRIX::Rx(double phi){
  double cs = cos(phi);
  double sn = sin(phi);
  M[0] = 1.0; M[1] = 0.0; M[2] = 0.0;
  M[3] = 0.0; M[4] = cs;  M[5] = -sn;
  M[6] = 0.0; M[7] = sn;  M[8] = cs;

}

void MATRIX::Ry(double phi){
  double cs = cos(phi);
  double sn = sin(phi);
  M[0] = cs;  M[1] = 0.0; M[2] = sn;
  M[3] = 0.0; M[4] = 1.0; M[5] = 0.0;
  M[6] = -sn; M[7] = 0.0; M[8] = cs;

}


void MATRIX::Rz(double phi){
  double cs = cos(phi);
  double sn = sin(phi);
  M[0] = cs;  M[1] = -sn; M[2] = 0.0;
  M[3] = sn;  M[4] = cs;  M[5] = 0.0;
  M[6] = 0.0; M[7] = 0.0; M[8] = 1.0;

}





void set_value(int& is_defined, MATRIX& value,boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<MATRIX>(obj.attr(attrName.c_str()));
      is_defined = 1;
  }
}


// ----------- Save --------------
void save(boost::property_tree::ptree& pt,std::string path,MATRIX& vt){
  pt.put(path+".n_rows", vt.n_rows);
  pt.put(path+".n_cols", vt.n_cols);
  pt.put(path+".n_elts", vt.n_elts);

  for(int i=0;i<vt.n_elts;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    pt.put(path+"."+rt,vt.M[i]);
  }
}

void save(boost::property_tree::ptree& pt,std::string path, char path_separator, MATRIX& vt){
  pt.put(path+".n_rows", vt.n_rows);
  pt.put(path+".n_cols", vt.n_cols);
  pt.put(path+".n_elts", vt.n_elts);

  for(int i=0;i<vt.n_elts;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    pt.put(boost::property_tree::ptree::path_type(path, path_separator),vt.M[i]);
  }
}


void save(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+"."+rt,vt[i]);
  }
}

void save(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<MATRIX>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+std::string(1, path_separator)+rt, path_separator,vt[i]);
  }
}


// ----------- Load --------------
void load(boost::property_tree::ptree& pt,std::string path, MATRIX& vt, int& status){ 
  std::cout<<"Sorry: load function for MATRIX object is not defined yet\n"; exit(0); 
}  /// Not useful

MATRIX load(boost::property_tree::ptree& pt,std::string path, int& status){ 

  int st;
  status = 0;

  int n_rows, n_cols, n_elts;

  libio::load(pt,path+".n_rows", n_rows, st); if(st==1) { status=1;}
  libio::load(pt,path+".n_cols", n_cols, st); if(st==1) { status=1;}
  libio::load(pt,path+".n_elts", n_elts, st); if(st==1) { status=1;}

  MATRIX vt(n_rows, n_cols);

  for(int i=0;i<vt.n_elts;i++){

    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;

    libio::load(pt,path+"."+rt, vt.M[i], st); if(st==1) { status=1;}
  }

  return vt;

}


void load(boost::property_tree::ptree& pt,std::string path, char path_separator, MATRIX& vt, int& status){ 
  std::cout<<"Sorry: load function for MATRIX object is not defined yet\n"; exit(0); 
}


void load(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt,int& status){ 

  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      MATRIX x( load(pt,path+"."+v.first, st)); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }

}

void load(boost::property_tree::ptree& pt,std::string path, char path_separator,vector<MATRIX>& vt,int& status){ 
  std::cout<<"Sorry: load function for vector<MATRIX> object is not defined yet\n"; exit(0); 
}







}// namespace liblinalg
}// namespace liblibra


