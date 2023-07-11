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
  \file VECTOR.h
  \brief The file describes the VECTOR class for representing 3D point
    
*/

#ifndef VECTOR_H
#define VECTOR_H

#if defined(USING_PCH)
#include "../pch.h"
#else

#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vector>
#endif 

#include "../io/libio.h"


/// liblibra namespace
namespace liblibra{

using namespace std;
using namespace libio;

/// liblinalg namespace
namespace liblinalg{

//============================================================
// Forward declared dependencies
class VECTOR;
class MATRIX3x3;
class MATRIX;
//==============================================================

class VECTOR{
/**
  This class is designed to hold the 3D point vector - very common in physical applications
*/

  public:
  double x;  ///< x component of the 3D vector
  double y;  ///< y component of the 3D vector
  double z;  ///< z component of the 3D vector

  int is_transposed; ///< transposition state: 1 = true = horizontal notation, 0 = false= vertical notation(default);         

  //---------------------- Constructor ---------------------------------


  VECTOR(){  is_transposed=0; x = 0.0; y = 0.0; z = 0.0; }  ///< Default constructor

  VECTOR(double xval,double yval,double zval){  is_transposed=0; x=xval; y=yval; z=zval; } ///< Constructor with the initialization

  // Copy constructor
  VECTOR(const VECTOR& obj){ is_transposed = obj.is_transposed; x=obj.x; y=obj.y; z=obj.z; }  ///< Copy constructor

  //--------------------- Destructor -----------------------------------


  ~VECTOR(){}

  //------------------- Initializations ---------------------------------


  void init(double val){x=y=z=val;}  ///< Initialize all component of the vecto to the same floating-point value 
  void init(double xval,double yval,double zval){x=xval;y=yval;z=zval; } ///< Initialize the vector components with the values provided
  void init(VECTOR& v){ is_transposed=v.is_transposed; x=v.x; y=v.y; z=v.z; }
  void init(const VECTOR& v){ is_transposed=v.is_transposed; x=v.x; y=v.y; z=v.z; }

  //------------------- Special vectors --------------------------------


  VECTOR unit() const{
  /**
    This function returns the normalized vector - same direction, only the magnitude is unity
  */

    double res; VECTOR U;
    res=sqrt(x*x+y*y+z*z);
    U.x=x/res;
     U.y=y/res;
    U.z=z/res;
    return U;
  } 

  //------------------ Small functions on vectors ----------------------


  inline double length(void) const{  return sqrt(x*x+y*y+z*z); } ///< returns the magnitude of the vector
  inline double length2(void) const{ return (x*x+y*y+z*z); }     ///< returns the square magnitude of the vector
  void normalize(void){
  /**
    Normalize the given vector. This function changes the oribinal object
  */

    double res;
    res=sqrt(x*x+y*y+z*z);
    if(res!=0.0){ x/=res; y/=res; z/=res; }
    else{ x = 0.0; y = 0.0; z = 0.0; }
  }
  inline void cross(const VECTOR &v1, const VECTOR &v2){
  /**
    Compute the cross (vector) product of two vectors, v1 and v2 and store the result in the caller object
  */

    x= v1.y*v2.z-v2.y*v1.z;
    y= v2.x*v1.z-v1.x*v2.z;
    z= v1.x*v2.y-v2.x*v1.y;
  }



  //---------------------- Operators -----------------------------------


  VECTOR operator-() const{  /** Negation operator. Returns the opposite vector */
    VECTOR tmp;
      tmp.x = -x;
      tmp.y = -y;
      tmp.z = -z;
      return tmp;
  }            

  VECTOR operator+(const VECTOR& v) const{  /** Addition operator. Allows v1 + v2 */
    VECTOR tmp;
    tmp.x=x+v.x;
    tmp.y=y+v.y;
    tmp.z=z+v.z;
    return tmp;
  }
  VECTOR operator+(VECTOR& v) const{  /** Another version of addition operator. Allows v1 + v2 */
    VECTOR tmp;
    tmp.x=x+v.x;
    tmp.y=y+v.y;
    tmp.z=z+v.z;
    return tmp;
  }

  VECTOR operator-(const VECTOR& v) const{  /** Subtraction operator. Allows v1 - v2 */
    VECTOR tmp;
    tmp.x=x-v.x;
    tmp.y=y-v.y;
    tmp.z=z-v.z;
    return tmp;
  }
  VECTOR operator-(VECTOR& v) const{  /** Another version of subtraction operator. Allows v1 - v2 */
    VECTOR tmp;
    tmp.x=x-v.x;
    tmp.y=y-v.y;
    tmp.z=z-v.z;
    return tmp;
  }

  VECTOR operator+(const double& v) const{
    /** Addition of vector an scalar operator. Allows v1 + A. This will return a new object. The scalar is added to all components */

    VECTOR tmp;
    tmp.x=x+v;
    tmp.y=y+v;
    tmp.z=z+v;
    return tmp;
  }
  VECTOR operator-(const double& v) const{    
    /** Subtraction of scalar from a vector operator. Allows v1 - A. This will return a new object. The scalar is subtracted from all components */

    VECTOR tmp;
    tmp.x=x-v;
    tmp.y=y-v;
    tmp.z=z-v;
    return tmp;
  }
  VECTOR operator=(const VECTOR &v){  
    /** Assignment operator. Important to return the object, not the reference!!! otherwise crap happens on the Python side   
    Allows v1 = v2 (copy v2 content to v1)
    */

    x=v.x;
    y=v.y;
    z=v.z;
    return *this;
  }
  VECTOR operator=(const double &v){  
    /** Assignment operator. Important to return the object, not the reference!!! otherwise crap happens on the Python side   
    Allows v1 = A (copy A scalar to all components of v1)
    */

    x=v;
    y=v;
    z=v;
    return *this;
  }
  void operator+=(const VECTOR &v){
    /** Increment operator. Allows: v1 = v1 + v2, where v1, v2 - are vectors. This operation changes v1 */

    x+=v.x;
    y+=v.y;
    z+=v.z;
  }
  void operator-=(const VECTOR &v){
    /** Decrement operator. Allows: v1 = v1 - v2, where v1, v2 - are vectors. This operation changes v1 */

    x-=v.x;
    y-=v.y;
    z-=v.z;
  }
  void operator+=(const double& v){
    /** Increment operator for scalar. Allows: v1 = v1 + A, where v1 - vector, A - scalar. This operation changes v1 */

    x+=v;
    y+=v;
    z+=v;
  }
  void operator-=(const double& v){
    /** Decrement operator for scalar. Allows: v1 = v1 - A, where v1 - vector, A - scalar. This operation changes v1 */

    x-=v;
    y-=v;
    z-=v;
  }

  void operator*=(const double &v){
    /** Scaling operator. Allows: v1 = v1 * A, where v1 - vector, A - scalar. This operation changes v1 */

    x*=v;
    y*=v;
    z*=v;
  }
  void operator/=(const double &v){
    /** Down-scaling operator. Allows: v1 = v1 / A, where v1 - vector, A - scalar. This operation changes v1 */

    x/=v;
    y/=v;
    z/=v;
  }


  //------------------ Friend functions -----------------------------     


  friend int operator == (const VECTOR& v1, const VECTOR& v2){
    return ((v1.x==v2.x)&&(v1.y==v2.y)&&(v1.z==v2.z));
  }
  friend int operator != (const VECTOR& v1, const VECTOR& v2){
    return ((v1.x!=v2.x)||(v1.y!=v2.y)||(v1.z!=v2.z));
  }
  friend double operator*(const VECTOR& v1, const VECTOR& v2){
  /**
    Returns the scalar (dot) product of two vectors
  */

    return (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z);
  }
  friend VECTOR operator*(const double& f,  const VECTOR& v1){
  /**
    Returns the vector multiplied from the left by a scalar. The original vector stays unchanged
  */

    VECTOR tmp;
    tmp.x=f*v1.x;
    tmp.y=f*v1.y;
    tmp.z=f*v1.z;
    return tmp;
  }
  friend VECTOR operator*(const VECTOR &v1, const double  &f){
  /**
    Returns the vector multiplied from the right by a scalar. The original vector stays unchanged
  */

    VECTOR tmp;
    tmp.x=v1.x*f;
    tmp.y=v1.y*f;
    tmp.z=v1.z*f;
    return tmp;
  }
  friend VECTOR operator/(const VECTOR &v1, const double &f){
  /**
    Returns the vector divided by a scalar. The original vector (numerator) stays unchanged
  */

    VECTOR tmp;
    tmp.x=v1.x/f;
    tmp.y=v1.y/f;
    tmp.z=v1.z/f;
    return tmp;
  }

  friend ostream& operator<<(ostream &strm,VECTOR ob){
  /** Formatted output of the vector to specified stream
  */

    int VECTOR_PRECISION = 8;
    int VECTOR_WIDTH = 15;
    strm.setf(ios::showpoint);
    strm.precision(VECTOR_PRECISION);
    strm.width(VECTOR_WIDTH);
    strm<<right;
    strm<<ob.x<<"  ";

    strm.setf(ios::showpoint);
    strm.precision(VECTOR_PRECISION);
    strm.width(VECTOR_WIDTH);
    strm<<right;
    strm<<ob.y<<"  ";

    strm.setf(ios::showpoint);
    strm.precision(VECTOR_PRECISION);
    strm.width(VECTOR_WIDTH);
    strm<<right;
    strm<<ob.z<<"  ";

    return strm;
  }
  friend istream& operator>>(istream& strm,VECTOR& ob){
  /** Formatted input of the vector from specified stream
  */

    strm>>ob.x>>ob.y>>ob.z;
    return strm;
  }


  //------------- Friend functions with other classes ---------------------------------


  friend VECTOR operator*(const MATRIX& m,  const VECTOR& v);    ///< Multiplication of vector and matrix: m * v
  friend VECTOR operator*(const MATRIX3x3& m,  const VECTOR& v); ///< Multiplication of vector and matrix: m * v




};


VECTOR cross(const double k, const VECTOR &v1, const VECTOR &v2);



typedef std::vector<VECTOR> VECTORList;  ///< Data type that holds a vector of VECTOR objects
typedef std::vector<vector<VECTOR> > VECTORMap;  ///< Data type that holds the table (grid) of VECTOR objects


//-------- IO functions --------
void set_value(int& is_defined, VECTOR& value, boost::python::object obj, std::string attrName);
void save(boost::property_tree::ptree& pt,std::string path,VECTOR& vt);
void save(boost::property_tree::ptree& pt,std::string path, char path_separator, VECTOR& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<VECTOR>& vt);
void save(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<VECTOR>& vt);

void load(boost::property_tree::ptree& pt,std::string path,VECTOR& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path, char path_separator, VECTOR& vt, int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<VECTOR>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<VECTOR>& vt,int& status);

//void export_VECTOR_objects();


}// namespace liblinalg
}// liblibra

#endif // VECTOR_H


