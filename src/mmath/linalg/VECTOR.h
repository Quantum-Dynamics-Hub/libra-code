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

#ifndef VECTOR_H
#define VECTOR_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vector>
using namespace std;

#include "../../io/libio.h"
using namespace libio;


namespace libmmath{
namespace liblinalg{

//============================================================
// Forward declared dependencies
class VECTOR;
class MATRIX3x3;
class MATRIX;
//==============================================================

class VECTOR
{
  public:
  double x,y,z;
  int is_transposed; // 1 = true = horizontal notation;
                     // 0 = false= vertical notation(default);         

  //---------------------- Constructor ---------------------------------


  VECTOR(){  is_transposed=0; x = 0.0; y = 0.0; z = 0.0; }
  VECTOR(double xval,double yval,double zval){  is_transposed=0; x=xval; y=yval; z=zval; }

  // Copy constructor
  VECTOR(const VECTOR& obj){ is_transposed = obj.is_transposed; x=obj.x; y=obj.y; z=obj.z; }

  //--------------------- Destructor -----------------------------------


  ~VECTOR(){}

  //------------------- Initializations ---------------------------------


  void init(double val){x=y=z=val;}
  void init(double xval,double yval,double zval){x=xval;y=yval;z=zval; }
  void init(VECTOR& v){ is_transposed=v.is_transposed; x=v.x; y=v.y; z=v.z; }
  void init(const VECTOR& v){ is_transposed=v.is_transposed; x=v.x; y=v.y; z=v.z; }

  //------------------- Special vectors --------------------------------


  VECTOR unit(){
    double res; VECTOR U;
    res=sqrt(x*x+y*y+z*z);
    U.x=x/res;
     U.y=y/res;
    U.z=z/res;
    return U;
  } 

  //------------------ Small functions on vectors ----------------------


  inline double length(void){  return sqrt(x*x+y*y+z*z); }
  inline double length2(void){ return (x*x+y*y+z*z); }
  void normalize(void){
    double res;
    res=sqrt(x*x+y*y+z*z);
    if(res!=0.0){ x/=res; y/=res; z/=res; }
    else{ x = 0.0; y = 0.0; z = 0.0; }
  }
  inline void cross(VECTOR &v1, VECTOR &v2) {
    x= v1.y*v2.z-v2.y*v1.z;
    y= v2.x*v1.z-v1.x*v2.z;
    z= v1.x*v2.y-v2.x*v1.y;
  }



  //---------------------- Operators -----------------------------------


  VECTOR operator-(){
    VECTOR tmp;
      tmp.x = -x;
      tmp.y = -y;
      tmp.z = -z;
      return tmp;
  }            
  VECTOR operator+(VECTOR& v){
    VECTOR tmp;
    tmp.x=x+v.x;
    tmp.y=y+v.y;
    tmp.z=z+v.z;
    return tmp;
  }
  VECTOR operator+(const VECTOR& v){
    VECTOR tmp;
    tmp.x=x+v.x;
    tmp.y=y+v.y;
    tmp.z=z+v.z;
    return tmp;
  }
  VECTOR operator-(VECTOR& v){
    VECTOR tmp;
    tmp.x=x-v.x;
    tmp.y=y-v.y;
    tmp.z=z-v.z;
    return tmp;
  }
  VECTOR operator-(const VECTOR& v){
    VECTOR tmp;
    tmp.x=x-v.x;
    tmp.y=y-v.y;
    tmp.z=z-v.z;
    return tmp;
  }
  VECTOR operator+(const double& v){    // Addition of vector and scalar
    VECTOR tmp;
    tmp.x=x+v;
    tmp.y=y+v;
    tmp.z=z+v;
    return tmp;
  }
  VECTOR operator-(const double& v){    // Subtraction of vector and scalar;
    VECTOR tmp;
    tmp.x=x-v;
    tmp.y=y-v;
    tmp.z=z-v;
    return tmp;
  }
  VECTOR operator=(const VECTOR &v){  // Important to return object, not the reference!!! otherwise crap happens on the Python side
    x=v.x;
    y=v.y;
    z=v.z;
    return *this;
  }
  VECTOR operator=(const double &v){  // same business with the reference
    x=v;
    y=v;
    z=v;
    return *this;
  }
  void operator+=(const VECTOR &v){
    x+=v.x;
    y+=v.y;
    z+=v.z;
  }
  void operator-=(const VECTOR &v){
    x-=v.x;
    y-=v.y;
    z-=v.z;
  }
  void operator+=(const double& v){
    x+=v;
    y+=v;
    z+=v;
  }
  void operator-=(const double& v){
    x-=v;
    y-=v;
    z-=v;
  }

  void operator*=(const double &v){
    x*=v;
    y*=v;
    z*=v;
  }
  void operator/=(const double &v){
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
    return (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z);
  }
  friend VECTOR operator*(const double& f,  const VECTOR& v1){
    VECTOR tmp;
    tmp.x=f*v1.x;
    tmp.y=f*v1.y;
    tmp.z=f*v1.z;
    return tmp;
  }
  friend VECTOR operator*(const VECTOR &v1, const double  &f){
    VECTOR tmp;
    tmp.x=v1.x*f;
    tmp.y=v1.y*f;
    tmp.z=v1.z*f;
    return tmp;
  }
  friend VECTOR operator/(const VECTOR &v1, const double &f){
    VECTOR tmp;
    tmp.x=v1.x/f;
    tmp.y=v1.y/f;
    tmp.z=v1.z/f;
    return tmp;
  }
  friend VECTOR cross(const double k, const VECTOR &v1, const VECTOR &v2) {
    VECTOR tmp;
    tmp.x=k*(v1.y*v2.z-v2.y*v1.z);
    tmp.y=k*(v1.z*v2.x-v2.z*v1.x);
    tmp.z=k*(v1.x*v2.y-v2.x*v1.y);
    return tmp;
  }
  friend ostream& operator<<(ostream &strm,VECTOR ob){
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
    strm>>ob.x>>ob.y>>ob.z;
    return strm;
  }


  //------------- Friend functions with other classes ---------------------------------


  friend VECTOR operator*(const MATRIX& m,  const VECTOR& v);    // Multiplication of vector and matrix;
  friend VECTOR operator*(const MATRIX3x3& m,  const VECTOR& v);    // Multiplication of vector and matrix;



};

typedef std::vector<VECTOR> VECTORList;
typedef std::vector<vector<VECTOR> > VECTORMap;


//-------- IO functions --------
void set_value(int& is_defined, VECTOR& value, boost::python::object obj, std::string attrName);
void save(boost::property_tree::ptree& pt,std::string path,VECTOR& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<VECTOR>& vt);
void load(boost::property_tree::ptree& pt,std::string path,VECTOR& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<VECTOR>& vt,int& status);


void export_VECTOR_objects();

}// namespace liblinalg
}// libmmath

#endif // VECTOR_H


