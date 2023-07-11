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

#ifndef MATRIX3x3_H
#define MATRIX3x3_H

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

//=============================================================================================

class MATRIX3x3{

  public:
  double xx, xy, xz,
         yx, yy, yz,
         zx, zy, zz;

  //---------------------- Constructor ---------------------------------


  MATRIX3x3(){ xx = xy = xz = yx = yy = yz = zx = zy = zz = 0.0;  }
  MATRIX3x3(const VECTOR& r1, const VECTOR& r2, const VECTOR& r3);

  // Copy constructor
  MATRIX3x3(const MATRIX3x3& ob){
  xx = ob.xx; xy = ob.xy; xz = ob.xz;
  yx = ob.yx; yy = ob.yy; yz = ob.yz;
  zx = ob.zx; zy = ob.zy; zz = ob.zz;
  }
  
  //--------------------- Destructor -----------------------------------


  ~MATRIX3x3(){ }
  
  //------------------- Initializations ---------------------------------


  void init(double x){ xx = xy = xz = yx = yy = yz = zx = zy = zz = x; }
  void init(VECTOR& r1, VECTOR& r2, VECTOR& r3);
  void init(const VECTOR& r1, const VECTOR& r2, const VECTOR& r3);
  void init(MATRIX3x3& m){
    xx = m.xx;  xy = m.xy;  xz = m.xz;
    yx = m.yx;  yy = m.yy;  yz = m.yz;
    zx = m.zx;  zy = m.zy;  zz = m.zz;
  }
  void init(const MATRIX3x3& m){
    xx = m.xx;  xy = m.xy;  xz = m.xz;
    yx = m.yx;  yy = m.yy;  yz = m.yz;
    zx = m.zx;  zy = m.zy;  zz = m.zz;
  }
  
  //------------------- Special matrices ---------------------------------


  void identity(){ xy = xz = yx = yz = zx = zy = 0.0; xx = yy = zz = 1.0; }
  void diag(double x){ xy = xz = yx = yz = zx = zy = 0.0; xx = yy = zz = x; }
  void diag(double x, double y, double z){ xy = xz = yx = yz = zx = zy = 0.0; xx = x; yy = y; zz = z; }
  void skew(VECTOR v);                      // Skew-symmetric matrix
  void Rx(double phi){
    double cs = cos(phi);
    double sn = sin(phi);
    xx = 1.0; xy = 0.0; xz = 0.0;
    yx = 0.0; yy = cs;  yz = -sn;
    zx = 0.0; zy = sn;  zz = cs;
  }
  void Ry(double phi){
    double cs = cos(phi);
    double sn = sin(phi);
    xx = cs;  xy = 0.0; xz = sn;
    yx = 0.0; yy = 1.0; yz = 0.0;
    zx = -sn; zy = 0.0; zz = cs;
  }
  void Rz(double phi){
    double cs = cos(phi);
    double sn = sin(phi);
    xx = cs;  xy = -sn; xz = 0.0;
    yx = sn;  yy = cs;  yz = 0.0;
    zx = 0.0; zy = 0.0; zz = 1.0;
  }
  void Rotation(const VECTOR& r);           // Create a rotation matrix for
                                            // rotation around direction given by direction
                                            // of the argument vector, and on the
                                            // amount given by the norm of the vector


  //------------------ Basic matrix operations --------------------------


  MATRIX3x3 inverse();
  void transpose(){
    double tmp;
    tmp = xy; xy = yx; yx = tmp;
    tmp = xz; xz = zx; zx = tmp;
    tmp = yz; yz = zy; yz = tmp;
  }
  MATRIX3x3 T(){
    MATRIX3x3 m;
    m.xx = xx; m.xy = yx; m.xz = zx;
    m.yx = xy; m.yy = yy; m.yz = zy;
    m.zx = xz; m.zy = yz; m.zz = zz;
    return m;
  }
  void eigen(MATRIX3x3&, MATRIX3x3&);
  void eigen(MATRIX3x3&, MATRIX3x3&,double);



  //------------------ Small functions on matrices ----------------------


  double Determinant(){ return (xx*(yy*zz - yz*zy) - yx*(xy*zz-zy*xz) + zx*(xy*yz - yy*xz));  }
  double tr(){ return (xx+yy+zz); }
  void get_vectors(VECTOR& r1,VECTOR& r2,VECTOR& r3);
  void tensor_product(VECTOR v1,VECTOR v2);



  //---------------------- Operators -----------------------------------


  MATRIX3x3 operator-() const{
    MATRIX3x3 m;
    m.xx = -xx;  m.xy = -xy;  m.xz = -xz;
    m.yx = -yx;  m.yy = -yy;  m.yz = -yz;
    m.zx = -zx;  m.zy = -zy;  m.zz = -zz;
    return m;
  }
  MATRIX3x3 operator*(MATRIX3x3 ob) const{
    MATRIX3x3 m;
    m.xx = xx*ob.xx + xy*ob.yx + xz*ob.zx; m.xy = xx*ob.xy + xy*ob.yy + xz*ob.zy; m.xz = xx*ob.xz + xy*ob.yz + xz*ob.zz;
    m.yx = yx*ob.xx + yy*ob.yx + yz*ob.zx; m.yy = yx*ob.xy + yy*ob.yy + yz*ob.zy; m.yz = yx*ob.xz + yy*ob.yz + yz*ob.zz;
    m.zx = zx*ob.xx + zy*ob.yx + zz*ob.zx; m.zy = zx*ob.xy + zy*ob.yy + zz*ob.zy; m.zz = zx*ob.xz + zy*ob.yz + zz*ob.zz;
    return m;
  }
  MATRIX3x3 operator+(MATRIX3x3 ob) const{
    MATRIX3x3 m;
    m.xx = xx + ob.xx; m.xy = xy + ob.xy; m.xz = xz + ob.xz;
    m.yx = yx + ob.yx; m.yy = yy + ob.yy; m.yz = yz + ob.yz;
    m.zx = zx + ob.zx; m.zy = zy + ob.zy; m.zz = zz + ob.zz;
    return m;
  }
  MATRIX3x3 operator-(MATRIX3x3 ob) const{
    MATRIX3x3 m;
    m.xx = xx - ob.xx; m.xy = xy - ob.xy; m.xz = xz - ob.xz;
    m.yx = yx - ob.yx; m.yy = yy - ob.yy; m.yz = yz - ob.yz;
    m.zx = zx - ob.zx; m.zy = zy - ob.zy; m.zz = zz - ob.zz;
    return m;
  }
  MATRIX3x3 operator+=(MATRIX3x3 ob){
    xx = xx + ob.xx; xy = xy + ob.xy; xz = xz + ob.xz;
    yx = yx + ob.yx; yy = yy + ob.yy; yz = yz + ob.yz;
    zx = zx + ob.zx; zy = zy + ob.zy; zz = zz + ob.zz;
    return *this;
  }
  MATRIX3x3 operator-=(MATRIX3x3 ob){
    xx = xx - ob.xx; xy = xy - ob.xy; xz = xz - ob.xz;
    yx = yx - ob.yx; yy = yy - ob.yy; yz = yz - ob.yz;
    zx = zx - ob.zx; zy = zy - ob.zy; zz = zz - ob.zz;
    return *this;
  }

  MATRIX3x3 operator*=(double f){
    xx *= f; xy *= f; xz *= f;
    yx *= f; yy *= f; yz *= f;
    zx *= f; zy *= f; zz *= f;
    return *this;
  }
  MATRIX3x3 operator/=(double f){
    xx /= f; xy /= f; xz /= f;
    yx /= f; yy /= f; yz /= f;
    zx /= f; zy /= f; zz /= f;
    return *this;
  }


  MATRIX3x3 operator/(double num) const{
    MATRIX3x3 m;
    m.xx = xx/num; m.xy = xy/num; m.xz = xz/num;
    m.yx = yx/num; m.yy = yy/num; m.yz = yz/num;
    m.zx = zx/num; m.zy = zy/num; m.zz = zz/num;
    return m;
  }
  MATRIX3x3 operator=(MATRIX3x3 ob){
    xx = ob.xx; xy = ob.xy; xz = ob.xz;
    yx = ob.yx; yy = ob.yy; yz = ob.yz;
    zx = ob.zx; zy = ob.zy; zz = ob.zz;
    return *this;
  }
  MATRIX3x3 operator=(double num){
    xx = xy = xz = yx = yy = yz = zx = zy = zz = num;
    return *this;
  }

  //------------------ Friend functions -----------------------------     


  friend bool operator == (const MATRIX3x3& m1, const MATRIX3x3& m2){
    // Are matrices equal
    return  ( (m1.xx == m2.xx) && (m1.xy == m2.xy) && (m1.xz == m2.xz) &&
              (m1.yx == m2.yx) && (m1.yy == m2.yy) && (m1.yz == m2.yz) &&
              (m1.zx == m2.zx) && (m1.zy == m2.zy) && (m1.zz == m2.zz)
            );  
  }
  friend bool operator != (const MATRIX3x3& m1, const MATRIX3x3& m2){
    // Are matrices not equal
    return  (!( (m1.xx == m2.xx) && (m1.xy == m2.xy) && (m1.xz == m2.xz) &&
                (m1.yx == m2.yx) && (m1.yy == m2.yy) && (m1.yz == m2.yz) &&
                (m1.zx == m2.zx) && (m1.zy == m2.zy) && (m1.zz == m2.zz)
              )
            ); 
  }
  friend MATRIX3x3 operator*(const double& f,  const MATRIX3x3& m1){
    // Multiplication of matrix and double
    MATRIX3x3 m;
    m.xx = m1.xx * f; m.xy = m1.xy * f; m.xz = m1.xz * f;
    m.yx = m1.yx * f; m.yy = m1.yy * f; m.yz = m1.yz * f;
    m.zx = m1.zx * f; m.zy = m1.zy * f; m.zz = m1.zz * f;
    return m;
  }
  friend MATRIX3x3 operator*(const MATRIX3x3 &m1, const double  &f){
    // Multiplication of matrix and double
    MATRIX3x3 m;
    m.xx = m1.xx * f; m.xy = m1.xy * f; m.xz = m1.xz * f;
    m.yx = m1.yx * f; m.yy = m1.yy * f; m.yz = m1.yz * f;
    m.zx = m1.zx * f; m.zy = m1.zy * f; m.zz = m1.zz * f;
    return m;
  }

  friend ostream &operator<<(ostream &strm,MATRIX3x3 ob){
    strm.setf(ios::showpoint);
    int MATRIX_PRECISION = 8;
    int MATRIX_WIDTH = 15;

    strm.precision(MATRIX_PRECISION);
    strm.width(MATRIX_WIDTH);
    strm<<left;
    strm<<ob.xx<<"  "<<ob.xy<<"  "<<ob.xz<<endl;
    strm<<ob.yx<<"  "<<ob.yy<<"  "<<ob.yz<<endl;
    strm<<ob.zx<<"  "<<ob.zy<<"  "<<ob.zz<<endl;
    return strm;
  }

  friend istream& operator>>(istream& strm,MATRIX3x3 &ob){
    return strm;
  }


  //------------- Friend functions with other classes ---------------------------------

  
  friend VECTOR operator*(const MATRIX3x3& m,  const VECTOR& v); // Multiplication of vector and matrix;


};

typedef std::vector<MATRIX3x3> MATRIX3x3List;
typedef std::vector<vector<MATRIX3x3> > MATRIX3x3Map;

//-------- IO functions --------
void set_value(int& defined, MATRIX3x3& value, boost::python::object obj, std::string attrName);
void save(boost::property_tree::ptree& pt,std::string path,MATRIX3x3& vt);
void save(boost::property_tree::ptree& pt,std::string path, char path_separator, MATRIX3x3& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<MATRIX3x3>& vt);
void save(boost::property_tree::ptree& pt,std::string path, char path_separation, vector<MATRIX3x3>& vt);

void load(boost::property_tree::ptree& pt,std::string path,MATRIX3x3& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path, char path_separator, MATRIX3x3& vt, int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<MATRIX3x3>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<MATRIX3x3>& vt,int& status);


}// namespace liblinalg
}// namespace liblibra
//===============================================================================

#endif  // MATRIX3x3

