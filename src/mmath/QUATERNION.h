#ifndef QUATERNION_H
#define QUATERNION_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vector> 
  
using namespace std; 

//============================================================
// Forward declared dependencies
class VECTOR;
class QUATERNION;
class MATRIX3x3;
class MATRIX;
    
// ====================== Quaternion Class ====================

class QUATERNION{

  public:
  double Lt,Lx,Ly,Lz;
  
  //---------------------- Constructor ---------------------------------
  QUATERNION() { Lt = Lx = Ly = Lz =0; } 
  QUATERNION(double t,double x,double y,double z) { Lt = t; Lx = x; Ly =y; Lz = z; } 
  // Copy constructor
  QUATERNION(const QUATERNION& q){ Lt = q.Lt; Lx = q.Lx; Ly = q.Ly; Lz = q.Lz; }

  //--------------------- Destructor -----------------------------------
  ~QUATERNION(){ }

  //------------------- Initializations ---------------------------------
  void init(double t,double x,double y, double z) { Lt = t; Lx = x; Ly =y; Lz = z;  } 
  void init(QUATERNION& q){ Lt = q.Lt; Lx = q.Lx; Ly = q.Ly; Lz = q.Lz; } 
  void init(const QUATERNION& q){ Lt = q.Lt; Lx = q.Lx; Ly = q.Ly; Lz = q.Lz; } 

  //------------------ Basic quaternion operations ---------------------
  QUATERNION inverse(){
    QUATERNION tmp;
    double norm=(Lt*Lt+Lx*Lx+Ly*Ly+Lz*Lz);
    if(norm>0.0){
      tmp.Lt =  Lt/norm;
      tmp.Lx = -Lx/norm;
      tmp.Ly = -Ly/norm;
      tmp.Lz = -Lz/norm;
    }else{cout<<"Zero quaternion norm! Can not divide!"<<endl; }
    return tmp;
  }
  QUATERNION conj(){
    QUATERNION res;
    res.Lt =  Lt;
    res.Lx = -Lx;
    res.Ly = -Ly;
    res.Lz = -Lz;
    return res;
  }
  double sqal()  {return Lt;}
  VECTOR vect();

  //------------------ Small functions on quaternions -----------------
  double norm()       {return (Lt*Lt+Lx*Lx+Ly*Ly+Lz*Lz);}
  double mod()        {return sqrt((Lt*Lt+Lx*Lx+Ly*Ly+Lz*Lz));}
  void normalize()    {double norm = sqrt(Lt*Lt+Lx*Lx+Ly*Ly+Lz*Lz); Lt = Lt/norm; Lx = Lx/norm; Ly = Ly/norm; Lz = Lz/norm; }

  //---------------------- Operators -----------------------------------
  QUATERNION operator ~(){
    Lx = -Lx;
    Ly = -Ly;
    Lz = -Lz;
    return *this;
  }
  QUATERNION operator-(){
    QUATERNION tmp;
    tmp.Lt = -Lt;
    tmp.Lx = -Lx;
    tmp.Ly = -Ly;
    tmp.Lz = -Lz;
    return tmp;
  }
  QUATERNION operator+(QUATERNION& q){
    QUATERNION tmp;
    tmp.Lt = Lt + q.Lt;
    tmp.Lx = Lx + q.Lx;
    tmp.Ly = Ly + q.Ly;
    tmp.Lz = Lz + q.Lz;
    return tmp;
  }
  QUATERNION operator+(const QUATERNION& q){
    QUATERNION tmp;
    tmp.Lt = Lt + q.Lt;
    tmp.Lx = Lx + q.Lx;
    tmp.Ly = Ly + q.Ly;
    tmp.Lz = Lz + q.Lz;
    return tmp;
  }
  QUATERNION operator-(QUATERNION& q){
    QUATERNION tmp;
    tmp.Lt = Lt - q.Lt;
    tmp.Lx = Lx - q.Lx;
    tmp.Ly = Ly - q.Ly;
    tmp.Lz = Lz - q.Lz;
    return tmp;
  }
  QUATERNION operator-(const QUATERNION& q){
    QUATERNION tmp;
    tmp.Lt = Lt - q.Lt;
    tmp.Lx = Lx - q.Lx;
    tmp.Ly = Ly - q.Ly;
    tmp.Lz = Lz - q.Lz;
    return tmp;
  }
  QUATERNION operator*(QUATERNION& ob){
    QUATERNION tmp;
    tmp.Lt = Lt*ob.Lt - (Lx*ob.Lx+Ly*ob.Ly+Lz*ob.Lz);
    tmp.Lx = Lt*ob.Lx + ob.Lt*Lx + (Ly*ob.Lz-ob.Ly*Lz);
    tmp.Ly = Lt*ob.Ly + ob.Lt*Ly + (Lz*ob.Lx-ob.Lz*Lx);
    tmp.Lz = Lt*ob.Lz + ob.Lt*Lz + (Lx*ob.Ly-ob.Lx*Ly);
    return tmp;
  }
  QUATERNION operator*(const QUATERNION& ob){
    QUATERNION tmp;
    tmp.Lt = Lt*ob.Lt - (Lx*ob.Lx+Ly*ob.Ly+Lz*ob.Lz);
    tmp.Lx = Lt*ob.Lx + ob.Lt*Lx + (Ly*ob.Lz-ob.Ly*Lz);
    tmp.Ly = Lt*ob.Ly + ob.Lt*Ly + (Lz*ob.Lx-ob.Lz*Lx);
    tmp.Lz = Lt*ob.Lz + ob.Lt*Lz + (Lx*ob.Ly-ob.Lx*Ly);
    return tmp;
  }
  QUATERNION operator*(double f){
    QUATERNION tmp;
    tmp.Lt = f*Lt;
    tmp.Lx = f*Lx;
    tmp.Ly = f*Ly;
    tmp.Lz = f*Lz;
    return tmp;
  }
  QUATERNION& operator=(const QUATERNION &q){
    Lt = q.Lt;
    Lx = q.Lx;
    Ly = q.Ly;
    Lz = q.Lz;
    return *this;
  }
  QUATERNION& operator=(const double &q){
    Lt = q;
    Lx = q;
    Ly = q;
    Lz = q;
    return *this;
  }
  void operator+=(const QUATERNION &q){
    Lt += q.Lt;
    Lx += q.Lx;
    Ly += q.Ly;
    Lz += q.Lz;
  }
  void operator-=(const QUATERNION &q){
    Lt -= q.Lt;
    Lx -= q.Lx;
    Ly -= q.Ly;
    Lz -= q.Lz;
  }
  void operator*=(const double &q){
    Lt *= q;
    Lx *= q;
    Ly *= q;
    Lz *= q;
  }
  void operator/=(const double &q){
    Lt /= q;
    Lx /= q;
    Ly /= q;
    Lz /= q;
  }

  //------------------ Friend functions -----------------------------     
  friend QUATERNION operator*(QUATERNION& q,  const VECTOR& v);    // Multiplication of vector and quaternion;
  friend QUATERNION operator*(const MATRIX& m,QUATERNION& q);      // Multiplication of Matrix(4x4) and quaternion
                                                                   // defined similar to multiplication of 4D vector
  friend QUATERNION operator*(MATRIX& m,QUATERNION& q);            // Multiplication of Matrix(4x4) and quaternion
                                                                   // defined similar to multiplication of 4D vector
  friend QUATERNION operator*(MATRIX& m,const QUATERNION& q);      // Multiplication of Matrix(4x4) and quaternion
                                                                   // defined similar to multiplication of 4D vector
  friend QUATERNION operator*(const MATRIX& m,const QUATERNION& q);// Multiplication of Matrix(4x4) and quaternion
                                                                   // defined similar to multiplication of 4D vector

  friend QUATERNION operator*(const VECTOR& v, QUATERNION& q);    // Multiplication of vector and quaternion;
  friend QUATERNION operator*(double f,QUATERNION& q);
  friend double dot_prod(QUATERNION&, QUATERNION&);          // Normal dot product of 2 quaternions as 4D vectors
  friend double dot_prod(const QUATERNION&, const QUATERNION&);

  friend ostream& operator<<(ostream &strm,QUATERNION ob){
    int QUATERNION_PRECISION = 8;
    int QUATERNION_WIDTH = 15;
    strm.setf(ios::showpoint);
    strm.precision(QUATERNION_PRECISION);
    strm.width(QUATERNION_WIDTH);
    strm<<right;
    strm<<ob.Lt<<"  ";

    strm.setf(ios::showpoint);
    strm.precision(QUATERNION_PRECISION);
    strm.width(QUATERNION_WIDTH);
    strm<<right;
    strm<<ob.Lx<<"  ";

    strm.setf(ios::showpoint);
    strm.precision(QUATERNION_PRECISION);
    strm.width(QUATERNION_WIDTH);
    strm<<right;
    strm<<ob.Ly<<"  ";

    strm.setf(ios::showpoint);
    strm.precision(QUATERNION_PRECISION);
    strm.width(QUATERNION_WIDTH);
    strm<<right;
    strm<<ob.Lz<<"  ";

    return strm;
  }
  friend istream& operator>>(istream& strm,QUATERNION& ob){
    strm>>ob.Lt>>ob.Lx>>ob.Ly>>ob.Lz;
    return strm;
  }


};

#endif  // QUATERNION_H




