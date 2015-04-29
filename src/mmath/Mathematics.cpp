//#include "../Mathematics.h"
#include "VECTOR.h"
#include "QUATERNION.h"
#include "MATRIX3x3.h"
#include "MATRIX.h"

namespace libmmath{

// Here we define friend functions for all classes

VECTOR operator*(const MATRIX& m,  const VECTOR& v){
  // Multiplication of vector and matrix;
  VECTOR res;
  res.x = m.M[0]*v.x + m.M[1]*v.y + m.M[2]*v.z;
  res.y = m.M[3]*v.x + m.M[4]*v.y + m.M[5]*v.z;
  res.z = m.M[6]*v.x + m.M[7]*v.y + m.M[8]*v.z;
  return res;
}


VECTOR operator*(const MATRIX3x3& m,  const VECTOR& v){
  // Multiplication of vector and matrix;
  VECTOR res;
  res.x = m.xx*v.x + m.xy*v.y + m.xz*v.z;
  res.y = m.yx*v.x + m.yy*v.y + m.yz*v.z;
  res.z = m.zx*v.x + m.zy*v.y + m.zz*v.z;
  return res;
}


QUATERNION operator*(QUATERNION& q,  const VECTOR& v){
  // Multiplication of vector and quaternion;
  QUATERNION tmp;
  QUATERNION V;
  V.Lt = 0.0;
  V.Lx = v.x;
  V.Ly = v.y;
  V.Lz = v.z;
  tmp = q*V;
  return tmp;
}
QUATERNION operator*(const VECTOR& v,  QUATERNION& q){
  // Multiplication of vector and quaternion;
  QUATERNION tmp;
  QUATERNION V;
  V.Lt = 0.0;
  V.Lx = v.x;
  V.Ly = v.y;
  V.Lz = v.z;
  tmp = V*q;
  return tmp;
}
QUATERNION operator*(const MATRIX& m,QUATERNION& q){
  // Multiplication of Matrix(4x4) and quaternion
  // defined similar to multiplication of 4D vector
  QUATERNION tmp;
  if(m.num_of_cols==4 && m.num_of_rows==4){
    tmp.Lt = m.M[0]*q.Lt  + m.M[1]*q.Lx  + m.M[2]*q.Ly  + m.M[3]*q.Lz;
    tmp.Lx = m.M[4]*q.Lt  + m.M[5]*q.Lx  + m.M[6]*q.Ly  + m.M[7]*q.Lz;
    tmp.Ly = m.M[8]*q.Lt  + m.M[9]*q.Lx  + m.M[10]*q.Ly + m.M[11]*q.Lz;
    tmp.Lz = m.M[12]*q.Lt + m.M[13]*q.Lx + m.M[14]*q.Ly + m.M[15]*q.Lz;
  }else{
  cout<<"Can not multiply quaternion and matrix of dimention other than 4x4"<<endl;
  }
  return tmp;
}
QUATERNION operator*(MATRIX& m,QUATERNION& q){
  // Multiplication of Matrix(4x4) and quaternion
  // defined similar to multiplication of 4D vector
  QUATERNION tmp;
  if(m.num_of_cols==4 && m.num_of_rows==4){
    tmp.Lt = m.M[0]*q.Lt  + m.M[1]*q.Lx  + m.M[2]*q.Ly  + m.M[3]*q.Lz;
    tmp.Lx = m.M[4]*q.Lt  + m.M[5]*q.Lx  + m.M[6]*q.Ly  + m.M[7]*q.Lz;
    tmp.Ly = m.M[8]*q.Lt  + m.M[9]*q.Lx  + m.M[10]*q.Ly + m.M[11]*q.Lz;
    tmp.Lz = m.M[12]*q.Lt + m.M[13]*q.Lx + m.M[14]*q.Ly + m.M[15]*q.Lz;
  }else{
  cout<<"Can not multiply quaternion and matrix of dimention other than 4x4"<<endl;
  }
  return tmp;
}
QUATERNION operator*(MATRIX& m,const QUATERNION& q){
  // Multiplication of Matrix(4x4) and quaternion
  // defined similar to multiplication of 4D vector
  QUATERNION tmp;
  if(m.num_of_cols==4 && m.num_of_rows==4){
    tmp.Lt = m.M[0]*q.Lt  + m.M[1]*q.Lx  + m.M[2]*q.Ly  + m.M[3]*q.Lz;
    tmp.Lx = m.M[4]*q.Lt  + m.M[5]*q.Lx  + m.M[6]*q.Ly  + m.M[7]*q.Lz;
    tmp.Ly = m.M[8]*q.Lt  + m.M[9]*q.Lx  + m.M[10]*q.Ly + m.M[11]*q.Lz;
    tmp.Lz = m.M[12]*q.Lt + m.M[13]*q.Lx + m.M[14]*q.Ly + m.M[15]*q.Lz;
  }else{
  cout<<"Can not multiply quaternion and matrix of dimention other than 4x4"<<endl;
  }
  return tmp;
}
QUATERNION operator*(const MATRIX& m,const QUATERNION& q){
  // Multiplication of Matrix(4x4) and quaternion
  // defined similar to multiplication of 4D vector
  QUATERNION tmp;
  if(m.num_of_cols==4 && m.num_of_rows==4){
    tmp.Lt = m.M[0]*q.Lt  + m.M[1]*q.Lx  + m.M[2]*q.Ly  + m.M[3]*q.Lz;
    tmp.Lx = m.M[4]*q.Lt  + m.M[5]*q.Lx  + m.M[6]*q.Ly  + m.M[7]*q.Lz;
    tmp.Ly = m.M[8]*q.Lt  + m.M[9]*q.Lx  + m.M[10]*q.Ly + m.M[11]*q.Lz;
    tmp.Lz = m.M[12]*q.Lt + m.M[13]*q.Lx + m.M[14]*q.Ly + m.M[15]*q.Lz;
  }else{
  cout<<"Can not multiply quaternion and matrix of dimention other than 4x4"<<endl;
  }
  return tmp;
}


}// namespace libmmath

