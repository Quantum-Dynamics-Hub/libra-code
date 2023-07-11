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

#include "VECTOR.h"
#include "QUATERNION.h"
#include "MATRIX3x3.h"
#include "MATRIX.h"

/// liblibra namespace
namespace liblibra{

/// liblinalg namespace
namespace liblinalg{

// Here we define friend functions for all classes


MATRIX operator ^(const VECTOR& v1, const VECTOR& v2){
/** Tensor product of two vectors   */

  MATRIX res(3,3);

  res.M[0] = v1.x*v2.x;   res.M[1] = v1.x*v2.y;  res.M[2] = v1.x*v2.z;
  res.M[3] = v1.y*v2.x;   res.M[4] = v1.y*v2.y;  res.M[5] = v1.y*v2.z;
  res.M[6] = v1.z*v2.x;   res.M[7] = v1.z*v2.y;  res.M[8] = v1.z*v2.z;

  return res;
}



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


QUATERNION operator*(const QUATERNION& q,  const VECTOR& v){
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

QUATERNION operator*(const VECTOR& v,  const QUATERNION& q){
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

QUATERNION operator*(const MATRIX& m, const QUATERNION& q){
  // Multiplication of Matrix(4x4) and quaternion
  // defined similar to multiplication of 4D vector
  QUATERNION tmp;
  if(m.n_cols==4 && m.n_rows==4){
    tmp.Lt = m.M[0]*q.Lt  + m.M[1]*q.Lx  + m.M[2]*q.Ly  + m.M[3]*q.Lz;
    tmp.Lx = m.M[4]*q.Lt  + m.M[5]*q.Lx  + m.M[6]*q.Ly  + m.M[7]*q.Lz;
    tmp.Ly = m.M[8]*q.Lt  + m.M[9]*q.Lx  + m.M[10]*q.Ly + m.M[11]*q.Lz;
    tmp.Lz = m.M[12]*q.Lt + m.M[13]*q.Lx + m.M[14]*q.Ly + m.M[15]*q.Lz;
  }else{
  cout<<"Can not multiply quaternion and matrix of dimention other than 4x4"<<endl;
  }
  return tmp;
}




void MATRIX_TO_QUATERNION(MATRIX& M,QUATERNION& L){
/***********************************************************
    Convert orientation matrix to corresponding quaternion 
                     representation
***********************************************************/

//----------- Quaternion initialization ---------
  double a,b,c,d,e,f;
  //****************************************************
  // from www.gamedev.ru/users/wat/articles/quaternions
  // *****************************************************/
                      
  double a00,a01,a02,a10,a11,a12,a20,a21,a22;
  double m[3][3];
  double ind_q[4];

  m[0][0] = a00 = M.M[0];
  m[0][1] = a01 = M.M[1];
  m[0][2] = a02 = M.M[2];
  m[1][0] = a10 = M.M[3];
  m[1][1] = a11 = M.M[4];
  m[1][2] = a12 = M.M[5];
  m[2][0] = a20 = M.M[6];
  m[2][1] = a21 = M.M[7];
  m[2][2] = a22 = M.M[8];

  int ind_i,ind_j,ind_k;
  int nxt[3] = {1,2,0};
  double  tr,s;

  tr = a00 + a11 + a22;
  
  if(tr>0.0){
     s = sqrt(tr + 1.0);
     L.Lt = s/2.0; 
  
     s = 0.5/s;

     L.Lx = (a12 - a21)*s;
     L.Ly = (a20 - a02)*s;
     L.Lz = (a01 - a10)*s;
  }else{
     ind_i = 0;
     if(a11>a00)              ind_i = 1;
     if(a22>m[ind_i][ind_i])  ind_i = 2;

     ind_j  = nxt[ind_i];
     ind_k  = nxt[ind_j];

     s = sqrt(m[ind_i][ind_i] - (m[ind_j][ind_j]+m[ind_k][ind_k])+1.0 );

     ind_q[ind_i] = s * 0.5;

     if(s != 0.0)  s = 0.5/s;

     ind_q[3]      = (m[ind_j][ind_k] - m[ind_k][ind_j]) * s;
     ind_q[ind_j]  = (m[ind_i][ind_j] + m[ind_j][ind_i]) * s;
     ind_q[ind_k]  = (m[ind_i][ind_k] + m[ind_k][ind_i]) * s;

     L.Lx = ind_q[0];
     L.Ly = ind_q[1];
     L.Lz = ind_q[2];
     L.Lt = ind_q[3];

  }


}


void MATRIX_TO_QUATERNION(MATRIX3x3& M,QUATERNION& L){
/***********************************************************
    Convert orientation matrix to corresponding quaternion
                     representation
***********************************************************/

//----------- Quaternion initialization ---------
  double a,b,c,d,e,f;
  //****************************************************
  // from www.gamedev.ru/users/wat/articles/quaternions
  // *****************************************************/

  double a00,a01,a02,a10,a11,a12,a20,a21,a22;
  double m[3][3];
  double ind_q[4];

  m[0][0] = a00 = M.xx;
  m[0][1] = a01 = M.xy;
  m[0][2] = a02 = M.xz;
  m[1][0] = a10 = M.yx;
  m[1][1] = a11 = M.yy;
  m[1][2] = a12 = M.yz;
  m[2][0] = a20 = M.zx;
  m[2][1] = a21 = M.zy;
  m[2][2] = a22 = M.zz;

  int ind_i,ind_j,ind_k;
  int nxt[3] = {1,2,0};
  double  tr,s;

  tr = a00 + a11 + a22;

  if(tr>0.0){
     s = sqrt(tr + 1.0);
     L.Lt = s/2.0;

     s = 0.5/s;

     L.Lx = (a12 - a21)*s;
     L.Ly = (a20 - a02)*s;
     L.Lz = (a01 - a10)*s;
  }else{
     ind_i = 0;
     if(a11>a00)              ind_i = 1;
     if(a22>m[ind_i][ind_i])  ind_i = 2;

     ind_j  = nxt[ind_i];
     ind_k  = nxt[ind_j];

     s = sqrt(m[ind_i][ind_i] - (m[ind_j][ind_j]+m[ind_k][ind_k])+1.0 );

     ind_q[ind_i] = s * 0.5;

     if(s != 0.0)  s = 0.5/s;

     ind_q[3]      = (m[ind_j][ind_k] - m[ind_k][ind_j]) * s;
     ind_q[ind_j]  = (m[ind_i][ind_j] + m[ind_j][ind_i]) * s;
     ind_q[ind_k]  = (m[ind_i][ind_k] + m[ind_k][ind_i]) * s;

     L.Lx = ind_q[0];
     L.Ly = ind_q[1];
     L.Lz = ind_q[2];
     L.Lt = ind_q[3];
  }
}


void QUATERNION_TO_MATRIX(QUATERNION& L,MATRIX& M){
/**************************************************************
   Convert quaternion to corresponding orientation matrix 
                        representation
***************************************************************/
   double q0,q1,q2,q3;

   q0 = L.Lt;
   q1 = L.Lx;
   q2 = L.Ly;
   q3 = L.Lz;

   M.M[0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
   M.M[1] = 2.0*(q1*q2 + q0*q3);
   M.M[2] = 2.0*(q1*q3 - q0*q2);
   M.M[3] = 2.0*(q1*q2 - q0*q3);
   M.M[4] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
   M.M[5] = 2.0*(q2*q3 + q0*q1);
   M.M[6] = 2.0*(q1*q3 + q0*q2);
   M.M[7] = 2.0*(q2*q3 - q0*q1);
   M.M[8] = q0*q0 - q1*q1 - q2*q2 + q3*q3;


}

void QUATERNION_TO_MATRIX(QUATERNION& L,MATRIX3x3& M){
/**************************************************************
   Convert quaternion to corresponding orientation matrix
                        representation
***************************************************************/
   double q0,q1,q2,q3;

   q0 = L.Lt;
   q1 = L.Lx;
   q2 = L.Ly;
   q3 = L.Lz;

   M.xx = q0*q0 + q1*q1 - q2*q2 - q3*q3;
   M.xy = 2.0*(q1*q2 + q0*q3);
   M.xz = 2.0*(q1*q3 - q0*q2);
   M.yx = 2.0*(q1*q2 - q0*q3);
   M.yy = q0*q0 - q1*q1 + q2*q2 - q3*q3;
   M.yz = 2.0*(q2*q3 + q0*q1);
   M.zx = 2.0*(q1*q3 + q0*q2);
   M.zy = 2.0*(q2*q3 - q0*q1);
   M.zz = q0*q0 - q1*q1 - q2*q2 + q3*q3;


}



}// namespace liblinalg
}// namespace liblibra

