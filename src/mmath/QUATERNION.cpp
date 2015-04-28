#include "QUATERNION.h"
#include "VECTOR.h"

VECTOR QUATERNION::vect(){VECTOR v; v.x = Lx; v.y = Ly; v.z = Lz; return v;}

double dot_prod(QUATERNION& q1, QUATERNION& q2){
  return (q1.Lt*q2.Lt + q1.Lx*q2.Lx + q1.Ly*q2.Ly + q1.Lz*q2.Lz);
}
double dot_prod(const QUATERNION& q1, const QUATERNION& q2){
  return (q1.Lt*q2.Lt + q1.Lx*q2.Lx + q1.Ly*q2.Ly + q1.Lz*q2.Lz);
}

QUATERNION operator*(double f,QUATERNION& q){
  QUATERNION tmp;
  tmp.Lt = f*q.Lt;
  tmp.Lx = f*q.Lx;
  tmp.Ly = f*q.Ly;
  tmp.Lz = f*q.Lz;
  return tmp;
}



