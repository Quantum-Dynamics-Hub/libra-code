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
  \file RigidBody_methods.cpp
  \brief The file implements RigidBody class interfaces - initializers and setters
  So far these methods are private, so they are not reflected in the auto-generated documentation    

*/

#include "RigidBody.h"

#include "../math_meigen/libmeigen.h"

/// liblibra namespace
namespace liblibra{

using namespace libmeigen;


/// librigidbody namespace
namespace librigidbody{

// ======================= Internal methods ============================

void RigidBody::calc_mass(int sz, double* m){
/**
  \brief compute mass

  Calculates mass of the RB 

  \param[in] sz The number of material points included on the rigid body
  \param[in] m Pointer to the array of material point masses
*/
  rb_mass = 0.0;              is_rb_mass = 1;
  for(int i=0;i<sz;i++){ rb_mass += m[i];  }
  rb_iM = (rb_mass==0.0)?0.0:(1.0/rb_mass); is_rb_iM = 1;
}
void RigidBody::calc_center_of_mass(int sz, double* m, VECTOR* r){
  rb_cm = 0.0;                is_rb_cm = 1;
  for(int i=0;i<sz;i++){ rb_cm += m[i]*r[i]; }
  rb_cm = rb_cm/rb_mass;
}

void RigidBody::calc_inertia_tensors(int sz, double* m, VECTOR* r){
  rb_I_I = 0.0;               is_rb_I_I = 1;
  VECTOR d;
  for(int i=0;i<sz;i++){
    d = r[i] - rb_cm;

    // Inertia tensor I_I
    rb_I_I.xx += m[i]*(d.y*d.y + d.z*d.z);
    rb_I_I.xy -= m[i]*d.x*d.y;
    rb_I_I.xz -= m[i]*d.x*d.z;

    rb_I_I.yy += m[i]*(d.x*d.x + d.z*d.z);
    rb_I_I.yz -= m[i]*d.y*d.z;

    rb_I_I.zz += m[i]*(d.x*d.x + d.y*d.y);

  }

  rb_I_I.yx = rb_I_I.xy;
  rb_I_I.zx = rb_I_I.xz;
  rb_I_I.zy = rb_I_I.yz;

}

void RigidBody::calc_orientations(int sz, double* m, VECTOR* r){
  //!!!math: rb_I_I = rb_A_I_to_e_T * rb_I_e * rb_A_I_to_e
  MATRIX H(3,3); MATRIX E(3,3); MATRIX C(3,3);
  H.set(0,0, rb_I_I.xx); H.set(0,1, rb_I_I.xy); H.set(0,2, rb_I_I.xz);
  H.set(1,0, rb_I_I.yx); H.set(1,1, rb_I_I.yy); H.set(1,2, rb_I_I.yz);
  H.set(2,0, rb_I_I.zx); H.set(2,1, rb_I_I.zy); H.set(2,2, rb_I_I.zz);  

  // H * C = C * E
  //solve_eigen_nosort(H, E, C, 0);
  solve_eigen(H, E, C, 0);

  //rb_I_I.eigen(rb_I_e, rb_A_I_to_e_T, MAX_NO); 

  rb_I_e.xx = E.get(0,0);
  rb_I_e.yy = E.get(1,1);
  rb_I_e.zz = E.get(2,2);
  rb_I_e.xy = rb_I_e.xz = rb_I_e.yx = rb_I_e.yz = rb_I_e.zx = rb_I_e.zy = 0.0;

  rb_A_I_to_e_T.xx = C.get(0,0);  rb_A_I_to_e_T.xy = C.get(0,1); rb_A_I_to_e_T.xz = C.get(0,2);
  rb_A_I_to_e_T.yx = C.get(1,0);  rb_A_I_to_e_T.yy = C.get(1,1); rb_A_I_to_e_T.yz = C.get(1,2);
  rb_A_I_to_e_T.zx = C.get(2,0);  rb_A_I_to_e_T.zy = C.get(2,1); rb_A_I_to_e_T.zz = C.get(2,2);

  is_rb_I_e = is_rb_A_I_to_e_T = 1;


  rb_A_I_to_e = rb_A_I_to_e_T.T(); //inverse(); //T(); 
  is_rb_A_I_to_e = 1;

  MATRIX_TO_QUATERNION(rb_A_I_to_e,rb_L);             is_rb_L = 1;

}

void RigidBody::calc_inverse_tensors(int sz){

  //--------------- invI_e -----------------------------
  if(sz>1){
    if(fabs(rb_I_e.Determinant())>minDet){
      rb_invI_e = rb_I_e.inverse();
    }
    else{ // This is for 2-center rigid bodies
      rb_invI_e.init(0.0);

      if(rb_I_e.xx>IEPS){  rb_invI_e.xx = 1.0/rb_I_e.xx;  }
      else{                  rb_invI_e.xx = 0.0;              }

      if(rb_I_e.yy>IEPS){  rb_invI_e.yy = 1.0/rb_I_e.yy;  }
      else{                  rb_invI_e.yy = 0.0;              }

      if(rb_I_e.zz>IEPS){  rb_invI_e.zz = 1.0/rb_I_e.zz;  }
      else{                  rb_invI_e.zz = 0.0;              }

    }// else (I_i[i].Determinant()==0)
  }
  else if(sz==1){  rb_invI_e.init(0.0);   }

  // Additional post-processing for linear fragments
  if(fabs(rb_invI_e.xx)>BIG) {
    cout<<"Warning: The fragment most likely was linear or otherwise degenerate\n";
    cout<<"Setting rb_invI_e.xx to zero!\n";
    rb_invI_e.xx = 0.0;
  }// if > BIG
  if(fabs(rb_invI_e.xy)>BIG) {
    cout<<"Warning: The fragment most likely was linear or otherwise degenerate\n";
    cout<<"Setting rb_invI_e.xy to zero!\n";
    rb_invI_e.xy = 0.0;
  }// if > BIG
  if(fabs(rb_invI_e.xz)>BIG) {
    cout<<"Warning: The fragment most likely was linear or otherwise degenerate\n";
    cout<<"Setting rb_invI_e.xz to zero!\n";
    rb_invI_e.xz = 0.0;
  }// if > BIG
  if(fabs(rb_invI_e.yx)>BIG) {
    cout<<"Warning: The fragment most likely was linear or otherwise degenerate\n";
    cout<<"Setting rb_invI_e.yx to zero!\n";
    rb_invI_e.yx = 0.0;
  }// if > BIG
  if(fabs(rb_invI_e.yy)>BIG) {
    cout<<"Warning: The fragment most likely was linear or otherwise degenerate\n";
    cout<<"Setting rb_invI_e.yy to zero!\n";
    rb_invI_e.yy = 0.0;
  }// if > BIG
  if(fabs(rb_invI_e.yz)>BIG) {
    cout<<"Warning: The fragment most likely was linear or otherwise degenerate\n";
    cout<<"Setting rb_invI_e.yz to zero!\n";
    rb_invI_e.yz = 0.0;
  }// if > BIG
  if(fabs(rb_invI_e.zx)>BIG) {
    cout<<"Warning: The fragment most likely was linear or otherwise degenerate\n";
    cout<<"Setting rb_invI_e.zx to zero!\n";
    rb_invI_e.zx = 0.0;
  }// if > BIG
  if(fabs(rb_invI_e.zy)>BIG) {
    cout<<"Warning: The fragment most likely was linear or otherwise degenerate\n";
    cout<<"Setting rb_invI_e.zy to zero!\n";
    rb_invI_e.zy = 0.0;
  }// if > BIG
  if(fabs(rb_invI_e.zz)>BIG) {
    cout<<"Warning: The fragment most likely was linear or otherwise degenerate\n";
    cout<<"Setting rb_invI_e.zz to zero!\n";
    rb_invI_e.zz = 0.0;
  }// if > BIG

  // Make it diagonal
  rb_invI_e.xy = rb_invI_e.xz = rb_invI_e.yx = rb_invI_e.yz = rb_invI_e.zx = rb_invI_e.zy = 0.0;
  is_rb_invI_e = 1;

  //-------------------- invI_I -----------------------------
  rb_invI_I = rb_A_I_to_e_T * rb_invI_e * rb_A_I_to_e;
  is_rb_invI_I = 1;

}

void RigidBody::calc_rot_constants(){
  rb_A = rb_invI_e.xx;  is_rb_A = 1;
  rb_B = rb_invI_e.yy;  is_rb_B = 1;
  rb_C = rb_invI_e.zz;  is_rb_C = 1;
}

void RigidBody::calc_rb_centers(int sz, VECTOR* r){
  if(rb_centers.size()>0){ rb_centers.clear(); rb_centers_size = 0; }
  VECTOR r_I; // coordinates of the center in body frame
  for(int i=0;i<sz;i++){
    r_I = rb_A_I_to_e*(r[i] - rb_cm);
    rb_centers.push_back(r_I);
    rb_centers_size++;
  }
}


}// namespace librigidbody

}// liblibra

