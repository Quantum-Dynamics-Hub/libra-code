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

#include "RigidBody.h"

namespace libdyn{
namespace librigidbody{

int RigidBody::init(int sz,double* m,VECTOR* r){
  calc_mass(sz, m);
  calc_center_of_mass(sz, m, r);
  calc_inertia_tensors(sz, m, r);
  calc_orientations(sz, m, r);
  calc_inverse_tensors(sz);
  calc_rot_constants();
  calc_rb_centers(sz,r);
  return 0;
}

int RigidBody::init(int sz,boost::python::list masses, boost::python::list positions){

  int sz1 = len(masses);
  int sz2 = len(positions);

  if(sz1!=sz2){ 
    std::cout<<"In RigidBody::init function\n";
    std::cout<<"The lenth of masses and positions lists are not equal\n";
    std::cout<<"exiting...\n"; exit(0);
  }
   
  double*  tmp; tmp = new double[sz1];
  VECTOR* vtmp; vtmp= new VECTOR[sz1];

  for(int i=0;i<sz1;i++){
    tmp[i]  = boost::python::extract<double>(masses[i]);   //    cout<<"mass "<<i<<"  = "<<tmp[i]<<endl;
    vtmp[i] = boost::python::extract<VECTOR>(positions[i]);//    cout<<"v "<<i<<"  = "<<vtmp[i]<<endl;
  }

  init(sz1,tmp,vtmp);

  delete [] tmp;
  delete [] vtmp;

}


int RigidBody::set_mass(const double& m){
  rb_mass = m;                  is_rb_mass = 1;
  rb_iM = (m==0.0)?0.0:(1.0/m); is_rb_iM = 1;
  return 0;
}

int RigidBody::set_position(const VECTOR& cm){
  rb_cm.init(cm);   is_rb_cm = 1;
  return 0;
}

int RigidBody::set_momentum(const VECTOR& p){
  rb_p.init(p);     is_rb_p = 1;
  // Also set velocity
  if(is_rb_mass){ if(rb_mass>0.0){   rb_v.init(p/rb_mass); is_rb_v = 1;  } }

  return 0;
}

int RigidBody::set_velocity(const VECTOR& v){
  rb_v.init(v);     is_rb_v = 1;
  // Also set linear momentum
  if(is_rb_mass){ rb_p.init(v*rb_mass); is_rb_p = 1;  }
  return 0;
}

int RigidBody::set_force(const VECTOR& f){
  rb_force.init(f);  is_rb_force = 1;
  return 0;
}

int RigidBody::set_torque(const VECTOR& t){
  rb_torque_e.init(t);  is_rb_torque_e = 1;
  return 0;
}

int RigidBody::set_forces_and_torques(int sz,VECTOR* r,VECTOR* f){
/* r - are the coordinates (in body frame) of the centers
       to which the forces are applied
   f - external forces (in lab frame) acting on the rigid body
*/
  rb_force = 0.0;      is_rb_force = 1;
  rb_torque_e = 0.0;   is_rb_torque_e = 1;

  for(int i=0;i<sz;i++){
    rb_force += f[i];
    rb_torque_e += cross(1.0,r[i], rb_A_I_to_e*f[i]);
  }
  return 0;
}


int RigidBody::set_inertia(const MATRIX3x3& Ie){
  rb_I_e.init(Ie);          is_rb_I_e = 1;
  // Need to put update of dependent variables here
  return 0;
}

int RigidBody::set_orientation(const MATRIX3x3& at){
  rb_A_I_to_e = at;                      is_rb_A_I_to_e = 1;
  rb_A_I_to_e_T = rb_A_I_to_e.T();       is_rb_A_I_to_e_T = 1;
  MATRIX_TO_QUATERNION(rb_A_I_to_e,rb_L);is_rb_L = 1;
  return 0;
}

int RigidBody::set_orientation(const VECTOR& u1, const VECTOR& u2, const VECTOR& u3){
  rb_A_I_to_e.init(u1,u2,u3);            is_rb_A_I_to_e = 1;
  rb_A_I_to_e_T = rb_A_I_to_e.T();       is_rb_A_I_to_e_T = 1;
  MATRIX_TO_QUATERNION(rb_A_I_to_e,rb_L);is_rb_L = 1;
  return 0;
}

int RigidBody::set_orientation(const QUATERNION& q){
  rb_L.init(q);                          is_rb_L = 1;
  QUATERNION_TO_MATRIX(rb_L,rb_A_I_to_e);is_rb_A_I_to_e = 1;
  rb_A_I_to_e_T = rb_A_I_to_e.T();       is_rb_A_I_to_e_T = 1;
  return 0;
}


void RigidBody::init_S_matrix(MATRIX& S){
  S.M[0]  = rb_L.Lt; S.M[1]  =-rb_L.Lx; S.M[2]  =-rb_L.Ly; S.M[3]  =-rb_L.Lz;
  S.M[4]  = rb_L.Lx; S.M[5]  = rb_L.Lt; S.M[6]  =-rb_L.Lz; S.M[7]  = rb_L.Ly;
  S.M[8]  = rb_L.Ly; S.M[9]  = rb_L.Lz; S.M[10] = rb_L.Lt; S.M[11] =-rb_L.Lx;
  S.M[12] = rb_L.Lz; S.M[13] =-rb_L.Ly; S.M[14] = rb_L.Lx; S.M[15] = rb_L.Lt;

}

int RigidBody::set_angular_momentum(const VECTOR& l){
  rb_l_e.init(l);                         is_rb_l_e = 1;
  // Update dependent variables
  if(is_rb_A && is_rb_B && is_rb_C) {
    rb_w_e.x = rb_A * rb_l_e.x; rb_w_e.y = rb_B * rb_l_e.y; rb_w_e.z = rb_C * rb_l_e.z;
    is_rb_w_e = 1;
  }
  else if(is_rb_invI_e){  rb_w_e = rb_invI_e * rb_l_e; is_rb_w_e = 1; }

  MATRIX S(4,4); init_S_matrix(S);
  QUATERNION tmp(0.0, rb_l_e.x, rb_l_e.y, rb_l_e.z);
  rb_p_r = 2.0*S*tmp;  is_rb_p_r = 1;

  return 0;
}

int RigidBody::set_angular_momentum(const double& lx, const double& ly, const double& lz){
  rb_l_e.init(lx,ly,lz);                  is_rb_l_e = 1;
  // Update dependent variables
  if(is_rb_A && is_rb_B && is_rb_C) {
    rb_w_e.x = rb_A * rb_l_e.x; rb_w_e.y = rb_B * rb_l_e.y; rb_w_e.z = rb_C * rb_l_e.z;
    is_rb_w_e = 1;
  }
  else if(is_rb_invI_e){  rb_w_e = rb_invI_e * rb_l_e; is_rb_w_e = 1; }
  MATRIX S(4,4); init_S_matrix(S);
  QUATERNION tmp(0.0, rb_l_e.x, rb_l_e.y, rb_l_e.z);
  rb_p_r = 2.0*S*tmp;  is_rb_p_r = 1;

  return 0;
}

int RigidBody::set_angular_velocity(const VECTOR& w){
  rb_w_e.init(w);                          is_rb_w_e = 1;
  // Update dependent variables
  if(is_rb_I_e) { rb_l_e = rb_I_e * rb_w_e; is_rb_l_e = 1;

  MATRIX S(4,4); init_S_matrix(S);
  QUATERNION tmp(0.0, rb_l_e.x, rb_l_e.y, rb_l_e.z);
  rb_p_r = 2.0*S*tmp;  is_rb_p_r = 1;
  }else{ exit(9); }
  return 0;
}

int RigidBody::set_angular_velocity(const double& wx, const double& wy, const double& wz){
  rb_w_e.init(wx,wy,wz);                    is_rb_w_e = 1;
  // Update dependent variables
  if(is_rb_I_e) { rb_l_e = rb_I_e * rb_w_e; is_rb_l_e = 1;

  MATRIX S(4,4); init_S_matrix(S);
  QUATERNION tmp(0.0, rb_l_e.x, rb_l_e.y, rb_l_e.z);
  rb_p_r = 2.0*S*tmp;  is_rb_p_r = 1;
  }
  return 0;
}

int RigidBody::set_quaternion_momentum(const QUATERNION& qm){

  rb_p_r = qm;                     is_rb_p_r = 1;
  MATRIX S(4,4); init_S_matrix(S); S = S.T();
  QUATERNION tmp;tmp = S*qm*0.5;

  rb_l_e.init(tmp.Lx,tmp.Ly,tmp.Lz);  is_rb_l_e = 1;
  rb_w_e.init(rb_A*rb_l_e.x, rb_B*rb_l_e.y, rb_C*rb_l_e.z); is_rb_w_e = 1;

  return 0;
}


}// namespace librigidbody
}// namespace libdyn

