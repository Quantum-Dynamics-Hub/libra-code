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
  \file RigidBody_methods5_3.cpp
  \brief The file implements the SQUISH propagation method based on quaternion rotations. This is a symplectic scheme.

*/

#include "RigidBody.h"

/// liblibra namespace
namespace liblibra{


/// librigidbody
namespace librigidbody{


void RigidBody::rotate_no_squish(int k,double dt){
/** 
  \brief A single SQUISH rotation.
  \param[in] k The selectron of the plane in which the rotation occurs
  \param[in] dt Duration of the propagation for this step 
*/
  double zeta,cos_zeta,sin_zeta, rot_const;
  MATRIX P(4,4);

  if(k==1) { P = *P1; rot_const = rb_A;}
  else if(k==2) { P = *P2; rot_const = rb_B;}
  else if(k==3) { P = *P3; rot_const = rb_C;}


  zeta = 0.25*rot_const*dt*dot_prod(rb_p_r,P*rb_L);
  cos_zeta = cos(zeta);
  sin_zeta = sin(zeta);

  rb_L = cos_zeta * rb_L + sin_zeta * P*rb_L;
  rb_p_r = cos_zeta * rb_p_r + sin_zeta * P*rb_p_r;

}

void RigidBody::propagate_no_squish(double t){
/**
  \brief The eventual SQUISH propagation scheme
  \param[in] t The propagation timestep

  Implemented according to:
  Ikeguchi, M. Partial Rigid-Body Dynamics in NPT, NPAT and NP$\gamma$T Ensembles for Proteins and Membranes. J. Comput. Chem. 2004, 25 (4), 529–541.

*/

  double t_half = 0.5 * t;

  rotate_no_squish(3,t_half);
  rotate_no_squish(2,t_half);
  rotate_no_squish(1,t);
  rotate_no_squish(2,t_half);
  rotate_no_squish(3,t_half);

  set_orientation(rb_L);
  set_quaternion_momentum(rb_p_r);

}

}// namespace librigidbody
}// liblibra
