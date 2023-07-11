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
  \file RigidBody_methods5_5.cpp
  \brief The file implements the Terec and qTerec propagation method based on recursive Taylor series expansion
  
  See more in:
  Akimov, A. V.; Kolomeisky, A. B. Recursive Taylor Series Expansion Method for Rigid-Body Molecular Dynamics. J. Chem. Theory Comput. 2011, 7, 3062–3071.

*/

#include "RigidBody.h"

/// liblibra namespace
namespace liblibra{


/// librigidbody namespace
namespace librigidbody{


void RigidBody::initialize_terec(int series_expansion){
/**
  \brief This is an auxiliary function that initializes auxiliary variables and parameters
  \param[in] series_expansion The length of the Taylor expansion chosen for the integration
*/

  //--------------- Calculate expansion coefficients ---------------------
  int i,n,k,indx,Coeffs_size;

  SERIES_EXPANSION = series_expansion;
  Coeffs_size = (SERIES_EXPANSION*(1+SERIES_EXPANSION)/2);

  if(Coeffs!=NULL){ delete [] Coeffs; }
  Coeffs = new double[Coeffs_size];

  for(i=0;i<Coeffs_size;i++){  Coeffs[i] = 0.0; }
  for(n=2;n<SERIES_EXPANSION;n++){
    for(k=0;k<n;k++){
      indx = ((1+n)*n/2)-1+k;
      Coeffs[indx] = BINOM(k,n-1);
    }
  }

}

void RigidBody::propagate_terec(double t){
/**
  \brief Terec propagation scheme
  \param[in] t Time of propagation

  This is my TEylor RECursive algorithm for free rigid body
  problem. See more in:
  Akimov, A. V.; Kolomeisky, A. B. Recursive Taylor Series Expansion Method for Rigid-Body Molecular Dynamics. J. Chem. Theory Comput. 2011, 7, 3062–3071.

  The following convention is used:

              u1[0].x  u2[0].x  u3[0].x
  A_I_to_e =  u1[0].y  u2[0].y  u3[0].y
              u1[0].z  u2[0].z  u3[0].z
*/
  int n;
  VECTOR *u1,*u2,*u3,*deriv;
  double alpha,beta,gamma,dtn,bin;
  int indx;

  u1 = new VECTOR[SERIES_EXPANSION];
  u2 = new VECTOR[SERIES_EXPANSION];
  u3 = new VECTOR[SERIES_EXPANSION];
  deriv = new VECTOR[SERIES_EXPANSION];

  rb_A_I_to_e.get_vectors(u1[0],u2[0],u3[0]);
  deriv[0] = rb_l_e;

  VECTOR U1(u1[0]);
  VECTOR U2(u2[0]);
  VECTOR U3(u3[0]);

  //----------------- Orientation matrix representation ---------------------

  alpha = (rb_C-rb_B);
  beta  = (rb_A-rb_C);
  gamma = (rb_B-rb_A);

  deriv[1].x = alpha*rb_l_e.y*rb_l_e.z;
  deriv[1].y = beta *rb_l_e.x*rb_l_e.z;
  deriv[1].z = gamma*rb_l_e.x*rb_l_e.y;

  VECTOR omega(rb_A*rb_l_e.x,rb_B*rb_l_e.y,rb_C*rb_l_e.z);

  u1[1].x = omega.z*u1[0].y - omega.y*u1[0].z;
  u1[1].y = omega.x*u1[0].z - omega.z*u1[0].x;
  u1[1].z = omega.y*u1[0].x - omega.x*u1[0].y;

  u2[1].x = omega.z*u2[0].y - omega.y*u2[0].z;
  u2[1].y = omega.x*u2[0].z - omega.z*u2[0].x;
  u2[1].z = omega.y*u2[0].x - omega.x*u2[0].y;

  u3[1].x = omega.z*u3[0].y - omega.y*u3[0].z;
  u3[1].y = omega.x*u3[0].z - omega.z*u3[0].x;
  u3[1].z = omega.y*u3[0].x - omega.x*u3[0].y;

  U1 += u1[1]*t;
  U2 += u2[1]*t;
  U3 += u3[1]*t;
  rb_l_e += deriv[1]*t;

  dtn = t;

  // Clear vectors
  for(n=2;n<SERIES_EXPANSION;n++){
    deriv[n] = 0.0; u1[n] = 0.0; u2[n] = 0.0; u3[n] = 0.0;
  }
  for(n=2;n<SERIES_EXPANSION;n++){
    for(int k=0;k<n;k++){
    //bin = BINOM(k,n-1);
    //-----------------------------
    indx = ((1+n)*n/2)-1+k;
    bin  = Coeffs[indx];
    //-----------------------------
    deriv[n].x+= bin*deriv[k].y*deriv[n-k-1].z;
    deriv[n].y+= bin*deriv[k].x*deriv[n-k-1].z;
    deriv[n].z+= bin*deriv[k].x*deriv[n-k-1].y;

    u1[n].x += bin*(rb_C*deriv[k].z*u1[n-k-1].y - rb_B*deriv[k].y*u1[n-k-1].z);
    u1[n].y += bin*(rb_A*deriv[k].x*u1[n-k-1].z - rb_C*deriv[k].z*u1[n-k-1].x);
    u1[n].z += bin*(rb_B*deriv[k].y*u1[n-k-1].x - rb_A*deriv[k].x*u1[n-k-1].y);

    u2[n].x += bin*(rb_C*deriv[k].z*u2[n-k-1].y - rb_B*deriv[k].y*u2[n-k-1].z);
    u2[n].y += bin*(rb_A*deriv[k].x*u2[n-k-1].z - rb_C*deriv[k].z*u2[n-k-1].x);
    u2[n].z += bin*(rb_B*deriv[k].y*u2[n-k-1].x - rb_A*deriv[k].x*u2[n-k-1].y);

    u3[n].x += bin*(rb_C*deriv[k].z*u3[n-k-1].y - rb_B*deriv[k].y*u3[n-k-1].z);
    u3[n].y += bin*(rb_A*deriv[k].x*u3[n-k-1].z - rb_C*deriv[k].z*u3[n-k-1].x);
    u3[n].z += bin*(rb_B*deriv[k].y*u3[n-k-1].x - rb_A*deriv[k].x*u3[n-k-1].y);

    } // for k

    deriv[n].x *= alpha;
    deriv[n].y *= beta;
    deriv[n].z *= gamma;
    dtn *= t/((double)n);

    rb_l_e += dtn*deriv[n];

    U1   += dtn*u1[n];
    U2   += dtn*u2[n];
    U3   += dtn*u3[n];

  } // for n

  U1.normalize();
  U2.normalize();
  U3.normalize();

  rb_A_I_to_e.init(U1,U2,U3);

  // Update dependent variables
  set_angular_momentum(rb_l_e);
  set_orientation(rb_A_I_to_e);

  // Clear the temporary memory
  delete [] u1;
  delete [] u2;
  delete [] u3;
  delete [] deriv;

}

void RigidBody::propagate_qterec(double t){
/**
  \brief Terec propagation scheme - quaternion version
  \param[in] t Time of propagation

  This is my TEylor RECursive algorithm for free rigid body
  problem. See more in:
  Akimov, A. V.; Kolomeisky, A. B. Recursive Taylor Series Expansion Method for Rigid-Body Molecular Dynamics. J. Chem. Theory Comput. 2011, 7, 3062–3071.

  The following convention is used:

              u1[0].x  u2[0].x  u3[0].x
  A_I_to_e =  u1[0].y  u2[0].y  u3[0].y
              u1[0].z  u2[0].z  u3[0].z
*/

  int n;
  QUATERNION *derivq;
  VECTOR *deriv;
  double alpha,beta,gamma,dtn,bin;
  int indx;

  derivq = new QUATERNION[SERIES_EXPANSION];
  deriv = new VECTOR[SERIES_EXPANSION];

  deriv[0] = rb_l_e;
  QUATERNION L(rb_L.Lt,rb_L.Lx,rb_L.Ly,rb_L.Lz);
  derivq[0] = L;

  //----------------- Quaternion representation ---------------------
  alpha = (rb_C-rb_B);
  beta  = (rb_A-rb_C);
  gamma = (rb_B-rb_A);

  deriv[1].x = alpha*rb_l_e.y*rb_l_e.z;
  deriv[1].y = beta *rb_l_e.x*rb_l_e.z;
  deriv[1].z = gamma*rb_l_e.x*rb_l_e.y;

  VECTOR omega(rb_A*rb_l_e.x,rb_B*rb_l_e.y,rb_C*rb_l_e.z);

  derivq[1].Lt = 0.5*( -L.Lx*omega.x - L.Ly*omega.y - L.Lz*omega.z);
  derivq[1].Lx = 0.5*(  L.Lt*omega.x - L.Lz*omega.y + L.Ly*omega.z);
  derivq[1].Ly = 0.5*(  L.Lz*omega.x + L.Lt*omega.y - L.Lx*omega.z);
  derivq[1].Lz = 0.5*( -L.Ly*omega.x + L.Lx*omega.y + L.Lt*omega.z);
  L += derivq[1]*t;

  rb_l_e += deriv[1]*t;

  dtn = t;

  // Clear vectors
  for(n=2;n<SERIES_EXPANSION;n++){
    deriv[n] = 0.0; derivq[n] = 0.0;
  }
  for(n=2;n<SERIES_EXPANSION;n++){
    for(int k=0;k<n;k++){
    //bin = BINOM(k,n-1);
    //-----------------------------
    indx = ((1+n)*n/2)-1+k;
    bin  = Coeffs[indx];
    //-----------------------------
    deriv[n].x+= bin*deriv[k].y*deriv[n-k-1].z;
    deriv[n].y+= bin*deriv[k].x*deriv[n-k-1].z;
    deriv[n].z+= bin*deriv[k].x*deriv[n-k-1].y;

    derivq[n].Lt += bin*0.5*( -derivq[k].Lx*deriv[n-k-1].x*rb_A - derivq[k].Ly*deriv[n-k-1].y*rb_B - derivq[k].Lz*deriv[n-k-1].z*rb_C );
    derivq[n].Lx += bin*0.5*(  derivq[k].Lt*deriv[n-k-1].x*rb_A - derivq[k].Lz*deriv[n-k-1].y*rb_B + derivq[k].Ly*deriv[n-k-1].z*rb_C );
    derivq[n].Ly += bin*0.5*(  derivq[k].Lz*deriv[n-k-1].x*rb_A + derivq[k].Lt*deriv[n-k-1].y*rb_B - derivq[k].Lx*deriv[n-k-1].z*rb_C );
    derivq[n].Lz += bin*0.5*( -derivq[k].Ly*deriv[n-k-1].x*rb_A + derivq[k].Lx*deriv[n-k-1].y*rb_B + derivq[k].Lt*deriv[n-k-1].z*rb_C );

    } // for k

    deriv[n].x *= alpha;
    deriv[n].y *= beta;
    deriv[n].z *= gamma;
    dtn *= t/((double)n);

    rb_l_e += dtn*deriv[n];
    L += dtn*derivq[n];

  } // for n

  L.normalize();

  rb_L = L;
  set_angular_momentum(rb_l_e);
  set_orientation(rb_L);

  // Clear the temporary memory
  delete [] derivq;
  delete [] deriv;

}


}// namespace librigidbody
}// liblibra
