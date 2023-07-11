/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Overlaps.h"
#include "Moments.h"

/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{

/// Here we will just properly wrap Moments - for conveniency, although may be more expensive

/// ----------- 1D --------------------------
double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs, double& dI_dXa, double& dI_dXb,
  vector<double*>& aux,int n_aux 
){
/*******************************************************
 This is basically <G(A)| x |G(B)> 
*******************************************************/
  double dI_dX;                                             
  return gaussian_moment(nxa,alp_a,Xa, 1,0.0,0.0, nxb,alp_b,Xb, is_normalize, is_derivs,dI_dXa,dI_dX,dI_dXb,aux,n_aux);

}// transition_dipole_moment - most genral version



double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs, double& dI_dXa, double& dI_dXb
){

  // Allocate working memory
  int i;
  int n_aux = 20; 
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = transition_dipole_moment(nxa,alp_a,Xa,  nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dXb, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;

}// no external memory


boost::python::list transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs
){

  double dI_dXa, dI_dXb;
  double I = transition_dipole_moment(nxa,alp_a,Xa, nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dXb);

  boost::python::list res;

  res.append(I); 
  if(is_derivs){
    res.append(dI_dXa);
    res.append(dI_dXb);
  }

  return res;
 
}// python version


double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize
){

  double dI_dXa, dI_dXb;
  double res = transition_dipole_moment(nxa,alp_a,Xa, nxb,alp_b,Xb, is_normalize, 0, dI_dXa, dI_dXb);

  return res;

}// with optional normalization

double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb
){

  double res = transition_dipole_moment(nxa,alp_a,Xa, nxb,alp_b,Xb, 1);
  return res;

}// default version








/// --------------3D -----------------------

VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
){
/********************************************************************************************
 This function computes moments:

 <g_a(x,y,z)|  (x, y, z)^T  |g_b(x,y,z)>   - transition dipole moments in 3D-Gaussian basis

  Rether than calling 3D moments, we will be calling 1D functions - for efficiency
********************************************************************************************/

  // 1D overlaps
  double dIx_dXa, dIx_dX, dIx_dXb;
  double dIy_dYa, dIy_dY, dIy_dYb;
  double dIz_dZa, dIz_dZ, dIz_dZb;

  double Ix = gaussian_overlap(nxa, alp_a, Ra.x, nxb, alp_b, Rb.x, is_normalize, is_derivs, dIx_dXa, dIx_dXb, auxd, n_aux);
  double Iy = gaussian_overlap(nya, alp_a, Ra.y, nyb, alp_b, Rb.y, is_normalize, is_derivs, dIy_dYa, dIy_dYb, auxd, n_aux);
  double Iz = gaussian_overlap(nza, alp_a, Ra.z, nzb, alp_b, Rb.z, is_normalize, is_derivs, dIz_dZa, dIz_dZb, auxd, n_aux);


  // 1D transition dipole moments: Mx = <G(A)| x |G(B)> , etc.
  double dMx_dXa, dMx_dX, dMx_dXb;
  double dMy_dYa, dMy_dY, dMy_dYb;
  double dMz_dZa, dMz_dZ, dMz_dZb;

  double Mx = gaussian_moment(nxa, alp_a, Ra.x,  1, 0.0, 0.0,  nxb, alp_b, Rb.x, is_normalize, is_derivs, dMx_dXa, dMx_dX, dMx_dXb, auxd, n_aux);
  double My = gaussian_moment(nya, alp_a, Ra.y,  1, 0.0, 0.0,  nyb, alp_b, Rb.y, is_normalize, is_derivs, dMy_dYa, dMy_dY, dMy_dYb, auxd, n_aux);
  double Mz = gaussian_moment(nza, alp_a, Ra.z,  1, 0.0, 0.0,  nzb, alp_b, Rb.z, is_normalize, is_derivs, dMz_dZa, dMz_dZ, dMz_dZb, auxd, n_aux);


  VECTOR M;
  M.x = Mx * Iy * Iz;
  M.y = Ix * My * Iz;
  M.z = Ix * Iy * Mz;

  dMdA = 0.0;
  dMdB = 0.0;

  if(is_derivs){

    dMdA.xx = dMx_dXa * Iy * Iz;    dMdA.xy = Mx * dIy_dYa * Iz;    dMdA.xz = Mx * Iy * dIz_dZa;
    dMdA.yx = dIx_dXa * My * Iz;    dMdA.yy = Ix * dMy_dYa * Iz;    dMdA.yz = Ix * My * dIz_dZa;
    dMdA.zx = dIx_dXa * Iy * Mz;    dMdA.zy = Ix * dIy_dYa * Mz;    dMdA.zz = Ix * Iy * dMz_dZa;

    dMdB.xx = dMx_dXb * Iy * Iz;    dMdB.xy = Mx * dIy_dYb * Iz;    dMdB.xz = Mx * Iy * dIz_dZb;
    dMdB.yx = dIx_dXb * My * Iz;    dMdB.yy = Ix * dMy_dYb * Iz;    dMdB.yz = Ix * My * dIz_dZb;
    dMdB.zx = dIx_dXb * Iy * Mz;    dMdB.zy = Ix * dIy_dYb * Mz;    dMdB.zz = Ix * Iy * dMz_dZb;

  }

  return M;

}// 3D transition_dipole_moment - most general version


VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
){

  // Allocate working memory
  int i;
  int n_aux = 20; 
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  VECTOR M; M = transition_dipole_moment(nxa,nya,nza,alp_a,Ra,  nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, dMdA, dMdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return M;

}// no external memory



boost::python::list transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize,int is_derivs
){

  MATRIX3x3 dMdA, dMdB;
  boost::python::list res;

  // Do computations
  VECTOR M; M = transition_dipole_moment(nxa,nya,nza,alp_a,Ra,  nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, dMdA, dMdB);

  res.append(M); 
  if(is_derivs){
    res.append(dMdA);
    res.append(dMdB);
  }

  return res;

}// python version


VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize
){

  MATRIX3x3 dMdA, dMdB;

  VECTOR M; M = transition_dipole_moment(nxa,nya,nza,alp_a,Ra,  nxb,nyb,nzb,alp_b,Rb, is_normalize, 0, dMdA, dMdB);
  return M;

}// optional normalization


VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb
){

  VECTOR M; M = transition_dipole_moment(nxa,nya,nza,alp_a,Ra,  nxb,nyb,nzb,alp_b,Rb, 1);
  return M;

}// default version



}//namespace libmolint
}//namespace liblibra
