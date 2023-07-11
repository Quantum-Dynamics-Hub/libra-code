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
#include "Integral_Derivative_Couplings.h"


/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;


namespace libmolint{


///==================== 1D derivative coupling integrals ==============

double derivative_coupling_integral(
  int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs, double& dI_dXa,double& dI_dXb,
  vector<double*>& aux,int n_aux
){

/****************************************************************************
 This function computes the 1D kinetic energy:

 <g_a(x)| d/dxb |g_b(x)>  - 1D

 Note that <g_a(x)| d/dxa |g_b(x)>  = 0.0

 Also note that this is not the same as electronic nubla (electronic momentum or kinetic energy) operator
 Here, the operator actually acts on the "nuclear" degrees of freedom
 
 Hey, we still may need derivatives of the derivative coupling vector
 (which will make it to rank-2 tensor: 3x3 matrix in 3D case or 1x1 in 1D case)

*****************************************************************************/

  // Operator = (d/dx)
  double Ix,w,dIx_dXa,dIx_dXb; 
  dI_dXa = 0.0;
  dI_dXb = 0.0;


  w = 2.0*alp_b;
  Ix += w * gaussian_overlap(nxa, alp_a, Xa, nxb+1, alp_b, Xb, 0, is_derivs, dIx_dXa, dIx_dXb, aux, n_aux);
  dI_dXa += w * dIx_dXa;
  dI_dXb += w * dIx_dXb;


  if(nxb>=1){ 
    w = -nxb;
    Ix += w * gaussian_overlap(nxa, alp_a, Xa, nxb-1, alp_b, Xb, 0, is_derivs, dIx_dXa, dIx_dXb, aux, n_aux);
    dI_dXa += w * dIx_dXa;
    dI_dXb += w * dIx_dXb;

  }

  // In case we need to normalize initial Gaussians
  if(is_normalize){

    double nrm = gaussian_normalization_factor(nxa,alp_a) * gaussian_normalization_factor(nxb,alp_b);
    Ix *= nrm;
    dI_dXa *= nrm; 
    dI_dXb *= nrm; 

  }



  return Ix;

}// 1D case - most general case


double derivative_coupling_integral(
  int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs, double& dI_dXa,double& dI_dXb
){


  // Allocate working memory
  int i;
  int n_aux = 20; 
  vector<double*> auxd(5);
  for(i=0;i<5;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = derivative_coupling_integral(nxa,alp_a,Xa, nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dXb, auxd, n_aux);

  // Clean working memory
  for(i=0;i<5;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;

}// version without external memory



boost::python::list derivative_coupling_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                                                 int is_normalize, int is_derivs ){
  double dI_dXa, dI_dXb;
  double I = derivative_coupling_integral(nxa,alp_a,Xa, nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dXb);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dI_dXa);
    res.append(dI_dXb);
  }

  return res;
 
}



double derivative_coupling_integral
( int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
  int is_normalize
){

  double dI_dxa, dI_dxb;
  double res = derivative_coupling_integral(nxa,alp_a,Xa, nxb,alp_b,Xb, is_normalize, 0, dI_dxa, dI_dxb);
  return res;

}// derivative_coupling_integral

double derivative_coupling_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb){

  double res = derivative_coupling_integral(nxa,alp_a,Xa, nxb,alp_b,Xb, 1);
  return res;

}// derivative_coupling_integral



///========================== 3D derivative coupling integrals ========================


VECTOR derivative_coupling_integral
( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
  int is_normalize, int is_derivs,
  MATRIX3x3& dDdA, MATRIX3x3& dDdB,
  vector<double*>& auxd,int n_aux
){
/****************************************************************************
 This function computes 3D derivative couplings

 D = <g_a(Ra)| ( d/dxb , d/dyb , d/dzb)^T |g_b(Rb)>  

 note that  <g_a(Ra)| ( d/dxa , d/dya , d/dza)^T |g_b(Rb)>  = 0.0

 but the derivatives of the D = DA with both A and B coordinates are not zero

*****************************************************************************/


  double dIx_dXa, dIx_dXb, dIy_dYa, dIy_dYb, dIz_dZa, dIz_dZb;
  double Ix = gaussian_overlap(nxa, alp_a, Ra.x, nxb, alp_b, Rb.x, is_normalize, is_derivs, dIx_dXa, dIx_dXb, auxd, n_aux);
  double Iy = gaussian_overlap(nya, alp_a, Ra.y, nyb, alp_b, Rb.y, is_normalize, is_derivs, dIy_dYa, dIy_dYb, auxd, n_aux);
  double Iz = gaussian_overlap(nza, alp_a, Ra.z, nzb, alp_b, Rb.z, is_normalize, is_derivs, dIz_dZa, dIz_dZb, auxd, n_aux);

  double dDx_dXa, dDx_dXb, dDy_dYa, dDy_dYb, dDz_dZa, dDz_dZb;
  double Dx = derivative_coupling_integral(nxa, alp_a, Ra.x, nxb, alp_b, Rb.x, is_normalize, is_derivs, dDx_dXa, dDx_dXb, auxd, n_aux);
  double Dy = derivative_coupling_integral(nya, alp_a, Ra.y, nyb, alp_b, Rb.y, is_normalize, is_derivs, dDy_dYa, dDy_dYb, auxd, n_aux);
  double Dz = derivative_coupling_integral(nza, alp_a, Ra.z, nzb, alp_b, Rb.z, is_normalize, is_derivs, dDz_dZa, dDz_dZb, auxd, n_aux);


  VECTOR D;
  D.x = Dx * Iy * Iz;
  D.y = Ix * Dy * Iz;
  D.z = Ix * Iy * Dz;

  dDdA = 0.0;
  dDdB = 0.0;

  if(is_derivs){

    dDdA.xx = dDx_dXa * Iy * Iz;    dDdA.xy = Dx * dIy_dYa * Iz;    dDdA.xz = Dx * Iy * dIz_dZa;
    dDdA.yx = dIx_dXa * Dy * Iz;    dDdA.yy = Ix * dDy_dYa * Iz;    dDdA.yz = Ix * Dy * dIz_dZa;
    dDdA.zx = dIx_dXa * Iy * Dz;    dDdA.zy = Ix * dIy_dYa * Dz;    dDdA.zz = Ix * Iy * dDz_dZa;

    dDdB.xx = dDx_dXb * Iy * Iz;    dDdB.xy = Dx * dIy_dYb * Iz;    dDdB.xz = Dx * Iy * dIz_dZb;
    dDdB.yx = dIx_dXb * Dy * Iz;    dDdB.yy = Ix * dDy_dYb * Iz;    dDdB.yz = Ix * Dy * dIz_dZb;
    dDdB.zx = dIx_dXb * Iy * Dz;    dDdB.zy = Ix * dIy_dYb * Dz;    dDdB.zz = Ix * Iy * dDz_dZb;

  }

  return D;

}// derivative_coupling_integral



VECTOR derivative_coupling_integral
( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
  int is_normalize, int is_derivs,
  MATRIX3x3& dDdA, MATRIX3x3& dDdB
){


  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(5);
  for(i=0;i<5;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  VECTOR D; D = derivative_coupling_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, dDdA, dDdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<5;i++){ delete [] auxd[i]; }  
  auxd.clear();
 

  return D;


}// version without external memory


boost::python::list derivative_coupling_integral
( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
  int is_normalize, int is_derivs
){

  MATRIX3x3 dDdA, dDdB;
  VECTOR D; D = derivative_coupling_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, dDdA, dDdB);

  boost::python::list res;
  res.append(D);
 
  if(is_derivs){
    res.append(dDdA);
    res.append(dDdB);
  }

  return res;
 
}// version for python


VECTOR derivative_coupling_integral
( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
  int is_normalize
){

  MATRIX3x3 dDdA, dDdB;
  VECTOR D; D = derivative_coupling_integral(nxa,nya,nza,alp_a,Ra,nxb,nyb,nzb,alp_b,Rb, is_normalize, 0, dDdA, dDdB);
  return D;
}

VECTOR derivative_coupling_integral
( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
){

  VECTOR D; D = derivative_coupling_integral(nxa,nya,nza,alp_a,Ra,nxb,nyb,nzb,alp_b,Rb, 1);
  return D;

}



}// namespace libmolint
}// namespace liblibra


