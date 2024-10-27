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

#include "Overlaps.h"
#include "Moments.h"

/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{



double gaussian_moment_ref(int nxc, double alp_c, double Xc,
                           int nxa,double alp_a, double Xa,
                           int nxb,double alp_b, double Xb
                      ){
/****************************************************************************
 This function computes moments:

 <g_a(x)|  (x-Xc)^nxc * exp(-alp_c*(x-Xc)^2)  |g_b(x)> = <g_a | g_p> 

 Where g_b(x) = (x-Xb)^nxb * exp(-alp_b*(x-Xb)^2)

 See the revised derivations in the Moments.docx / Moments.pdf

****************************************************************************/

  double gamma_cb = alp_c + alp_b;
  double X_cb = (alp_c * Xc + alp_b * Xb)/gamma_cb;
  double dX = X_cb - Xc;
  double res = 0.0;
  double pwk = 1.0;
  for(int k=0;k<=nxc;k++){
    double bink = BINOM(k,nxc);
    res += bink * pwk * gaussian_overlap_ref(nxa, alp_a, Xa, (nxc-k), gamma_cb, X_cb );
    pwk *= dX;
  }// for i
  res *= exp(-alp_c * alp_b * (Xc - Xb) * (Xc - Xb)/gamma_cb);

  return res;


}


double gaussian_moment(int nxa, double alp_a, double Xa, 
                       int nxc, double alp_c, double Xc, 
                       int nxb, double alp_b, double Xb,
                       int is_normalize, int is_derivs, 
                       double& dI_dXa, double& dI_dXc, double& dI_dXb,
                       vector<double*>& aux,int n_aux ){
/****************************************************************************

****************************************************************************/


  int i,k;
  double gamma_cb = alp_c + alp_b;
  double Xcb = (alp_c * Xc + alp_b * Xb )/gamma_cb;
  double dX = Xcb - Xc;
  double dx = Xc - Xb;
  double pref = exp( -(alp_c * alp_b / gamma_cb) * dx * dx );

  // Aliaces: since places [0,1,2,3] will be used in the gaussian_overlap function, 
  // we reserve next spots:

  double Iabc = 0.0;
  dI_dXa = 0.0;
  dI_dXb = 0.0;
  dI_dXc = 0.0;
  double pwk = 1.0;
  double pwk_prev = 0.0;
  // Compute Gaussian integrals
  for(k=0; k<=nxc; k++){
    double bink = BINOM(k,nxc);
    double g, dSdXa, dSdXcb; 
    //============= Function itself ===============
    // In this call we do not do normalization since we are working with temporary ("virtual") Gaussians
    // normalization will be applied later to Ga an Gb
    g = gaussian_overlap(nxa, alp_a, Xa,   (nxc-k),  gamma_cb, Xcb,  0, is_derivs, dSdXa, dSdXcb, aux, n_aux );
    Iabc += bink * pwk * g;

    //============ Derivatives =======================
    if(is_derivs){ 

      dI_dXa += bink * pwk * dSdXa;
      dI_dXc += bink * (k * pwk_prev * (alp_c/gamma_cb - 1.0) * g + pwk * (alp_c/gamma_cb) * dSdXcb );
      dI_dXb += bink * (k * pwk_prev * (alp_b/gamma_cb) * g       + pwk * (alp_b/gamma_cb) * dSdXcb );

    }// if is_derivs
   
    pwk_prev = pwk;
    pwk *= dX;
  }// for k

  Iabc *= pref;

  if(is_derivs){
    dI_dXa *= pref; 
    dI_dXc *= pref; dI_dXc -= 2.0*(alp_c * alp_b / gamma_cb) * dx * Iabc; 
    dI_dXb *= pref; dI_dXc += 2.0*(alp_c * alp_b / gamma_cb) * dx * Iabc;
  }

  // In case we need to normalize initial Gaussians
  if(is_normalize){
    double nrm = gaussian_normalization_factor(nxa,alp_a) * gaussian_normalization_factor(nxb,alp_b);
    Iabc *= nrm;
    dI_dXa *= nrm; 
    dI_dXb *= nrm; 
    dI_dXc *= nrm;
  }

  return Iabc;

}// gaussian_moment - very general version



double gaussian_moment(int nxa, double alp_a, double Xa, 
                       int nxc, double alp_c, double Xc, 
                       int nxb, double alp_b, double Xb,
                       int is_normalize, int is_derivs, 
                       double& dI_dXa, double& dI_dXc,double& dI_dXb
                      ){

  // Allocate working memory
  int i;
  int n_aux = 20; //nxa+nxb+1;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_moment(nxa, alp_a, Xa,  nxc, alp_c,Xc,  nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dXc, dI_dXb, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;


}// gaussian_moment - without external memory allocation


boost::python::list gaussian_moment(int nxa, double alp_a, double Xa, 
                                    int nxc, double alp_c, double Xc, 
                                    int nxb, double alp_b, double Xb,
                                    int is_normalize, int is_derivs ){
  double dI_dXa, dI_dXb, dI_dXc;
  double I = gaussian_moment(nxa,alp_a,Xa, nxc,alp_c,Xc, nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dXc, dI_dXb);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dI_dXa);
    res.append(dI_dXc);
    res.append(dI_dXb);
  }

  return res;
 
}// gaussian_moment - boost::python version



double gaussian_moment(int nxa, double alp_a, double Xa, 
                       int nxc, double alp_c, double Xc, 
                       int nxb, double alp_b, double Xb,
                       int is_normalize
                      ){

  double dI_dXa, dI_dXc, dI_dXb;
  double res = gaussian_moment(nxa,alp_a,Xa, nxc,alp_c,Xc, nxb,alp_b,Xb, is_normalize, 0, dI_dXa, dI_dXc, dI_dXb);
  return res;


}// gaussian_moment - without derivatives


double gaussian_moment(int nxa, double alp_a, double Xa, 
                       int nxc, double alp_c, double Xc, 
                       int nxb, double alp_b, double Xb ){

  double res = gaussian_moment(nxa,alp_a,Xa, nxc,alp_c,Xc, nxb,alp_b,Xb, 1);
  return res;


}// gaussian_moment - default version: no derivatives, with normalization







double gaussian_moment(int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                       int nx, int ny,  int nz,  double alp, const VECTOR& R,
                       int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                       int is_normalize, 
                       int is_derivs, VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB,
                       vector<double*>& auxd,int n_aux
                      ){
/********************************************************************************************
 This function computes moments:

 <g_a(x)|  (x-X)^nx * (x-Y)^ny * (x-Z)^nz * exp(-alp*(x-X)^2)  |g_b(x)>   - 3D verions of the above

********************************************************************************************/

  double dIx_dXa, dIx_dX, dIx_dXb;
  double dIy_dYa, dIy_dY, dIy_dYb;
  double dIz_dZa, dIz_dZ, dIz_dZb;

  double Ix = gaussian_moment(nxa, alp_a, Ra.x,  nx, alp, R.x,  nxb, alp_b, Rb.x, is_normalize, is_derivs, dIx_dXa, dIx_dX, dIx_dXb, auxd, n_aux);
  double Iy = gaussian_moment(nya, alp_a, Ra.y,  ny, alp, R.y,  nyb, alp_b, Rb.y, is_normalize, is_derivs, dIy_dYa, dIy_dY, dIy_dYb, auxd, n_aux);
  double Iz = gaussian_moment(nza, alp_a, Ra.z,  nz, alp, R.z,  nzb, alp_b, Rb.z, is_normalize, is_derivs, dIz_dZa, dIz_dZ, dIz_dZb, auxd, n_aux);

  double I = Ix * Iy * Iz;

  dIdA = 0.0;
  dIdR = 0.0;
  dIdB = 0.0;

  if(is_derivs){
    dIdA.x = dIx_dXa * Iy * Iz;
    dIdA.y = Ix * dIy_dYa * Iz;
    dIdA.z = Ix * Iy * dIz_dZa;

    dIdR.x = dIx_dX  * Iy * Iz;
    dIdR.y = Ix * dIy_dY  * Iz;
    dIdR.z = Ix * Iy * dIz_dZ ;

    dIdB.x = dIx_dXb * Iy * Iz;
    dIdB.y = Ix * dIy_dYb * Iz;
    dIdB.z = Ix * Iy * dIz_dZb;
  }


  return I;

}// gaussian_moment - very general 3D version


double gaussian_moment(int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                       int nx, int ny,  int nz,  double alp, const VECTOR& R,
                       int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                       int is_normalize, 
                       int is_derivs, VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB
                      ){
  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_moment(nxa,nya,nza,alp_a,Ra, nx,ny,nz,alp,R, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, dIdA, dIdR, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 

  return res;

}// 3D version without external memory

boost::python::list gaussian_moment(int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                                    int nx, int ny,  int nz,  double alp, const VECTOR& R,
                                    int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                                    int is_normalize, int is_derivs
                                    ){
  VECTOR dIdA, dIdR, dIdB;
  double I = gaussian_moment(nxa,nya,nza,alp_a,Ra, nx,ny,nz,alp,R, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, dIdA, dIdR, dIdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dIdA);
    res.append(dIdR);
    res.append(dIdB);
  }

  return res;
 
}// 3D version for python


double gaussian_moment(int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                       int nx, int ny,  int nz,  double alp, const VECTOR& R,
                       int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                       int is_normalize
                      ){

  VECTOR dIdA, dIdR, dIdB;
  double res = gaussian_moment(nxa,nya,nza,alp_a,Ra, nx,ny,nz,alp,R, nxb,nyb,nzb,alp_b,Rb, is_normalize, 0,dIdA,dIdR,dIdB);
  return res;

}// 3D version without derivatives

double gaussian_moment(int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                       int nx, int ny,  int nz,  double alp, const VECTOR& R,
                       int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb
                      ){
  double res = gaussian_moment(nxa,nya,nza,alp_a,Ra, nx,ny,nz,alp,R, nxb,nyb,nzb,alp_b,Rb, 1);
  return res;

}// 3D - default version with normalization and no derivatives





}// namespace libmolint
}// namespace liblibra

