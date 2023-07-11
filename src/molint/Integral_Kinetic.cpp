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
#include "Integral_Kinetic.h"

/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{


///==================== 1D kinetic energy integrals ==============

double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize,
                        int is_derivs, double& dI_dXa,double& dI_dXb,
                        vector<double*>& aux,int n_aux ){

/****************************************************************************
 This function computes the 1D kinetic energy:

 <g_a(x)| -1/2 * d^2/dx^2 |g_b(x)>  - just the case of 1D Gaussians

*****************************************************************************/

  // Operator = (-1/2 * d^2/dx^2)
  double Ix,w,dIx_dXa,dIx_dXb; 
  dI_dXa = 0.0;
  dI_dXb = 0.0;

  w = (2.0*nxb + 1.0)*alp_b; 
  Ix = w * gaussian_overlap(nxa, alp_a, Xa, nxb, alp_b, Xb, 0, is_derivs, dIx_dXa, dIx_dXb, aux, n_aux);
  dI_dXa = w * dIx_dXa;
  dI_dXb = w * dIx_dXb;


  w = -2.0*alp_b*alp_b;
  Ix += w * gaussian_overlap(nxa, alp_a, Xa, nxb+2, alp_b, Xb, 0, is_derivs, dIx_dXa, dIx_dXb, aux, n_aux);
  dI_dXa += w * dIx_dXa;
  dI_dXb += w * dIx_dXb;


  if(nxb>=2){ 
    w = -0.5*nxb*(nxb-1.0);
    Ix += w * gaussian_overlap(nxa, alp_a, Xa, nxb-2, alp_b, Xb, 0, is_derivs, dIx_dXa, dIx_dXb, aux, n_aux);
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

}


double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize,
                        int is_derivs, double& dI_dXa,double& dI_dXb
                       ){

  // Allocate working memory
  int i;
  int n_aux = 20; //nxa+nxb+1;
  vector<double*> auxd(5);
  for(i=0;i<5;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = kinetic_integral(nxa,alp_a,Xa, nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dXb, auxd, n_aux);

  // Clean working memory
  for(i=0;i<5;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;

}// kinetic_integral


boost::python::list kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                                     int is_normalize, int is_derivs ){
  double dI_dXa, dI_dXb;
  double I = kinetic_integral(nxa,alp_a,Xa, nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dXb);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dI_dXa);
    res.append(dI_dXb);
  }

  return res;
 
}



double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize
                       ){
  double dI_dxa, dI_dxb;
  double res = kinetic_integral(nxa,alp_a,Xa, nxb,alp_b,Xb, is_normalize, 0, dI_dxa, dI_dxb);
  return res;

}// kinetic_integral

double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb){

  double res = kinetic_integral(nxa,alp_a,Xa, nxb,alp_b,Xb, 1);
  return res;

}// kinetic_integral






///========================== 3D kinetic energy integrals ========================



double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
                        int is_derivs,
                        VECTOR& dIdA, VECTOR& dIdB,
                        vector<double*>& auxd,int n_aux
                       ){
/****************************************************************************
 This function computes the 1D kinetic energy:

 <g_a(Ra)| -1/2 * ( d^2/dx^2 + d^2/dy^2 + d^2/dz^2 ) |g_b(Rb)>  - just the case of 1D Gaussians

// is_normalize - controls if the (original) Gaussians should be normalized:
//   is_normalize = 0 - use Gaussians as they are
//   is_normalize = 1 - use normalize the overlap as if the Gaussians are normalized

*****************************************************************************/


  double dIx_dXa, dIx_dXb, dIy_dYa, dIy_dYb, dIz_dZa, dIz_dZb;

  double Ix = gaussian_overlap(nxa, alp_a, Ra.x, nxb, alp_b, Rb.x, is_normalize, is_derivs, dIx_dXa, dIx_dXb, auxd, n_aux);
  double Iy = gaussian_overlap(nya, alp_a, Ra.y, nyb, alp_b, Rb.y, is_normalize, is_derivs, dIy_dYa, dIy_dYb, auxd, n_aux);
  double Iz = gaussian_overlap(nza, alp_a, Ra.z, nzb, alp_b, Rb.z, is_normalize, is_derivs, dIz_dZa, dIz_dZb, auxd, n_aux);

  double dKx_dXa, dKx_dXb, dKy_dYa, dKy_dYb, dKz_dZa, dKz_dZb;

  double Kx = kinetic_integral(nxa, alp_a, Ra.x, nxb, alp_b, Rb.x, is_normalize, is_derivs, dKx_dXa, dKx_dXb, auxd, n_aux);
  double Ky = kinetic_integral(nya, alp_a, Ra.y, nyb, alp_b, Rb.y, is_normalize, is_derivs, dKy_dYa, dKy_dYb, auxd, n_aux);
  double Kz = kinetic_integral(nza, alp_a, Ra.z, nzb, alp_b, Rb.z, is_normalize, is_derivs, dKz_dZa, dKz_dZb, auxd, n_aux);


  double I = (Kx * Iy * Iz + Ix * Ky * Iz + Ix * Iy * Kz);

  dIdA = 0.0;
  dIdB = 0.0;

  if(is_derivs){
    dIdA.x = (dKx_dXa * Iy * Iz + dIx_dXa * Ky * Iz + dIx_dXa * Iy * Kz);
    dIdA.y = (Kx * dIy_dYa * Iz + Ix * dKy_dYa * Iz + Ix * dIy_dYa * Kz);
    dIdA.z = (Kx * Iy * dIz_dZa + Ix * Ky * dIz_dZa + Ix * Iy * dKz_dZa);
    dIdB = -dIdA;
  }

  return I;

}// kinetic_integral


double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, int is_derivs,
                        VECTOR& dIdA, VECTOR& dIdB
                       ){
  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(5);
  for(i=0;i<5;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = kinetic_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<5;i++){ delete [] auxd[i]; }  
  auxd.clear();
 

  return res;

}

boost::python::list kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                     int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                     int is_normalize, int is_derivs
                                    ){
  VECTOR dIdA, dIdB;
  double I = kinetic_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, dIdA, dIdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dIdA);
    res.append(dIdB);
  }

  return res;
 
}


double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize
                       ){

  VECTOR dIdA, dIdB;
  double res = kinetic_integral(nxa,nya,nza,alp_a,Ra,nxb,nyb,nzb,alp_b,Rb, is_normalize, 0,dIdA,dIdB);
  return res;
}

double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                       ){
  double res = kinetic_integral(nxa,nya,nza,alp_a,Ra,nxb,nyb,nzb,alp_b,Rb, 1);
  return res;
}





}// namespace libmolint
}// namespace liblinalg



