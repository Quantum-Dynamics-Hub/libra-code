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



double gaussian_moment_ref(int nx, double alp, double X,
                           int nxa,double alp_a, double Xa,
                           int nxb,double alp_b, double Xb
                      ){
/****************************************************************************
 This function computes moments:

 <g_a(x)|  (x-X)^nx * exp(-alp*(x-X)^2)  |g_b(x)> = <g_a | g_p>  - note, here we deal with 1D Gaussians

 We want to express the middle expression in terms of the right Gaussian:

 Remember also:

 exp(-alp*(x-X)^2) * exp(-alp_b*(x-Xb)^2) = exp(-alp*alp_b*(X-Xb)^2/gamma) * exp(-gamma*(x-Xp)^2)

 where Xp = (alp*X + alp_b * Xb)/gamma

 gamma = alp + alp_b

 So:
                                     nx
 (x-X)^nx = (x - Xp + Xp - X)^nx =  sum  [  C_nx^i * (x-Xp)^i * (Xp-X)^(nx-i) ]
                                   i = 0

                                        nxb
 (x-Xb)^nxb = (x - Xp + Xp - Xb)^nxa =  sum  [  C_nxb^j * (x-Xp)^j * (Xp-Xb)^(nxb-j) ]
                                       j = 0

                                                                    nx      nxb
 So: | g_p > = (x-X)^nx * |g_b> = exp(-alp*alp_b*(X-Xb)^2/gamma) * sum     sum [ C_nx^i * C_nxb^j * (Xp-X)^i * (Xp-Xb)^j * { (x-Xp)^(nx-i+nxb-j) * exp(-gamma*(x-Xp)^2) } ]
                                                                    i = 0   j=0
 =


****************************************************************************/


  double gamma = alp + alp_b;
  double Xp = (alp*X + alp_b*Xb)/gamma;
  double Xp_ = Xp - X;
  double Xpb = Xp - Xb;


  double res = 0.0;

  double pwi = 1.0;
  for(int i=0;i<=nx;i++){
    double bini = BINOM(i,nx);

    double pwj = 1.0;
    for(int j=0;j<=nxb;j++){
      double binj = BINOM(j,nxb);


      res += bini * binj * pwi * pwj * gaussian_overlap(nxa, alp_a, Xa, (nx-i+nxb-j), gamma, Xp );


      pwj *= Xpb;
    }// for j

    pwi *= Xp_;
  }// for i

  res *= exp(-alp*alp_b*(X-Xb)*(X-Xb)/gamma);

  return res;

}


double gaussian_moment(int nxa,double alp_a, double Xa, int nx, double alp, double X, int nxb,double alp_b, double Xb,
                       int is_normalize,
                       int is_derivs, double& dI_dXa, double& dI_dX,double& dI_dXb,
                       vector<double*>& aux,int n_aux ){
/****************************************************************************

****************************************************************************/


  int i;
  double gamma = alp_a + alp_b;
  double ag = alp_a/gamma;
  double bg = alp_b/gamma;

  double Xp = ag*Xa + bg*Xb;
  double Xpa = Xp - Xa;
  double Xpb = Xp - Xb;

  // Jacobian elements:
  double dXpa_dXa = ag - 1.0;
  double dXpb_dXa = ag;
  double dXpa_dXb = bg;
  double dXpb_dXb = bg - 1.0;
  double dXp_dXa = ag;
  double dXp_dXb = bg;


  // Aliaces: since places [0,1,2,3] will be used in the gaussian_overlap function, 
  // we reserve next spots:
  double* f;  f = aux[4];
  double* dfdXpa; dfdXpa = aux[5];
  double* dfdXpb; dfdXpb = aux[6];
  double* g;  g = aux[7];       // will be containing S integral
  double* dgdX;   dgdX = aux[8];
  double* dgdXp;  dgdXp = aux[9];

  // Compute binomial expansion and derivatives, if necessary - as usual
  binomial_expansion(nxa, nxb, Xpa, Xpb, f, dfdXpa, dfdXpb, is_derivs); // (x+x1)^n1 * (x+x2)^n2 = summ_i { x^i * f_i (x1,x2,n1,n2) }

  // Compute Gaussian integrals
  for(i=0;i<=(nxa+nxb+1);i++){ 

    double dS_dX, dS_dXp; 

    // In this call we do not do normalization since we are working with temporary ("virtual") Gaussians
    // normalization will be applied later to Ga an Gb
    g[i] = gaussian_overlap(nx, alp, X,   i,gamma, Xp,  0, is_derivs, dS_dX, dS_dXp, aux, n_aux );

    dgdX[i] = dS_dX;
    dgdXp[i] = dS_dXp; 

  }// for i


  // Now we are ready to put everything together and compute the overall integral, and its derivatives w.r.t. original variables
  double I = 0.0;
  dI_dXa = 0.0;
  dI_dXb = 0.0;
  dI_dX = 0.0;


  for(int i=0;i<=(nxa+nxb);i++){
    I += f[i] * g[i];

    if(is_derivs){
      // Now, unlike in overlap, the dg_dX and dg_dXp will be non-zero
      // this is because g is the overal itself

      dI_dXa += ( dfdXpa[i] * dXpa_dXa + dfdXpb[i] * dXpb_dXa ) * g[i] + f[i] * dgdXp[i] * dXp_dXa;
      dI_dXb += ( dfdXpa[i] * dXpa_dXb + dfdXpb[i] * dXpb_dXb ) * g[i] + f[i] * dgdXp[i] * dXp_dXb;
      dI_dX += f[i] * dgdX[i];

    }// is_derivs
  }

  double pref = exp(-alp_a*alp_b*(Xa-Xb)*(Xa-Xb)/gamma);
  double res = I * pref;

  if(is_derivs){  
    dI_dXa *= pref;
    dI_dXb *= pref;
    dI_dX *= pref;

    dI_dXa += (-2.0*alp_a*alp_b*(Xa-Xb)/gamma)*res;
    dI_dXb -= (-2.0*alp_a*alp_b*(Xa-Xb)/gamma)*res;
  }


  // In case we need to normalize initial Gaussians
  if(is_normalize){

    double nrm = gaussian_normalization_factor(nxa,alp_a) * gaussian_normalization_factor(nxb,alp_b);
    res *= nrm;
    dI_dXa *= nrm; 
    dI_dXb *= nrm; 
    dI_dX *= nrm;

  }

  return res;

}// gaussian_moment - very general version



double gaussian_moment(int nxa,double alp_a, double Xa, int nx, double alp, double X, int nxb,double alp_b, double Xb,
                       int is_normalize,
                       int is_derivs, double& dI_dXa, double& dI_dX,double& dI_dXb
                      ){

  // Allocate working memory
  int i;
  int n_aux = 20; //nxa+nxb+1;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_moment(nxa,alp_a,Xa,  nx,alp,X,  nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dX, dI_dXb, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;


}// gaussian_moment - without external memory allocation


boost::python::list gaussian_moment(int nxa,double alp_a, double Xa, int nx, double alp, double X, int nxb,double alp_b, double Xb,
                                    int is_normalize, int is_derivs ){
  double dI_dXa, dI_dXb, dI_dX;
  double I = gaussian_moment(nxa,alp_a,Xa, nx,alp,X, nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dX, dI_dXb);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dI_dXa);
    res.append(dI_dX);
    res.append(dI_dXb);
  }

  return res;
 
}// gaussian_moment - boost::python version



double gaussian_moment(int nxa,double alp_a, double Xa, int nx, double alp, double X, int nxb,double alp_b, double Xb,
                       int is_normalize
                      ){

  double dI_dXa, dI_dX, dI_dXb;
  double res = gaussian_moment(nxa,alp_a,Xa, nx,alp,X, nxb,alp_b,Xb, is_normalize, 0, dI_dXa, dI_dX, dI_dXb);
  return res;


}// gaussian_moment - without derivatives


double gaussian_moment(int nxa,double alp_a, double Xa, int nx, double alp, double X, int nxb,double alp_b, double Xb ){

  double res = gaussian_moment(nxa,alp_a,Xa, nx,alp,X, nxb,alp_b,Xb, 1);
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

