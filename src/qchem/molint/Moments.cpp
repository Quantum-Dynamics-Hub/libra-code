#include "Overlaps.h"
#include "Moments.h"

namespace libqchem{
namespace libmolint{



double gaussian_moment(int nx, double alp, double X,
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


double gaussian_moment(int nx, int ny,  int nz,  double alp, VECTOR& R,
                       int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                       int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                      ){

  VECTOR dIdA, dIdB;

  double res = gaussian_moment(nx,ny,nz,alp,R,  nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, 1, dIdA, dIdB);
  return res;
}


double gaussian_moment(int nx, int ny,  int nz,  double alp, VECTOR& R,
                       int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                       int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                       int is_normalize, 
                       VECTOR& dIdA, VECTOR& dIdB
                      ){
  // Allocate working memory
  int i;
  int n_aux = 40;
  vector<double*> auxd(20);
  for(i=0;i<20;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_moment(nx,ny,nz,alp,R,  nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<20;i++){ delete [] auxd[i]; }  
  auxd.clear();
 

  return res;
}



double gaussian_moment(int nx, int ny,  int nz,  double alp, VECTOR& R,
                       int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                       int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                       int is_normalize, 
                       VECTOR& dIdA, VECTOR& dIdB,
                       vector<double*>& auxd,int n_aux
                      ){
/********************************************************************************************
 This function computes moments:

 <g_a(x)|  (x-X)^nx * (x-Y)^ny * (x-Z)^nz * exp(-alp*(x-X)^2)  |g_b(x)> = <g_a | g_p>  - 3D verions of the above

********************************************************************************************/

  double Ix = gaussian_moment(nx, alp, R.x,  nxa, alp_a, Ra.x,  nxb, alp_b, Rb.x);
  double Iy = gaussian_moment(ny, alp, R.y,  nya, alp_a, Ra.y,  nyb, alp_b, Rb.y);
  double Iz = gaussian_moment(nz, alp, R.z,  nza, alp_a, Ra.z,  nzb, alp_b, Rb.z);

  dIdA = 0.0;
  dIdB = 0.0;

  if(is_normalize){
    Ix *= (gaussian_norm(nxa,alp_a) * gaussian_norm(nxb,alp_b));
    Iy *= (gaussian_norm(nya,alp_a) * gaussian_norm(nyb,alp_b));
    Iz *= (gaussian_norm(nza,alp_a) * gaussian_norm(nzb,alp_b));
  }

  double I = Ix * Iy * Iz;
  if(fabs(I)<1e-15){ I = 0.0; }


  return I;

}



}// namespace libmolint
}// namespace libqchem

