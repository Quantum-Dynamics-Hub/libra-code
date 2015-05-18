#include "Overlaps.h"

namespace libqchem{
namespace libmolint{


double gaussian_int(int n, double alp){
/****************************************************************************
 This function computes the elementary integral:  m = 2*n
   +inf
  int   {   x^2m * exp(-alp*x^2) dx }   = ( (2m-1)! /( a^(m/2) * 2^(m/2) ) ) * sqrt(pi/a)
   -inf

  http://en.wikipedia.org/wiki/Gaussian_integral
****************************************************************************/

  double res = 0.0;

  if(n%2==0){  // even    
    if(n==0){  res = sqrt(M_PI/alp); }
    else{ 
      res = (FACTORIAL(n-1)/pow(2.0*alp, n/2))*sqrt(M_PI/alp);
    }
  }
  else{        // odd
    res = 0.0;
  }

  return res;
}

double gaussian_norm(int n,double alp){

  return (1.0/sqrt(gaussian_int(2.0*n,2.0*alp)));
}




double gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb ){

/****************************************************************************
 This function computes the 1D overlap:

 <g_a(x)|g_b(x)>  - just the case of 1D Gaussians

 Remember also:

 exp(-alp_a*(x-Xa)) * exp(-alp_b*(x-Xb)) = exp(-alp_a*alp_b*(Xa-Xb)^2/gamma) * exp(-gamma*(x-Xp)^2)

 where Xp = (alp_a*Xa + alp_b * Xb)/gamma

 gamma = alp_a + alp_b

 So:
                                         nx
 (x-Xa)^nxa = (x - Xp + Xp - Xa)^nxa =  sum  [  C_nxa^i * (x-Xp)^i * (Xp-Xa)^(nxa-i) ]
                                       i = 0

                                        nxb
 (x-Xb)^nxb = (x - Xp + Xp - Xb)^nxa =  sum  [  C_nxb^j * (x-Xp)^j * (Xp-Xb)^(nxb-j) ]
                                       j = 0

                                                 nxa     nxb
 <g_a|g_b> = exp(-alp_a*alp_b*(Xa-Xb)^2/gamma) * sum    sum [ C_nxa^i * C_nxb^j * (Xp-Xa)^i * (Xp-Xb)^j * { (x-Xp)^(nxa-i+nxb-j) * exp(-gamma*(x-Xp)^2) } ]
                                                 i = 0   j=0

****************************************************************************/

  double gamma = alp_a + alp_b;
  double Xp = (alp_a*Xa + alp_b*Xb)/gamma;
  double Xpa = Xp - Xa;
  double Xpb = Xp - Xb;


  double res = 0.0;

  double pwi = 1.0;
  for(int i=0;i<=nxa;i++){
    double bini = BINOM(i,nxa);

    double pwj = 1.0;
    for(int j=0;j<=nxb;j++){
      double binj = BINOM(j,nxb);


      res += bini * binj * pwi * pwj * gaussian_int((nxa-i+nxb-j), gamma);


      pwj *= Xpb;
    }// for j

    pwi *= Xpa;
  }// for i

  res *= exp(-alp_a*alp_b*(Xa-Xb)*(Xa-Xb)/gamma);

  return res;

}


double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                       ){

  // This works even if Gaussians are not normalized
  VECTOR dIdA, dIdB;
  double res = gaussian_overlap(nxa,nya,nza,alp_a,Ra,nxb,nyb,nzb,alp_b,Rb,1,dIdA,dIdB);

  return res;
}

double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
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
  double res = gaussian_overlap(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<20;i++){ delete [] auxd[i]; }  
  auxd.clear();
 

  return res;

}


double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
                        VECTOR& dIdA, VECTOR& dIdB,
                        vector<double*>& auxd,int n_aux
                       ){
// Compute the overlap integral of two Gaussians GA and GB (each with arbitrary angular momentum)
// is_normalize - controls if the Gaussians should be normalized:
//   is_normalize = 0 - use Gaussians as they are
//   is_normalize = 1 - use normalize the overlap as if the Gaussians are normalized


  double Ix = gaussian_overlap(nxa, alp_a, Ra.x, nxb, alp_b, Rb.x);
  double Iy = gaussian_overlap(nya, alp_a, Ra.y, nyb, alp_b, Rb.y);
  double Iz = gaussian_overlap(nza, alp_a, Ra.z, nzb, alp_b, Rb.z);


  if(is_normalize){
    Ix *= (gaussian_norm(nxa,alp_a) * gaussian_norm(nxb,alp_b));
    Iy *= (gaussian_norm(nya,alp_a) * gaussian_norm(nyb,alp_b));
    Iz *= (gaussian_norm(nza,alp_a) * gaussian_norm(nzb,alp_b));
  }

  double I = Ix * Iy * Iz;
  if(fabs(I)<1e-15){ I = 0.0; }

  dIdA = 0.0;
  dIdB = 0.0;

  return I;

}





}// namespace libmolint
}// namespace libqchem


