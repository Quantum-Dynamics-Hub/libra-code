#include "Overlaps.h"
#include "Integral_Kinetic.h"

namespace libqchem{
namespace libmolint{


double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb ){

/****************************************************************************
 This function computes the 1D kinetic energy:

 <g_a(x)| -1/2 * d^2/dx^2 |g_b(x)>  - just the case of 1D Gaussians

*****************************************************************************/

  // Operator = (-1/2 * d^2/dx^2)
  double Ix,w; 

  w = (2.0*nxb + 1.0)*alp_b; 
  Ix = w * gaussian_overlap(nxa, alp_a, Xa, nxb, alp_b, Xb);

  w = 2.0*alp_b*alp_b;
  Ix -= w * gaussian_overlap(nxa, alp_a, Xa, nxb+2, alp_b, Xb);

  if(nxb>=2){ 
    w = 0.5*nxb*(nxb-1.0);
    Ix -= w * gaussian_overlap(nxa, alp_a, Xa, nxb-2, alp_b, Xb);
  }

  return Ix;

}// kinetic_integral


double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                       ){

  // This works even if Gaussians are not normalized
  VECTOR dIdA, dIdB;
  double res = kinetic_integral(nxa,nya,nza,alp_a,Ra,nxb,nyb,nzb,alp_b,Rb,1,dIdA,dIdB);

  return res;
}

double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
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
  double res = kinetic_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<20;i++){ delete [] auxd[i]; }  
  auxd.clear();
 

  return res;

}



double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
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

  double Ix = gaussian_overlap(nxa, alp_a, Ra.x, nxb, alp_b, Rb.x);
  double Iy = gaussian_overlap(nya, alp_a, Ra.y, nyb, alp_b, Rb.y);
  double Iz = gaussian_overlap(nza, alp_a, Ra.z, nzb, alp_b, Rb.z);

  double Kx = kinetic_integral(nxa, alp_a, Ra.x, nxb, alp_b, Rb.x);
  double Ky = kinetic_integral(nya, alp_a, Ra.y, nyb, alp_b, Rb.y);
  double Kz = kinetic_integral(nza, alp_a, Ra.z, nzb, alp_b, Rb.z);


  if(is_normalize){
    Ix *= (gaussian_norm(nxa,alp_a) * gaussian_norm(nxb,alp_b));
    Iy *= (gaussian_norm(nya,alp_a) * gaussian_norm(nyb,alp_b));
    Iz *= (gaussian_norm(nza,alp_a) * gaussian_norm(nzb,alp_b));

    Kx *= (gaussian_norm(nxa,alp_a) * gaussian_norm(nxb,alp_b));
    Ky *= (gaussian_norm(nya,alp_a) * gaussian_norm(nyb,alp_b));
    Kz *= (gaussian_norm(nza,alp_a) * gaussian_norm(nzb,alp_b));
  }

  double K = (Kx * Iy * Iz + Ix * Ky * Iz + Ix * Iy * Kz);
  if(fabs(K)<1e-15){ K = 0.0; }

  dIdA = 0.0;
  dIdB = 0.0;

  return K;

}// kinetic_integral


}// namespace libmolint
}// namespace libqchem



