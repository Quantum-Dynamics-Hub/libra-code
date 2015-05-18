#include "Pseudopotential.h"
#include "Moments.h"

namespace libqchem{
namespace libmolint{


double pseudopot02(double C0, double C2, double alp, VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                  ){

  VECTOR dIdA, dIdB;

  double res = pseudopot02(C0, C2, alp, R, nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, 1, dIdA, dIdB);

  return res;

}

double pseudopot02(double C0, double C2, double alp, VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                   int is_normalize, 
                   VECTOR& dIdA, VECTOR& dIdB
                  ){
/*****************************************************************************

 <g_a | [C0 + C2*(r-R)^2]*exp(-alp*(r-R)^2) | g_b>

******************************************************************************/

  // Allocate working memory
  int i;
  int n_aux = 40;
  vector<double*> auxd(20);
  for(i=0;i<20;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = pseudopot02(C0, C2, alp, R, nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<20;i++){ delete [] auxd[i]; }  
  auxd.clear();
 

  return res;


}


double pseudopot02(double C0, double C2, double alp, VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                   int is_normalize, 
                   VECTOR& dIdA, VECTOR& dIdB,
                   vector<double*>& auxd,int n_aux
                  ){
/*****************************************************************************

 <g_a | [C0 + C2*(r-R)^2]*exp(-alp*(r-R)^2) | g_b>

******************************************************************************/

    double res = 0.0;

    res = C0 * gaussian_moment(0,0,0,alp,R, nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize,dIdA,dIdB,auxd,n_aux);

    res += C2 * gaussian_moment(2,0,0,alp,R, nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize,dIdA,dIdB,auxd,n_aux);
    res += C2 * gaussian_moment(0,2,0,alp,R, nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize,dIdA,dIdB,auxd,n_aux);
    res += C2 * gaussian_moment(0,0,2,alp,R, nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize,dIdA,dIdB,auxd,n_aux);


    return res;
}


}// namespace libmolint
}// namespace libqchem


