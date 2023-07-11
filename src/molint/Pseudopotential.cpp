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

#include "Pseudopotential.h"
#include "Moments.h"

/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{



double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB,
                   vector<double*>& auxd,int n_aux
                  ){
/*****************************************************************************

 <g_a | [C0 + C2*(r-R)^2]*exp(-alp*(r-R)^2) | g_b>

******************************************************************************/

    double res = 0.0;
    VECTOR didR, didA, didB;
    dIdR = 0.0;
    dIdA = 0.0;
    dIdB = 0.0;

    res = C0 * gaussian_moment(nxa,nya,nza,alp_a,Ra,  0,0,0,alp,R, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, didA,didR,didB,auxd,n_aux);
    if(is_derivs){
      dIdR += C0 * didR;      dIdA += C0 * didA;      dIdB += C0 * didB;
    }

    res += C2 * gaussian_moment(nxa,nya,nza,alp_a,Ra, 2,0,0,alp,R, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, didA,didR,didB,auxd,n_aux);
    if(is_derivs){      dIdR += C2 * didR;      dIdA += C2 * didA;      dIdB += C2 * didB;    }

    res += C2 * gaussian_moment(nxa,nya,nza,alp_a,Ra, 0,2,0,alp,R, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, didA,didR,didB,auxd,n_aux);
    if(is_derivs){      dIdR += C2 * didR;      dIdA += C2 * didA;      dIdB += C2 * didB;    }

    res += C2 * gaussian_moment(nxa,nya,nza,alp_a,Ra, 0,0,2,alp,R, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, didA,didR,didB,auxd,n_aux);
    if(is_derivs){      dIdR += C2 * didR;      dIdA += C2 * didA;      dIdB += C2 * didB;    }


    return res;

}// the most general


double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB
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
  double res = pseudopot02(C0, C2, alp, R, nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, dIdR, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<20;i++){ delete [] auxd[i]; }  
  auxd.clear();
 

  return res;

}// no external memory


boost::python::list pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                                int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                                int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                                int is_normalize, int is_derivs
                               ){

  VECTOR dIdA, dIdR, dIdB;
  double I = pseudopot02(C0, C2, alp, R, nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, is_derivs, dIdR, dIdA, dIdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dIdR);
    res.append(dIdA);
    res.append(dIdB);
  }

  return res;
 
}// version for python




double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                   int is_normalize
                  ){

  VECTOR dIdR,dIdA,dIdB;

  double res = pseudopot02(C0, C2, alp, R, nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, 0, dIdR, dIdA, dIdB);

  return res;

} 

double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb
                  ){

  VECTOR dIdR,dIdA,dIdB;

  double res = pseudopot02(C0, C2, alp, R, nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, 1);

  return res;
} 





}// namespace libmolint
}// namespace liblibra


