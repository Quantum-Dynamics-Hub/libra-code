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

#ifndef MOMENTS_H
#define MOMENTS_H

#include "../math_specialfunctions/libspecialfunctions.h"
#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{


// Essentially the generalized overlaps
// 1D version
double gaussian_moment_ref(int nx, double alp, double X,
                           int nxa,double alp_a, double Xa,
                           int nxb,double alp_b, double Xb
                          );

double gaussian_moment(int nxa,double alp_a, double Xa, int nx, double alp, double X, int nxb,double alp_b, double Xb,
                       int is_normalize,
                       int is_derivs, double& dI_dXa, double& dI_dX,double& dI_dXb,
                       vector<double*>& aux,int n_aux );

double gaussian_moment(int nxa,double alp_a, double Xa, int nx, double alp, double X, int nxb,double alp_b, double Xb,
                       int is_normalize,
                       int is_derivs, double& dI_dXa, double& dI_dX,double& dI_dXb
                      );
boost::python::list gaussian_moment(int nxa,double alp_a, double Xa, int nx, double alp, double X, int nxb,double alp_b, double Xb,
                                    int is_normalize, int is_derivs );
double gaussian_moment(int nxa,double alp_a, double Xa, int nx, double alp, double X, int nxb,double alp_b, double Xb,
                       int is_normalize
                      );
double gaussian_moment(int nxa,double alp_a, double Xa, int nx, double alp, double X, int nxb,double alp_b, double Xb );




// 3D version
double gaussian_moment(int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                       int nx, int ny,  int nz,  double alp, const VECTOR& R,
                       int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                       int is_normalize, 
                       int is_derivs, VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB,
                       vector<double*>& auxd,int n_aux
                      );
double gaussian_moment(int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                       int nx, int ny,  int nz,  double alp, const VECTOR& R,
                       int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                       int is_normalize, 
                       int is_derivs, VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB
                       );
boost::python::list gaussian_moment(int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                                    int nx, int ny,  int nz,  double alp, const VECTOR& R,
                                    int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                                    int is_normalize, int is_derivs
                                   );
double gaussian_moment(int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                       int nx, int ny,  int nz,  double alp, const VECTOR& R,
                       int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                       int is_normalize
                      );
double gaussian_moment(int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                       int nx, int ny,  int nz,  double alp, const VECTOR& R,
                       int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb
                      );



}// namespace libmolint
}// namespace liblibra


#endif // MOMENTS_H
