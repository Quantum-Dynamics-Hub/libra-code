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

#ifndef INTEGRAL_KINETIC_H
#define INTEGRAL_KINETIC_H

#include "../math_specialfunctions/libspecialfunctions.h"
#include "../math_linalg/liblinalg.h"


/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{


// 1D Gaussian kinetic energy
double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize, int is_derivs, double& dI_dxa,double& dI_dxb,
                        vector<double*>& aux,int n_aux );
double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize, int is_derivs, double& dI_dxa,double& dI_dxb
                       );
boost::python::list kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                                     int is_normalize, int is_derivs );
double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize );
double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb );




// 3D Gaussians kinetic energy
double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
                        int is_derivs,
                        VECTOR& dIdA, VECTOR& dIdB,
                        vector<double*>& auxd,int n_aux
                       );
double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, int is_derivs,
                        VECTOR& dIdA, VECTOR& dIdB
                       );
boost::python::list kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                     int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                     int is_normalize, int is_derivs
                                    );
double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize
                       );
double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                       );




}// namespace libmolint
}// namespace liblinalg


#endif // INTEGRAL_KINETIC_H
