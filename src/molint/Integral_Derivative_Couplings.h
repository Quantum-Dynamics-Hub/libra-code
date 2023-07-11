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

#ifndef INTEGRAL_DERIVATIVE_COUPLINGS_H
#define INTEGRAL_DERIVATIVE_COUPLINGS_H

#include "../math_specialfunctions/libspecialfunctions.h"
#include "../math_linalg/liblinalg.h"


/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{


///==================== 1D derivative coupling integrals ==============

double derivative_coupling_integral
( int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs, double& dI_dXa,double& dI_dXb,
  vector<double*>& aux,int n_aux
);

boost::python::list derivative_coupling_integral
( int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs 
);


double derivative_coupling_integral
( int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
  int is_normalize
);

double derivative_coupling_integral
( int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb);



///========================== 3D derivative coupling integrals ========================


VECTOR derivative_coupling_integral
( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
  int is_normalize, int is_derivs,
  MATRIX3x3& dDdA, MATRIX3x3& dDdB,
  vector<double*>& auxd,int n_aux
);

VECTOR derivative_coupling_integral
( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
  int is_normalize, int is_derivs,
  MATRIX3x3& dDdA, MATRIX3x3& dDdB
);


boost::python::list derivative_coupling_integral
( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
  int is_normalize, int is_derivs
);

VECTOR derivative_coupling_integral
( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
  int is_normalize
);

VECTOR derivative_coupling_integral
( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
);



}// namespace libmolint
}// namespace liblinalg


#endif // INTEGRAL_DERIVATIVE_COUPLINGS_H

