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

#ifndef MULTIPOLES_H
#define MULTIPOLES_H

#include "../math_specialfunctions/libspecialfunctions.h"
#include "../math_linalg/liblinalg.h"


/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{

///=========== 1D =====================
double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs, double& dI_dXa, double& dI_dXb,
  vector<double*>& aux,int n_aux 
);

double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs, double& dI_dXa, double& dI_dXb
);

boost::python::list transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs
);

double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize
);

double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb
);





///=========== 3D =====================
VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
);

VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
);

boost::python::list transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize,int is_derivs
);

VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize
);

VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb
);





}// namespace libmolint
}// namespace liblibra


#endif // MULTIPOLES_H
