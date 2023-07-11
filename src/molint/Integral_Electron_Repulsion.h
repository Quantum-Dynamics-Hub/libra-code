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

#ifndef ERI_H
#define ERI_H

#include "../math_specialfunctions/libspecialfunctions.h"
#include "../math_linalg/liblinalg.h"


/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{

double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize, 
   int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD,
    vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
);

double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize, 
   int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD
);

boost::python::list electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize, int is_derivs
);

double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize
);

double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd
);






}// namespace libmolint
}// namespace liblinalg

#endif //ERI_H
