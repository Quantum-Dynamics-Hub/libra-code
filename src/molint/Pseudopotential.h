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

#ifndef PSEUDOPOTENTIAL_H
#define PSEUDOPOTENTIAL_H

#include "../math_specialfunctions/libspecialfunctions.h"
#include "../math_linalg/liblinalg.h"


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
                  );
double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB
                  );

boost::python::list pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                                int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                                int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                                int is_normalize, int is_derivs
                               );

double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                   int is_normalize
                  );


double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb
                  );



}// namespace libmolint
}// namespace liblibra


#endif // PSEUDOPOTENTIAL_H
