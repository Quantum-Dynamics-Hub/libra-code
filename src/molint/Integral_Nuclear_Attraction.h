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

#ifndef INTEGRAL_NUCLEAR_ATTRACTION_H
#define INTEGRAL_NUCLEAR_ATTRACTION_H

#include "../math_specialfunctions/libspecialfunctions.h"
#include "../math_linalg/liblinalg.h"


/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{


double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                   VECTOR& Rc,int is_normalize, 
                                   int is_derivs, VECTOR& DA,VECTOR& DB, VECTOR& DC,
                                   vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv                                   
                                  );

double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                   VECTOR& Rc,int is_normalize, 
                                   int is_derivs, VECTOR& DA, VECTOR& DB, VECTOR& DC
                                  );

boost::python::list nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                     int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                     VECTOR& Rc, int is_normalize, int is_derivs
                                    );

double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        VECTOR& Rc, int is_normalize
                       );

double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb, VECTOR& Rc
                       );


}// namespace libmolint
}// namespace liblinalg


#endif // INTEGRAL_NUCLEAR_ATTRACTION_H
