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

#ifndef POTENTIALS_BONDS_H
#define POTENTIALS_BONDS_H

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libpot{

//------------------ Bond potentials ------------------------------

double Bond_Harmonic(VECTOR& ri,VECTOR& rj,  /*Inputs*/
                     VECTOR& fi,VECTOR& fj,  /*Outputs*/
                      double K, double r0);  /*Parameters*/
double Bond_Harmonic(VECTOR& ri,VECTOR& rj,  /*Inputs*/
                     VECTOR& fi,VECTOR& fj,  /*Outputs*/
                     MATRIX& Hess,
                      double K, double r0, int opt);  /*Parameters*/

boost::python::list Bond_Harmonic(VECTOR ri,VECTOR rj,    /*Inputs*/
                                  double K, double r0);   /*Parameters*/
boost::python::list Bond_Harmonic(VECTOR ri,VECTOR rj,    /*Inputs*/
                                  double K, double r0, int opt);  /*Parameters*/



double Bond_Quartic(VECTOR& ri,VECTOR& rj,  /*Inputs*/
                    VECTOR& fi,VECTOR& fj,  /*Outputs*/
                    double K, double r0);   /*Parameters*/
boost::python::list Bond_Quartic(VECTOR ri,VECTOR rj,    /*Inputs*/
                                 double K, double r0);   /*Parameters*/


double Bond_Morse(VECTOR& ri,VECTOR& rj,            /*Inputs*/
                  VECTOR& fi,VECTOR& fj,            /*Outputs*/
                  double D, double r0,double alp);  /*Parameters*/
boost::python::list Bond_Morse(VECTOR ri,VECTOR rj,              /*Inputs*/
                               double D, double r0,double alp);  /*Parameters*/



}// namespace libpot
}// liblibra

#endif // POTENTIALS_BONDS_H
