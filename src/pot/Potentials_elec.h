#ifndef POTENTIALS_ELEC_H
#define POTENTIALS_ELEC_H

#include "../mmath/libmmath.h"
using namespace libmmath;
using namespace libmmath::liblinalg;

namespace libpot{



//------------------------- Electrostatic potentials ----------------------------------------

double Elec_Coulomb(VECTOR& ri,VECTOR& rj,     /*Inputs*/
                    VECTOR& fi,VECTOR& fj,     /*Outputs*/
                    double qi,double qj,
                    double eps,double delta);  /*Parameters*/


} //namespace libpot

#endif //POTENTIALS_ELEC_H
