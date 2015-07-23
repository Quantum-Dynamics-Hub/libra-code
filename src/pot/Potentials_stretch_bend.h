#ifndef POTENTIALS_STRETCH_BEND_H
#define POTENTIALS_STRETCH_BEND_H

#include "../mmath/libmmath.h"
using namespace libmmath;
using namespace libmmath::liblinalg;

namespace libpot{



//------------------- Stretch-bend potentials --------------------------------
double Stretch_Bend_Harmonic(VECTOR& r1,VECTOR& r2,VECTOR& r3,         /* Inputs */
                             VECTOR& f1,VECTOR& f2,VECTOR& f3,         /* Outputs*/
                             double k_ijk,double k_kji, double theta_0,
                             double r_ij0,double r_kj0                 /* Parameters*/
                            );


}// namespace libpot

#endif //POTENTIALS_STRETCH_BEND_H
