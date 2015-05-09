#ifndef POTENTIALS_VDW_H
#define POTENTIALS_VDW_H

#include "../mmath/libmmath.h"
#include "Switching_functions.h"
using namespace libmmath;

namespace libpot{


//--------------------------- Vdw potentials -------------------------------------------

double Vdw_LJ(VECTOR& ri,VECTOR& rj,          /*Inputs*/
              VECTOR& fi,VECTOR& fj,          /*Outputs*/
              double sigma, double espilon);  /*Parameters*/

double Vdw_Buffered14_7(VECTOR& ri,VECTOR& rj,          /*Inputs*/
                        VECTOR& fi,VECTOR& fj,          /*Outputs*/
                        double sigma, double espilon);  /*Parameters*/

double Vdw_Morse(VECTOR& ri,VECTOR& rj,            /*Inputs*/
                 VECTOR& fi,VECTOR& fj,            /*Outputs*/
                 double D, double r0,double alp);  /*Parameters*/


} // namespace libpot

#endif //POTENTIALS_VDW_H
