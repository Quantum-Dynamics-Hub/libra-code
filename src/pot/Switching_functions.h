#ifndef SWITCHING_FUNCTIONS_H
#define SWITCHING_FUNCTIONS_H

#include "../mmath/libmmath.h"
using namespace libmmath;


namespace libpot{

//------------------ Switching functions --------------------------
void SWITCH(VECTOR& r1,VECTOR&r2, double R_on,double R_off,double& SW,VECTOR& dSW);
boost::python::list SWITCH(VECTOR r1,VECTOR r2, double R_on,double R_off);

void DOUBLE_SWITCH(double x,double a,double eps,double& SW,double& dSW);
boost::python::list DOUBLE_SWITCH(double x,double a,double eps);



}// namespace libpot


#endif // SWITCHING_FUNCTIONS_H
