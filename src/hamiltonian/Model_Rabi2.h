#ifndef MODEL_RABI2_H
#define MODEL_RABI2_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libhamiltonian{


void Rabi2_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params_);
boost::python::list Rabi2_Ham(double x, boost::python::list params_);


}//namespace libhamiltonian

#endif // MODEL_RABI2
