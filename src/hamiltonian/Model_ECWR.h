#ifndef MODEL_ECWR_H
#define MODEL_ECWR_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libhamiltonian{


void ECWR_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params_);
boost::python::list ECWR_Ham(double x, boost::python::list params_);


}//namespace libhamiltonian

#endif // MODEL_ECWR
