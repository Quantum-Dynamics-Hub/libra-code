#ifndef MODEL_DAC_H
#define MODEL_DAC_H

#include "../../mmath/libmmath.h"


namespace libhamiltonian{
namespace libhamiltonian_model{

//using namespace libmmath;
//using namespace libmmath::liblinalg;

void DAC_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params_);
boost::python::list DAC_Ham(double x, boost::python::list params_);


}// namespace libhamiltonian_model
}// namespace libhamiltonian

#endif // MODEL_DAC
