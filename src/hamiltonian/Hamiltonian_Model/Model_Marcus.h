#ifndef MODEL_MARCUS_H
#define MODEL_MARCUS_H

#include "../../mmath/libmmath.h"
using namespace libmmath;

namespace libhamiltonian{
namespace libhamiltonian_model{


void Marcus_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params_);
boost::python::list Marcus_Ham(double x, boost::python::list params_);


}// namespace libhamiltonian_model
}// namespace libhamiltonian

#endif // MODEL_MARCUS
