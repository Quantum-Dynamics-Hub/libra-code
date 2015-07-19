#ifndef SURFACE_HOPPING_METHOD1_H
#define SURFACE_HOPPING_METHOD1_H

// External dependencies
#include "../mmath/libmmath.h"
#include "../hamiltonian/libhamiltonian.h"

// Dynamics classes
#include "nuclear/libnuclear.h"
#include "electronic/libelectronic.h"
#include "ensemble/libensemble.h"

using namespace libhamiltonian;

namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
using namespace libensemble;


void compute_hopping_probabilities_esh(Ensemble& ens, MATRIX* g, double dt, int use_boltz_factor,double T);

}// namespace libdyn

#endif // SURFACE_HOPPING_METHOD1_H
