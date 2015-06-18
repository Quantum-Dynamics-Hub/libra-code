#ifndef DYNAMICS_NUCLEAR_H
#define DYNAMICS_NUCLEAR_H

// External dependencies
#include "../mmath/libmmath.h"
#include "../hamiltonian/libhamiltonian.h"

// Dynamics classes
#include "nuclear/libnuclear.h"
#include "electronic/libelectronic.h"
#include "thermostat/libthermostat.h"

using namespace libhamiltonian;

namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
using namespace libthermostat;


// Dynamics_Nuclear.cpp
void propagate_nuclear(double dt,Nuclear* mol,Electronic* el,Hamiltonian* ham,int opt);
void propagate_nuclear(double dt,Nuclear* mol,Electronic* el,Hamiltonian* ham,Thermostat* therm, int opt);


}// namespace libdyn

#endif // DYNAMICS_NUCLEAR_H

