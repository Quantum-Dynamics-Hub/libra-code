#ifndef ENERGY_AND_FORCES_H
#define ENERGY_AND_FORCES_H

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


double compute_kinetic_energy(Nuclear* mol);
double compute_kinetic_energy(Nuclear& mol);
double compute_kinetic_energy(Ensemble& ens);

double compute_potential_energy(Nuclear* mol, Electronic* el, Hamiltonian* ham, int opt);
double compute_potential_energy(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt);
double compute_potential_energy(Ensemble& ens,int opt);

void compute_forces(Nuclear* mol, Electronic* el, Hamiltonian* ham, int opt);
void compute_forces(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt);
void compute_forces(Ensemble& ens,int opt);

void compute_energies(Ensemble* ens, double& Epot, double& Ekin, double& Etot,int opt);


}// namespace libdyn

#endif // ENERGY_AND_FORCES_H
