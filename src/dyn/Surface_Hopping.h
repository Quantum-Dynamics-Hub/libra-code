#ifndef SURFACE_HOPPING_H
#define SURFACE_HOPPING_H

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


void compute_hopping_probabilities_fssh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_fssh(Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_fssh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T);


void compute_hopping_probabilities_gfsh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_gfsh(Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_gfsh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T);


void compute_hopping_probabilities_mssh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_mssh(Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_mssh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T);



void hop(int& initstate, Nuclear* mol, Hamiltonian* ham, double ksi, MATRIX* g, int do_rescaling, int rep, int do_reverse);
int hop(int initstate, Nuclear& mol, Hamiltonian& ham, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse);
int hop(int initstate, Ensemble& ens, int i, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse);


void rescale_velocities_adiabatic(Nuclear* mol, Hamiltonian* ham, int& new_st,int& old_st, int do_reverse);
int rescale_velocities_adiabatic(Nuclear& mol, Hamiltonian& ham, int old_st, int do_reverse);

void rescale_velocities_diabatic(Nuclear* mol, Hamiltonian* ham, int& new_st,int& old_st);
int rescale_velocities_diabatic(Nuclear& mol, Hamiltonian& ham, int old_st);



}// namespace libdyn

#endif // SURFACE_HOPPING_H
