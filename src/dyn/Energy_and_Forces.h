/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Energy_and_Forces.h
  \brief The file describes the functions for dynamics-immediate energy calculations

  The "dynamics-immediate" means the energies and forces computed and organized to 
  be used in dynamical calculations of different type - classical, quantum, quantum-classical.
    
*/

#ifndef ENERGY_AND_FORCES_H
#define ENERGY_AND_FORCES_H


// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../nhamiltonian/libnhamiltonian.h"

// Dynamics classes
#include "nuclear/libnuclear.h"
#include "electronic/libelectronic.h"
//#include "ensemble/libensemble.h"
#include "dyn_control_params.h"
#include "dyn_variables.h"


/// liblibra namespace
namespace liblibra{

using namespace libnhamiltonian;

/// libdyn namespace
namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
//using namespace libensemble;


double compute_kinetic_energy(MATRIX& p, MATRIX& invM, vector<int>& which_dofs);
double compute_kinetic_energy(MATRIX& p, MATRIX& invM);

vector<double> compute_kinetic_energies(MATRIX& p, MATRIX& invM, vector<int>& which_dofs);
vector<double> compute_kinetic_energies(MATRIX& p, MATRIX& invM);



CMATRIX tsh_indx2ampl(vector<int>& res, int nstates);


double average_potential_energy(dyn_control_params& prms, dyn_variables& dyn_vars, nHamiltonian& ham);
double average_potential_energy(bp::dict prms, dyn_variables& dyn_vars, nHamiltonian& ham);

vector<double> potential_energies(dyn_control_params& prms, dyn_variables& dyn_vars, nHamiltonian& ham);
vector<double> potential_energies(bp::dict prms, dyn_variables& dyn_vars, nHamiltonian& ham);


//MATRIX aux_get_forces(dyn_control_params& prms, dyn_variables& dynvars,  nHamiltonian& ham);
//MATRIX aux_get_forces(bp::dict prms, dyn_variables& dynvars, nHamiltonian& ham);
void update_forces(dyn_control_params& prms, dyn_variables& dynvars,  nHamiltonian& ham);
void update_forces(bp::dict prms, dyn_variables& dynvars, nHamiltonian& ham);

/* old way
MATRIX aux_get_forces(dyn_control_params& prms, CMATRIX& amplitudes, vector<CMATRIX>& projectors, vector<int>& act_states, 
                      nHamiltonian& ham);
MATRIX aux_get_forces(bp::dict prms, CMATRIX& amplitudes, vector<CMATRIX>& projectors, vector<int>& act_states, 
                      nHamiltonian& ham);
*/
vector<CMATRIX> get_Eadi(nHamiltonian& ham);
vector<CMATRIX> get_Eadi(nHamiltonian* ham);

/*
double compute_kinetic_energy(Nuclear* mol);
double compute_kinetic_energy(Nuclear& mol);
double compute_kinetic_energy(Ensemble& ens);

double compute_potential_energy(Nuclear* mol, Electronic* el, Hamiltonian* ham, int opt);
double compute_potential_energy(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt);
double compute_potential_energy(Ensemble& ens,int opt);

double compute_forces(Nuclear* mol, Electronic* el, Hamiltonian* ham, int opt);
double compute_forces(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt);
double compute_forces(Ensemble& ens,int opt);

*/

}// namespace libdyn
}// liblibra

#endif // ENERGY_AND_FORCES_H
