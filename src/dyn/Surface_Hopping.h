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
  \file Surface_Hopping.h
  \brief The file describes the functions used in surface hopping methods
    
*/

#ifndef SURFACE_HOPPING_H
#define SURFACE_HOPPING_H


// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../hamiltonian/libhamiltonian.h"

// Dynamics classes
#include "nuclear/libnuclear.h"
#include "electronic/libelectronic.h"
#include "ensemble/libensemble.h"

/// liblibra namespace
namespace liblibra{

using namespace libhamiltonian;


/// libdyn namespace
namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
using namespace libensemble;


void compute_hopping_probabilities_fssh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_fssh(Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_fssh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T);
MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, CMATRIX& Hvib, double dt);


void compute_hopping_probabilities_gfsh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_gfsh(Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_gfsh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T);
MATRIX compute_hopping_probabilities_gfsh(CMATRIX& Coeff, CMATRIX& Hvib, double dt);


void compute_hopping_probabilities_mssh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_mssh(Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
                                        double dt, int use_boltz_factor,double T);
void compute_hopping_probabilities_mssh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T);
MATRIX compute_hopping_probabilities_mssh(CMATRIX& Coeff);



void hop(int& initstate, Nuclear* mol, Hamiltonian* ham, double ksi, MATRIX* g, int do_rescaling, int rep, int do_reverse);
int hop(int initstate, Nuclear& mol, Hamiltonian& ham, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse);
int hop(int initstate, Ensemble& ens, int i, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse);
int hop(int initstate, MATRIX& g, double ksi);


void rescale_velocities_adiabatic(Nuclear* mol, Hamiltonian* ham, int& new_st,int& old_st, int do_reverse);
int rescale_velocities_adiabatic(Nuclear& mol, Hamiltonian& ham, int old_st, int do_reverse);

void rescale_velocities_diabatic(Nuclear* mol, Hamiltonian* ham, int& new_st,int& old_st);
int rescale_velocities_diabatic(Nuclear& mol, Hamiltonian& ham, int old_st);


int ida(CMATRIX& Coeff, int old_st, int new_st, double E_old, double E_new, double T, double ksi);



MATRIX coherence_intervals(CMATRIX& Coeff, MATRIX& rates);
int dish(Electronic& el, MATRIX& t_m, const MATRIX& tau_m, const CMATRIX& Hvib,
          int use_boltz_flag, double Ekin, double T, double ksi1, double ksi2);
int dish(Electronic& el, Nuclear& mol, Hamiltonian& ham, 
          MATRIX& t_m, const MATRIX& tau_m, int use_boltz_flag, double T, double ksi1, double ksi2);


}// namespace libdyn
}// liblibra

#endif // SURFACE_HOPPING_H
