/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

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

void hop(int ntraj, vector<int>& initstate, vector<Nuclear*>& mol, vector<Hamiltonian*>& ham, 
         vector<double> ksi, vector<MATRIX*>& g, int do_rescaling, int rep, int do_reverse);
vector<int>
hop(int ntraj, vector<int> initstate, vector<Nuclear>& mol, vector<Hamiltonian>& ham, 
    vector<double> ksi, vector<MATRIX>& g, int do_rescaling, int rep, int do_reverse);

boost::python::list
hop(int ntraj, boost::python::list initstate, boost::python::list mol, boost::python::list ham, 
    boost::python::list ksi, boost::python::list g, int do_rescaling, int rep, int do_reverse);


void rescale_velocities_adiabatic(int ntraj, vector<Nuclear*>& mol, vector<Hamiltonian*>& ham,
      vector<int>& new_st, vector<int>& old_st, int do_reverse);
void rescale_velocities_adiabatic(int ntraj, vector<Nuclear>& mol, vector<Hamiltonian>& ham,
      vector<int>& new_st, vector<int>& old_st, int do_reverse);






}// namespace libdyn

#endif // SURFACE_HOPPING_METHOD1_H
