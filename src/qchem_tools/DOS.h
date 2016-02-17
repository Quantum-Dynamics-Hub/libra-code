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
/**
  \file DOS.h
  \brief This file declares functions for computing Densities of States (atomically and orbitally-resolved)
*/

#ifndef DENSITY_OF_STATES_H
#define DENSITY_OF_STATES_H

#include "../qchem/libqchem.h"
using namespace libqchem;
using namespace libqchem::libqobjects;

#include "../chemobjects/libchemobjects.h"
using namespace libchemobjects;
using namespace libchemobjects::libchemsys;

#include "../hamiltonian/Hamiltonian_Atomistic/Hamiltonian_QM/Basis_Setups/libbasis_setups.h"
using namespace libhamiltonian::libhamiltonian_atomistic::libhamiltonian_qm::libbasis_setups;

#include "../hamiltonian/Hamiltonian_Atomistic/Hamiltonian_QM/Control_Parameters/libcontrol_parameters.h"
using namespace libhamiltonian::libhamiltonian_atomistic::libhamiltonian_qm::libcontrol_parameters;

#include "../hamiltonian/Hamiltonian_Atomistic/Hamiltonian_QM/Model_Parameters/libmodel_parameters.h"
using namespace libhamiltonian::libhamiltonian_atomistic::libhamiltonian_qm::libmodel_parameters;

#include "../hamiltonian/Hamiltonian_Atomistic/Hamiltonian_QM/Electronic_Structure.h"
using namespace libhamiltonian::libhamiltonian_atomistic::libhamiltonian_qm;


/// libqchem_tools namespace
namespace libqchem_tools{


void compute_dos
( Electronic_Structure& el, vector<AO>& basis_ao, Control_Parameters& prms,
  vector<int>& fragment, vector< vector<int> >& atom_to_ao_map
);

void compute_dos
( Electronic_Structure& el, vector<AO>& basis_ao, Control_Parameters& prms,
  boost::python::list fragment, vector< vector<int> >& atom_to_ao_map
);

void compute_dos
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, vector< vector<int> >& atom_to_ao_map
);



}// namespace libqchem_tools


#endif // DENSITY_OF_STATES_H
