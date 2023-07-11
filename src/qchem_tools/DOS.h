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
  \file DOS.h
  \brief This file declares functions for computing Densities of States (atomically and orbitally-resolved)
*/

#ifndef DENSITY_OF_STATES_H
#define DENSITY_OF_STATES_H

#include "../qobjects/libqobjects.h"
#include "../chemobjects/libchemobjects.h"
#include "../basis_setups/libbasis_setups.h"
#include "../control_parameters/libcontrol_parameters.h"
#include "../model_parameters/libmodel_parameters.h"
#include "../atomistic/Hamiltonian_QM/Electronic_Structure.h"
#include "../atomistic/Hamiltonian_QM/Hamiltonian_QM.h"

/// liblibra namespace
namespace liblibra{

using namespace libqobjects;
using namespace libchemobjects;
using namespace libchemobjects::libchemsys;
using namespace libbasis_setups;
using namespace libcontrol_parameters;
using namespace libmodel_parameters;
using namespace libatomistic::libhamiltonian_qm;




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
}// liblibra

#endif // DENSITY_OF_STATES_H
