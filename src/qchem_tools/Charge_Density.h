/*********************************************************************************
* Copyright (C) 2016-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Charge_Density.h
  \brief This file describes the functions for computing and printing charge density information

*/

#ifndef CHARGE_DENSITY_H
#define CHARGE_DENSITY_H


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


void charge_density( Electronic_Structure& el, System& syst, vector<AO>& basis_ao, Control_Parameters& prms);
void charge_density(MATRIX& C, vector<listHamiltonian_QM>& ham, System& syst, vector<vector<int> >& active_orb, Control_Parameters& prms);
void charge_density(MATRIX& C, boost::python::list ham, System& syst, boost::python::list active_orb, Control_Parameters& prms);


}// namespace libqchem_tools
}// liblibra

#endif // CHARGE_DENSITY_H
