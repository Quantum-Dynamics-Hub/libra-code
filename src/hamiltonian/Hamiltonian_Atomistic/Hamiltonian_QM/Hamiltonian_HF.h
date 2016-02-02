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
  \file Hamiltonian_HF.h
  \brief The file describes functions for Hartree-Fock (HF) calculations
*/

#ifndef HAMILTONIAN_HF_H
#define HAMILTONIAN_HF_H

#include "../../../calculators/libcalculators.h"
using namespace libcalculators;

#include "../../../qchem/libqchem.h"
using namespace libqchem;
using namespace libqchem::libqobjects;

#include "../../../chemobjects/libchemobjects.h"
using namespace libchemobjects;
using namespace libchemobjects::libchemsys;

#include "Control_Parameters/libcontrol_parameters.h"
using namespace libhamiltonian::libhamiltonian_atomistic::libhamiltonian_qm::libcontrol_parameters;

#include "Model_Parameters/libmodel_parameters.h"
using namespace libhamiltonian::libhamiltonian_atomistic::libhamiltonian_qm::libmodel_parameters;

#include "Electronic_Structure.h"


/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_atomistic namespace
namespace libhamiltonian_atomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{



// Hamiltonian_HF.cpp
void Hamiltonian_core_hf
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF
);


void Hamiltonian_core_hf
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF
);


void Hamiltonian_Fock_hf(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
                         Control_Parameters& prms,Model_Parameters& modprms,
                         vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                        );

void Hamiltonian_Fock_hf(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
                         Control_Parameters& prms,Model_Parameters& modprms,
                         vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                        );




}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian

#endif // HAMILTONIAN_HF_H
