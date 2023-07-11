/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
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


#include "Electronic_Structure.h"


/// liblibra namespace
namespace liblibra{

namespace libatomistic{

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
}// namespace libatomistic
}// liblibra

#endif // HAMILTONIAN_HF_H
