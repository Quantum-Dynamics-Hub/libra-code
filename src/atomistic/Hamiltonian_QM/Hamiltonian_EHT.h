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
  \file Hamiltonian_EHT.h
  \brief The file describes functions for extended Huckel theory (EHT) calculations
*/

#ifndef HAMILTONIAN_EHT_H
#define HAMILTONIAN_EHT_H

#include "Electronic_Structure.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{



// Hamiltonian_EHT.cpp
void Hamiltonian_core_eht
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF
);

void Hamiltonian_core_eht
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF
);

void Hamiltonian_core_deriv_eht
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF,
  int c,
  MATRIX* dHao_dx, MATRIX* dHao_dy, MATRIX* dHao_dz, 
  MATRIX* dSao_dx, MATRIX* dSao_dy, MATRIX* dSao_dz
);

void Hamiltonian_core_deriv_eht
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF,
  int c,
  MATRIX& dHao_dx, MATRIX& dHao_dy, MATRIX& dHao_dz, 
  MATRIX& dSao_dx, MATRIX& dSao_dy, MATRIX& dSao_dz
);



void Hamiltonian_Fock_eht
( Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
);

void Hamiltonian_Fock_eht
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
);



/*


void Hamiltonian_Fock_eht(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
                           Control_Parameters& prms,Model_Parameters& modprms,
                           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                          );

void Hamiltonian_Fock_eht(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
                           Control_Parameters& prms,Model_Parameters& modprms,
                           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                          );

void Hamiltonian_Fock_derivs_eht
( Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int c, 
  MATRIX* dHao_dx, MATRIX* dHao_dy, MATRIX* dHao_dz,
  MATRIX* dFao_alp_dx, MATRIX* dFao_alp_dy, MATRIX* dFao_alp_dz,
  MATRIX* dFao_bet_dx, MATRIX* dFao_bet_dy, MATRIX* dFao_bet_dz
);

void Hamiltonian_Fock_derivs_eht
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int c, 
  MATRIX& dHao_dx,     MATRIX& dHao_dy,     MATRIX& dHao_dz,
  MATRIX& dFao_alp_dx, MATRIX& dFao_alp_dy, MATRIX& dFao_alp_dz,
  MATRIX& dFao_bet_dx, MATRIX& dFao_bet_dy, MATRIX& dFao_bet_dz
);

*/



}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra

#endif // HAMILTONIAN_EHT_H
