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
#ifndef HAMILTONIAN_EHT_H
#define HAMILTONIAN_EHT_H

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


namespace libhamiltonian{
namespace libhamiltonian_atomistic{
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
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian

#endif // HAMILTONIAN_EHT_H
