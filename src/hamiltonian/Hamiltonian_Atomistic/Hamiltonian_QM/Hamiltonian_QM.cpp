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

#include "Hamiltonian_QM.h"


namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_qm{


/****************************************************************************

****************************************************************************/

void Hamiltonian_core(
  System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF
){

  if(prms.hamiltonian=="hf"){

    Hamiltonian_core_hf(syst, basis_ao,  prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, DF);

  }
  else if(prms.hamiltonian=="indo"){

    Hamiltonian_core_indo(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, DF);

  }

}


void Hamiltonian_Fock(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
                      Control_Parameters& prms,Model_Parameters& modprms,
                      vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                     ){

  if(prms.hamiltonian=="hf"){

    Hamiltonian_Fock_hf(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

  }
  else if(prms.hamiltonian=="indo"){

    Hamiltonian_Fock_indo(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

  }


}

void Hamiltonian_Fock(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
                      Control_Parameters& prms,Model_Parameters& modprms,
                      vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                     ){
  Hamiltonian_Fock(&el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

}


}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian



