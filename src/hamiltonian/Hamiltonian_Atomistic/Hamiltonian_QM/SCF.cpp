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

/****************************************************************************//**
 \file SCF.cpp
 \brief Implementation of the self-consistent field methods

 Details:    

  This file contains the following functions:

  double scf_oda(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                 vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                 Electronic* el,Electronic* el0, Memory*)

  double scf_diis_fock(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                       vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                       Electronic* el,Electronic* el0, Memory*)

  double scf_diis_dm(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                     vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                     Electronic* el,Electronic* el0, Memory*)




****************************************************************************/
#include "SCF.h"

namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_qm{



double scf(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM
){

  double res = 0.0;

  if(prms.scf_algo=="oda"){
    if(prms.use_disk){
      res = scf_oda_disk(el,syst,basis_ao, prms,modprms, atom_to_ao_map,ao_to_atom_map, BM);
    } 
    else{
      res = scf_oda(el,syst,basis_ao, prms,modprms, atom_to_ao_map,ao_to_atom_map, BM);
    }
  }
/*
  else if(prms.scf_algo=="diis_fock"){
    res = scf_diis_fock(el,syst,basis_ao, prms,modprms, atom_to_ao_map,ao_to_atom_map);
  }
  else if(prms.scf_algo=="diis_dm"){
    res = scf_diis_dm(el,syst,basis_ao, prms,modprms, atom_to_ao_map,ao_to_atom_map);
  }
*/
  return res;

}

double scf(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM
){

  scf(&el,syst,basis_ao,  prms,modprms,  atom_to_ao_map,ao_to_atom_map, BM);
}



}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian




