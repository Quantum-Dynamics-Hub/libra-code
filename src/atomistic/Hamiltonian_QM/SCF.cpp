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
  \file SCF.cpp
  \brief The file implements the self-consistent field (SCF) algorithm for solving 
  stationary Schrodinger's equation - the particular selection of the method is controlled
  by the input parameters. Options: SCF_none, SCF_oda, SCF_oda_disk
    
*/

#include "SCF.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{



double scf(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM
){
/**
  This function implements the SCF with the choice of algorithms: SCF_none, SCF_oda, SCF_oda_disk
  The choice is controlled by the parameter prms

  \param[in,out] el The pointer to the object containing all the electronic structure information (MO-LCAO coefficients, 
  density matrix, Fock, etc)
  \param[in,out] syst The reference to the object containing all the nuclear information - geometry and atomic types
  \param[in] basis_ao The vector of AO objects - the AO basis for given calculations
  \param[in] prms The object that contains all the parameters controlling the simulation - all settings, flags, etc.
  \param[in,out] modprms The object that contains all the Hamiltonian parameters for given system and method choice
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] BM Benchmark flag - if =1 - do some benchmarking, if =0 - don't do it

  Returns the converged total electronic energy 
*/



  double res = 0.0;

  if(prms.scf_algo=="none"){
    res = scf_none(el,syst,basis_ao, prms,modprms, atom_to_ao_map,ao_to_atom_map, BM);
  }
  else if(prms.scf_algo=="oda"){
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
/**
  Python-friendly version
  This function implements the SCF with the choice of algorithms: SCF_none, SCF_oda, SCF_oda_disk
  The choice is controlled by the parameter prms

  \param[in,out] el The object containing all the electronic structure information (MO-LCAO coefficients, 
  density matrix, Fock, etc)
  \param[in,out] syst The reference to the object containing all the nuclear information - geometry and atomic types
  \param[in] basis_ao The vector of AO objects - the AO basis for given calculations
  \param[in] prms The object that contains all the parameters controlling the simulation - all settings, flags, etc.
  \param[in,out] modprms The object that contains all the Hamiltonian parameters for given system and method choice
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] BM Benchmark flag - if =1 - do some benchmarking, if =0 - don't do it

  Returns the converged total electronic energy 
*/


  return scf(&el,syst,basis_ao,  prms,modprms,  atom_to_ao_map,ao_to_atom_map, BM);
}



}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra



