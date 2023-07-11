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
  \brief The file describes the functions for the self-consistent field (SCF) algorithm for solving 
  stationary Schrodinger's equation.
  Here, the generic as well as specific version of the SCF-implementing functions are summarized
    
*/

#ifndef SCF_H
#define SCF_H

#include "Hamiltonian_QM.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{



double scf(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM);
double scf(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM);


double scf_oda(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM);
double scf_oda(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM);

double scf_oda_disk(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM);
double scf_oda_disk(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM);

double scf_none(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM);
double scf_none(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM);


/*
double scf_diis_fock(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM);
double scf_diis_dm(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM);
*/


}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra



#endif // SCF_H
