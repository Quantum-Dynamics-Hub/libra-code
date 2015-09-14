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

#ifndef HAMILTONIAN_QM_H
#define HAMILTONIAN_QM_H


#include "Hamiltonian_INDO.h"
#include "Hamiltonian_HF.h"



namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_qm{


// Hamiltonian_QM.cpp
void Hamiltonian_core(
  System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF
);

void Hamiltonian_Fock(
  Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
);

void Hamiltonian_Fock(
  Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
);




class listHamiltonian_QM{

public:

    listHamiltonian_QM(){ ;; }
    listHamiltonian_QM(std::string ctrl_filename,System& syst);

    int Nelec;
    int Norb;
    Control_Parameters prms; 
    Model_Parameters modprms;

    std::vector<AO> basis_ao;
    std::vector<vector<int> > atom_to_ao_map;
    std::vector<int> ao_to_atom_map;


    Electronic_Structure* el;


    // Methods
    void init(std::string ctrl_filename,System& syst);
    double compute_scf(System& syst);

    void get_parameters_from_file(std::string filename){ libcontrol_parameters::get_parameters_from_file(filename, prms); }
    void set_electronic_structure(Electronic_Structure& el_);





};





}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


#endif // HAMILTONIAN_QM_H
