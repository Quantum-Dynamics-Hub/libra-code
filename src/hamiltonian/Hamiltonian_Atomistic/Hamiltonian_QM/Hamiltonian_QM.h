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
#include "Electronic_Structure.h"



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



double energy_and_forces
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
);


void derivative_couplings
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
  int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3,
  MATRIX& Dmo_a_x, MATRIX& Dmo_a_y, MATRIX& Dmo_a_z,
  MATRIX& Dmo_b_x, MATRIX& Dmo_b_y, MATRIX& Dmo_b_z
);


VECTOR force
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int Norb, int at_indx, 
  int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3
);



// Excitations.cpp
void excite(int Norb, excitation& ex, 
            int Nocc_alp, vector< pair<int,double> >& occ_alp,
            int Nocc_bet, vector< pair<int,double> >& occ_bet);




class listHamiltonian_QM{

public:

    listHamiltonian_QM(){ add_excitation(0,1,0,1); }
    listHamiltonian_QM(std::string ctrl_filename,System& syst);

    int Nelec;
    int Norb;
    Control_Parameters prms; 
    Model_Parameters modprms;

    //============= Ground state =================
    std::vector<AO> basis_ao;
    std::vector<vector<int> > atom_to_ao_map;
    std::vector<int> ao_to_atom_map;
    Electronic_Structure* el;

    //============ Excited states ================
    vector<excitation> basis_ex;  // may be the same as in prms, but may be different


    

    // Methods
    void init(std::string ctrl_filename,System& syst);
    double compute_scf(System& syst);

    void get_parameters_from_file(std::string filename){ libcontrol_parameters::get_parameters_from_file(filename, prms); }
    void set_electronic_structure(Electronic_Structure& el_);

    double energy_and_forces(System& syst);

    void add_excitation(int f_o, int f_s, int t_o, int t_s);

    //void set_excitonic_basis(boost::python::list basis_ex);
    void excite_alp(int I,int J);
    void excite_bet(int I,int J);





};





}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


#endif // HAMILTONIAN_QM_H
