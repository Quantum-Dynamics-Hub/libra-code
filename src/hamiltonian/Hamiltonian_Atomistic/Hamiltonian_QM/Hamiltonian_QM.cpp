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
#include "SCF.h"


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


listHamiltonian_QM::listHamiltonian_QM(std::string ctrl_filename,System& syst){

  init(ctrl_filename, syst);

}

void listHamiltonian_QM::init(std::string ctrl_filename,System& syst){

  //=========== STEP 1: Create control parameters (setting computation options) ================
  libcontrol_parameters::get_parameters_from_file(ctrl_filename, prms);


  //=========== STEP 2:  Create model parameters and load them from file (using control parameters options) ================
  // Initialize/read model parameters (need basis info)
  if(prms.hamiltonian=="eht" or prms.hamiltonian=="geht"){
    set_parameters_eht(prms, modprms);
  }
  else if(prms.hamiltonian=="indo"){
    set_parameters_indo(prms, modprms);
  }
  else if(prms.hamiltonian=="geht1"){
    set_parameters_geht1(prms, modprms); 
  }
  else if(prms.hamiltonian=="geht2"){
    set_parameters_geht2(prms, modprms); 
  }

  //=========== STEP 3: Set basis (STO-3G_DZ) ================
  //------- Input --------------
  vector<std::string> mol_at_types;
  vector<VECTOR> R;

  for(int i=0; i<syst.Number_of_atoms;i++){
    mol_at_types.push_back(syst.Atoms[i].Atom_element);
    R.push_back(syst.Atoms[i].Atom_RB.rb_cm);
  }

  //-------- Output -----------
  int verb = 0;
  set_basis_STO_3G_DZ(mol_at_types, R, modprms, verb, basis_ao, Nelec, Norb, atom_to_ao_map, ao_to_atom_map);


   //=========== STEP 4: Electronic structure ================
   el = new Electronic_Structure(Norb);
   el->Nelec = Nelec;
   el->Nocc_alp = Nelec/2;
   el->Nocc_bet = Nelec - el->Nocc_alp;


  //=========== STEP 5: Depending on hamiltonian to use, set internal parameters ================
  // this step runs after AO basis is set!!!

  if(prms.hamiltonian=="eht" || prms.hamiltonian=="geht" ||
     prms.hamiltonian=="geht1" || prms.hamiltonian=="geht2"){
      set_parameters_eht_mapping(modprms, basis_ao);
      set_parameters_eht_mapping1(modprms,syst.Number_of_atoms,mol_at_types);
  }

  //=========== STEP 6: Overlap matrix ================
  int x_period = 0;
  int y_period = 0;
  int z_period = 0;
  VECTOR t1, t2, t3;

  update_overlap_matrix(x_period, y_period, z_period, t1, t2, t3, basis_ao, *el->Sao); 

  //=========== STEP 7: Method-specific Parameters ================
  if(prms.hamiltonian=="indo"){
    int opt = 1; 
    el->Sao->Init_Unit_Matrix(1.0);  
    indo_core_parameters(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, opt,1);
  }
  else if(prms.hamiltonian=="cndo"||prms.hamiltonian=="cndo2"){
    int opt = 0; 
    el->Sao->Init_Unit_Matrix(1.0);  
    indo_core_parameters(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, opt,1);
  }


  //=========== STEP 8: Core Hamiltonian ================
  int debug = 1;
  Hamiltonian_core(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, *el->Hao,  *el->Sao, debug);


  //=========== STEP 9: Guess density matrix ================
  std::string eigen_method="generalized";
  vector<Timer> bench_t2(4);

  Fock_to_P(el->Norb, el->Nocc_alp, 1, el->Nocc_alp, eigen_method, prms.pop_opt,
            el->Fao_alp, el->Sao, el->C_alp, el->E_alp, el->bands_alp, el->occ_alp, el->P_alp,
            bench_t2);

  Fock_to_P(el->Norb, el->Nocc_bet, 1, el->Nocc_bet, eigen_method, prms.pop_opt,
            el->Fao_bet, el->Sao, el->C_bet, el->E_bet, el->bands_bet, el->occ_bet, el->P_bet,
            bench_t2);

  //=========== STEP 10: SCF solution ================
  scf(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, 0); 

    

}

void listHamiltonian_QM::set_electronic_structure(Electronic_Structure& el_){

  el = &el_;  // by reference

}


}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian



