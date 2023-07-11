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
  \file Hamiltonian_HF.cpp
  \brief The file implements functions for Hartree-Fock (HF) calculations
*/

#include "Hamiltonian_HF.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{


void Hamiltonian_core_hf
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[out] Hao The pointer to the matrix object in which the core Hamiltonian will be stored
  \param Sao The pointer to the AO overla matrix (not actually used here)
  \param[in] DF Debug flag - controlls how much of extra info is printed out
  
  Compute the core HF (Hartree-Fock) Hamiltonian
*/

  int i,j,n, I,J;
  VECTOR da,db,dc;
  
  int Norb = basis_ao.size(); // how many AOs

  if(Norb!=Hao->n_cols){  
    cout<<"In Hamiltonian_core_hf: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }


  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){

      Hao->M[i*Norb+j] = kinetic_integral(basis_ao[i],basis_ao[j]); //,da,db);

      for(n=0;n<syst.Number_of_atoms;n++){
        Hao->M[i*Norb+j] -= modprms.PT[syst.Atoms[n].Atom_element].Zeff 
                          * nuclear_attraction_integral(basis_ao[i],basis_ao[j], syst.Atoms[n].Atom_RB.rb_cm );// ,n,da,db,dc);
      }// for n

    }// for j
  }// for i


}

void Hamiltonian_core_hf
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[out] Hao The pointer to the matrix object in which the core Hamiltonian will be stored
  \param Sao The pointer to the AO overla matrix (not actually used here)
  \param[in] DF Debug flag - controlls how much of extra info is printed out
  
  Compute the core HF (Hartree-Fock) Hamiltonian - Python-friendly version
*/


  Hamiltonian_core_hf( syst, basis_ao, prms, modprms,  atom_to_ao_map, ao_to_atom_map, &Hao, &Sao, DF);

}


void Hamiltonian_Fock_hf(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
                         Control_Parameters& prms,Model_Parameters& modprms,
                         vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map){
/**
  \param[in,out] el The electronic structure of the system (some of the results will be printed into it)
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  
  Compute the Fock matrix of the HF (Hartree-Fock) Hamiltonian
*/

  int a,b,c,d,A,B,C,D;

  int Norb = basis_ao.size(); // how many AOs 
  if(Norb!=el->Hao->n_cols){  
    cout<<"In Hamiltonian_Fock_hf: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }

  // Update total density matrix
  *el->P = *el->P_alp + *el->P_bet;


  update_Mull_orb_pop(el->P, el->Sao, el->Mull_orb_pop_gross, el->Mull_orb_pop_net);


  vector<double> Zeff(syst.Number_of_atoms, 0.0);
  vector<double> Mull_charges_gross(syst.Number_of_atoms, 0.0);
  vector<double> Mull_charges_net(syst.Number_of_atoms, 0.0);

  for(a=0;a<syst.Number_of_atoms;a++){ Zeff[a] = modprms.PT[syst.Atoms[a].Atom_element].Zeff; } // e.g. 4 for STO-3G C

  update_Mull_charges(ao_to_atom_map, Zeff, el->Mull_orb_pop_gross, el->Mull_orb_pop_net, Mull_charges_gross, Mull_charges_net);

  for(a=0;a<syst.Number_of_atoms;a++){ 
    syst.Atoms[a].Atom_mull_charge_gross = Mull_charges_gross[a]; 
    syst.Atoms[a].Atom_mull_charge_net = Mull_charges_net[a]; 
  }



  // Compute Fock matrices: Core part
  *el->Fao_alp = *el->Hao;
  *el->Fao_bet = *el->Hao;



  // Formation of the Fock matrix: add Coulomb and Exchange parts
  for(a=0;a<Norb;a++){
    for(b=0;b<Norb;b++){
      for(c=0;c<Norb;c++){
        for(d=0;d<Norb;d++){

          //  (P_cd * (ab|cd) - P_alp_cd*(ad|cb))
          double J_abcd,K_adcb;
          modprms.hf_int.get_JK_values(a,b,c,d,J_abcd,K_adcb);

          if(prms.use_rosh){
            el->Fao_alp->M[a*Norb+b] += (el->P->M[c*Norb+d]*J_abcd - 0.5*el->P->M[c*Norb+d]*K_adcb);
            el->Fao_bet->M[a*Norb+b] += (el->P->M[c*Norb+d]*J_abcd - 0.5*el->P->M[c*Norb+d]*K_adcb);
          }
          else{
            el->Fao_alp->M[a*Norb+b] += (el->P->M[c*Norb+d]*J_abcd - el->P_alp->M[c*Norb+d]*K_adcb);
            el->Fao_bet->M[a*Norb+b] += (el->P->M[c*Norb+d]*J_abcd - el->P_bet->M[c*Norb+d]*K_adcb);
          }


        }// for d
      }// for c
    }// for b
  }// for a

}

void Hamiltonian_Fock_hf(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
                         Control_Parameters& prms,Model_Parameters& modprms,
                         vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                        ){
/**
  \param[in,out] el The electronic structure of the system (some of the results will be printed into it)
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  
  Compute the Fock matrix of the HF (Hartree-Fock) Hamiltonian - Python-friendly version
*/


  Hamiltonian_Fock_hf(&el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

}



}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra

