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
  \file SCF_none.cpp
  \brief The file implements the self-consistent field (SCF) algorithm for solving 
  stationary Schrodinger's equation. This version of the algorithm does not involve
  any additional convergence acceleration techniqes
    
*/

#include "SCF.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{



double scf_none(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM){
/**
  This function implements the SCF algorithm without any additional convergnece acceleration techniques.
  Just plain standard SCF

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

  int DF = 0;
  int Norb = el->Norb;
  int Nocc_alp = el->Nocc_alp;
  int Nocc_bet = el->Nocc_bet;
  std::string eigen_method="generalized";

  vector<Timer> bench_t(10); // timers for different type of operations
  vector<Timer> bench_t2(4);


  MATRIX* P_alp_old;        P_alp_old       = new MATRIX(Norb,Norb);
  MATRIX* P_bet_old;        P_bet_old       = new MATRIX(Norb,Norb);


  //===============  Initialization =======================

  Hamiltonian_Fock(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

  if(DF){
    cout<<"Fock matrix at first iteration (alp)\n";
    el->get_Fao_alp().show_matrix();
  }

  double E = (energy_elec(el->P_alp, el->Hao, el->Fao_alp) + energy_elec(el->P_bet, el->Hao, el->Fao_bet)  );

  cout<<"Initial energy = "<< E<<endl;

  double E_old = E;
  double e_err = 2.0*prms.etol;
  double d_err = 2.0*prms.den_tol;

  *P_alp_old = *el->P_alp;
  *P_bet_old = *el->P_bet;


  //===============  Now to SCF iterations =======================
  int i = 0;
  int run = 1;
  while(run){
    

    Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, prms.pop_opt, el->Fao_alp, el->Sao, el->C_alp, el->E_alp, el->bands_alp, el->occ_alp, el->P_alp, bench_t2);
    Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, prms.pop_opt, el->Fao_bet, el->Sao, el->C_bet, el->E_bet, el->bands_bet, el->occ_bet, el->P_bet, bench_t2);
    *el->P = *el->P_alp + *el->P_bet;

    Hamiltonian_Fock(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);


    d_err = fabs((*el->P_alp - *P_alp_old).max_elt()) + fabs((*el->P_bet - *P_bet_old).max_elt());

    *P_alp_old = *el->P_alp;
    *P_bet_old = *el->P_bet;

    E_old = E;
    E = (energy_elec(el->P_alp, el->Hao, el->Fao_alp) + energy_elec(el->P_bet, el->Hao, el->Fao_bet)  );

    e_err = fabs(E_old - E);

    cout<<"Iteration "<<i<<" e_err = "<<e_err<<" d_err = "<<d_err<<" E_el = "<<E<<endl;

    if(i>prms.Niter){
        run = 0;
        cout<<"Convergence is not achieved in "<<prms.Niter<<" iterations\n";
    }
    if(e_err<prms.etol && d_err<prms.den_tol){
        run = 0;
        cout<<"Success: Convergence is achieved\n";
        cout<<"Electronic energy = "<<E<<endl;
    }

    i = i + 1;
  }// while

  delete P_alp_old;
  delete P_bet_old;

  return E;
}


double scf_none(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM
){
/**
  Python-friendly version
  This function implements the SCF algorithm without any additional convergnece acceleration techniques.
  Just plain standard SCF

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


  return scf_none(&el,syst,basis_ao,  prms,modprms,  atom_to_ao_map,ao_to_atom_map, BM);
}




}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra

    

