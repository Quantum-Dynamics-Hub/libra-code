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
  \file Hamiltonian_INDO.h
  \brief The file describes functions for INDO calculations
*/


#ifndef HAMILTONIAN_INDO_H
#define HAMILTONIAN_INDO_H


#include "Electronic_Structure.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{



// Hamiltonian_INDO.cpp
vector<int> compute_sorb_indices
( int sz, vector<AO>& basis_ao, vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
);

void compute_indo_core_parameters
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  vector<int>& sorb_indx,
  int opt, int a, int b, double& eri, double& V_AB
);



void compute_indo_core_parameters_derivs
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  vector<int>& sorb_indx,
  int opt, int a, int b, int c, VECTOR& deri, VECTOR& dV_AB
);
void compute_indo_core_parameters_derivs
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  vector<int>& sorb_indx,
  int opt, int a, int b, int c, VECTOR& deri, VECTOR& dV_AB,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
);
void compute_indo_core_parameters_derivs
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  vector<int>& sorb_indx,
  int opt, int a, int b, vector<VECTOR>& deri, vector<VECTOR>& dV_AB,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
);




void compute_all_indo_core_parameters_derivs
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, 
  vector<int>& sorb_indx, int opt
);
void compute_all_indo_core_parameters_derivs1
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, 
  vector<int>& sorb_indx, int opt
);




void indo_core_parameters
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int opt, int DF
);




void Hamiltonian_core_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF
);

void Hamiltonian_core_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF
);

void Hamiltonian_core_deriv_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF,
  int c,
  MATRIX* dHao_dx, MATRIX* dHao_dy, MATRIX* dHao_dz, 
  MATRIX* dSao_dx, MATRIX* dSao_dy, MATRIX* dSao_dz
);

void Hamiltonian_core_deriv_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF,
  int c,
  MATRIX& dHao_dx, MATRIX& dHao_dy, MATRIX& dHao_dz, 
  MATRIX& dSao_dx, MATRIX& dSao_dy, MATRIX& dSao_dz
);



void get_integrals(int i,int j,vector<AO>& basis_ao, double eri_aa, double G1, double F2, double& ii_jj,double& ij_ij);



void Hamiltonian_Fock_indo(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
                           Control_Parameters& prms,Model_Parameters& modprms,
                           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                          );

void Hamiltonian_Fock_indo(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
                           Control_Parameters& prms,Model_Parameters& modprms,
                           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                          );

void Hamiltonian_Fock_derivs_indo
( Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int c, 
  MATRIX* dHao_dx, MATRIX* dHao_dy, MATRIX* dHao_dz,
  MATRIX* dFao_alp_dx, MATRIX* dFao_alp_dy, MATRIX* dFao_alp_dz,
  MATRIX* dFao_bet_dx, MATRIX* dFao_bet_dy, MATRIX* dFao_bet_dz
);

void Hamiltonian_Fock_derivs_indo
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int c, 
  MATRIX& dHao_dx,     MATRIX& dHao_dy,     MATRIX& dHao_dz,
  MATRIX& dFao_alp_dx, MATRIX& dFao_alp_dy, MATRIX& dFao_alp_dz,
  MATRIX& dFao_bet_dx, MATRIX& dFao_bet_dy, MATRIX& dFao_bet_dz
);






/*
void indo_core_parameters(vector<int>&, vector<int>&, vector<AO>&, Nuclear&, vector<double>&, vector<double>&, Memory*, int);

void Hamiltonian_core_indo(Control_Parameters&, Model_Parameters&, Nuclear&,
                           vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&, 
                           MATRIX*, MATRIX*, Memory*, vector<double>&, vector<double>&);


void get_integrals(int, int, vector<AO>&, double, double, double, double&, double&);

void Hamiltonian_Fock_indo(Control_Parameters&, Model_Parameters&, Nuclear&,
                           vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&,
                           Electronic*, Memory*, vector<double>&, vector<double>&);
*/




}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra

#endif // HAMILTONIAN_INDO_H
