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

/**
  \file Mulliken.cpp
  \brief The file implement functions for Mulliken population analysis    
    
*/

#include "Mulliken.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libcalculators namespace
namespace libcalculators{

void update_Mull_orb_pop(MATRIX* P, MATRIX* S, vector<double>& Mull_orb_pop_gross,vector<double>& Mull_orb_pop_net){
/**
  \brief Recompute(update) Mulliken-type orbital-resolved populations

  Compute Mulliken orbital-resolved populations from the density matrix
  memory for Mull_orb_pop_*  is assumed already allocated
  gross and net populations are computed

  \param[in] P Pointer to the density matrix 
  \param[in] S Pointer to the overlap matrix 
  \param[in,out] Mull_orb_pop_gross The vector which will collect orbital-resolved Mulliken gross populations
  \param[in,out] Mull_orb_pop_net The vector which will collect orbital-resolved Mulliken net populations
*/


  int a,b;
  double tmp_a, tmp_ab;

  int Norb = Mull_orb_pop_gross.size(); // is assumed be the same as for the other array

  for(a=0;a<Norb;a++){
    Mull_orb_pop_gross[a] = 0.0;
    Mull_orb_pop_net[a] = 0.0;
  }

  for(a=0;a<Norb;a++){

    tmp_a = 0.0;
    for(b=0;b<Norb;b++){

      tmp_ab = P->M[a*Norb+b] * S->M[a*Norb+b];
      tmp_a += tmp_ab;

      if(b==a){ Mull_orb_pop_net[a] = tmp_ab; }

    }// for b

    Mull_orb_pop_gross[a] = tmp_a;  

  }// a


}// void update_Mull_orb_pop(MATRIX* P, MATRIX* S, vector<double>& Mull_orb_pop_gross,vector<double>& Mull_orb_pop_net)


void update_Mull_charges(vector<int>& fragment, vector<int>& basis_fo, vector<vector<int> >& at_orbitals,vector<double>& Zeff,
                         vector<double>& Mull_orb_pop_gross, vector<double>& Mull_orb_pop_net,
                         vector<double>& Mull_charges_gross, vector<double>& Mull_charges_net){
/**
  \brief Update Mulliken charges on atoms

  This is older (and not very efficient) version

  \param[in] fragment  Contains indices of all nuclei, for which the Mull charges should be updated.
                        So that fragment[a] - is the global index of the a-th atom of the fragment.
  \param[in]  basis_fo  Contains the global indices of the basis functions centered on the chosen (by the "fragment" variable) fragment
  \param[in] at_orbitals Mapping bewteen the local indices of the basis functions in the atom set and the global indices of all
                        basis functions. For instance at_orbitals[n][i] is the global index of the i-th basis fucntion of n-th
                        atoms.
  \param[in] Zeff Effective charges of all nuclei
  \param[out] Mull_orb_pop_gross Mulliken gross populations on all orbitals
  \param[out] Mull_orb_pop_net Mulliken net populations on all orbitals
  \param[out] Mull_charges_gross Mulliken gross charges on all atoms 
  \param[out] Mull_charges_net Mulliken net charges on all atoms 

*/



  int Nat_frag = fragment.size(); 

  // Init arrays - only nuclear charges
  for(int a=0;a<Nat_frag;a++){  // O(Nfrag) 
    int n = fragment[a]; // index of a-th nucleus in global array
    Mull_charges_gross[n] = Zeff[n];
    Mull_charges_net[n] = Zeff[n];

    for(int i=0;i<at_orbitals[n].size();i++){  // O(Nnucl)
      int j = at_orbitals[n][i];

      for(int k=0;k<basis_fo.size();k++){      // O(Nbas)
        if(basis_fo[k]==j){  

          Mull_charges_gross[n] -= Mull_orb_pop_gross[k];
          Mull_charges_net[n] -= Mull_orb_pop_net[k];
        }// if
      }// for k
    }// for i
  }// for n

}// void update_Mull_charges(vector<double>& Mull_orb_pop_gross, vector<double>& Mull_orb_pop_net, ...


void update_Mull_charges(vector<int>& ao_to_atom_map, vector<double>& Zeff,
                         vector<double>& Mull_orb_pop_gross, vector<double>& Mull_orb_pop_net,
                         vector<double>& Mull_charges_gross, vector<double>& Mull_charges_net){
/**
  \brief Update Mulliken charges on atoms

  This is newer (and more efficient) version

  \param[in] ao_to_atom_map Mapping from the grobal indices of orbitals to the global indices nuclei:
                         ao_to_atom_map[i] - is the index of the atom on which i-th AO is localized.
  \param[in] Zeff Effective charges of all nuclei
  \param[out] Mull_orb_pop_gross Mulliken gross populations on all orbitals
  \param[out] Mull_orb_pop_net Mulliken net populations on all orbitals
  \param[out] Mull_charges_gross Mulliken gross charges on all atoms 
  \param[out] Mull_charges_net Mulliken net charges on all atoms 

*/



  int i, a;
  int Nat = Zeff.size();  // number of atoms

  for(a=0;a<Nat;a++){  // all atoms 

    Mull_charges_gross[a] = Zeff[a];
    Mull_charges_net[a] = Zeff[a];

  }// for n


  for(i=0;i<ao_to_atom_map.size();i++){  // O(AOs)

    a = ao_to_atom_map[i]; // index of atom on which i-th AO is located

    Mull_charges_gross[a] -= Mull_orb_pop_gross[i];
    Mull_charges_net[a] -= Mull_orb_pop_net[i];

  }// for i


}// void update_Mull_charges(vector<double>& Mull_orb_pop_gross, vector<double>& Mull_orb_pop_net, ...




}// namespace libcalculators

}// liblibra

