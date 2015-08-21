#include "Mulliken.h"

namespace libcalculators{

void update_Mull_orb_pop(MATRIX* P, MATRIX* S, vector<double>& Mull_orb_pop_gross,vector<double>& Mull_orb_pop_net){
// Compute Mulliken orbital-resolved populations from the density matrix
// memory for Mull_orb_pop_*  is assumed already allocated
// gross and net populations are computed

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
// fragment - contains indexes of all nuclei, for which the Mull charges should be updated
// at_orbitals - contains indexes of AO for each of all nuclei in the system  - size is Nnucl
// len(Mull_orb_pop_*) = size of fragment basis = len(basis_fo)
// len(Mull_charges_*) = size of global nuclear array = len(Zeff)


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
// fragment - contains indexes of all nuclei, for which the Mull charges should be updated
// at_orbitals - contains indexes of AO for each of all nuclei in the system  - size is Nnucl
// len(Mull_orb_pop_*) = size of fragment basis = len(basis_fo)
// len(Mull_charges_*) = size of global nuclear array = len(Zeff)


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
