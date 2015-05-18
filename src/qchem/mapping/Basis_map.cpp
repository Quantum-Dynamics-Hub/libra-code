#include "Basis.h"

/****************************************************************************

  This file contains following functions:

  void map_atoms_and_orbitals(int Nnucl, const vector<AO>&  basis_ao, vector<vector<int> >& at_orbitals)
  void show_mapping(const vector<vector<int> >& at_orbitals)

****************************************************************************/


void map_atoms_and_orbitals(int Nnucl, const vector<AO>&  basis_ao, vector<vector<int> >& at_orbitals){
//  Global mapping 
//  vector<vector<int> > at_orbitals; // mapping of index of atom on the set of AO indices, e.g. at_orbitals[n][i] - is the index of i-th AO on atom n
//  vector<AO>  basis_ao;   // AO basis

  // Clean the output storage
  if(at_orbitals.size()>0){  at_orbitals.clear(); }

  for(int n=0;n<Nnucl;n++){

    vector<int> x;
    for(int i=0;i<basis_ao.size();i++){

      if(basis_ao[i].at_indx==n){  x.push_back(i);   }

    }// for i

    at_orbitals.push_back(x);

  }// for n

}

void show_mapping(const vector<vector<int> >& at_orbitals){

  for(int n=0;n<at_orbitals.size();n++){
  
    cout<<"List n= "<<n<<" has following entries: ";
    for(int i=0;i<at_orbitals[n].size();i++){  cout<<at_orbitals[n][i]<<"  ";   }
    cout<<endl;

  }// for n

}


