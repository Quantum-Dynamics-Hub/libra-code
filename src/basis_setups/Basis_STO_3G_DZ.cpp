/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Basis_Setups.h"


/// liblibra namespace
namespace liblibra{


namespace libbasis_setups{



void set_basis_STO_3G_DZ(vector<std::string>& at_type, vector<VECTOR>& R,  Model_Parameters& modpar, int verb,
  vector<AO>& basis_ao, int& Nelec, int& Norb, 
  vector<vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map){
/** According to Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
// This is for any atoms with S-P shells only (no D shells!!!)
// But the exponents may be varied
// Exponent is the same for S and P orbitals
// KSI[Z] - is the exponent of the S and P valence orbitals of the atom with atomic charge Z (number in periodic table)
*/

  // Initialization
  int i;
  std::map<std::string,double>::iterator iter;

  Nelec = 0;
  Norb = 0;

  // Size of atomic system setup:
  if(at_type.size()!=R.size()){  
    cout<<"Error: size of the atom names array is not equal to the size of the coordinates array\n";
    exit(0);
  }
  int Nat = at_type.size();

  // Some cleaning
  if(basis_ao.size()>0){ basis_ao.clear(); }
  if(atom_to_ao_map.size()>0){ atom_to_ao_map.clear(); }
  if(ao_to_atom_map.size()>0){ ao_to_atom_map.clear(); }


  if(verb){ cout<<"Setting STO-3G_DZ basis:\n"; }


  for(i=0;i<Nat;i++){   

    // Compute the number of valence electrons
    Nelec += modpar.PT[at_type[i]].Nval;

    for(iter = modpar.PT[at_type[i]].IP.begin(); iter != modpar.PT[at_type[i]].IP.end(); iter++){ // Iterate over all shells of the given atom 

      std::string shell = iter->first;  
      int Nzeta = modpar.PT[at_type[i]].Nzeta[shell];
      int Nquant = modpar.PT[at_type[i]].Nquant[shell];
      double IP = modpar.PT[at_type[i]].IP[shell];
      double exp1 = modpar.PT[at_type[i]].zetas[shell][0];
      double exp2 = modpar.PT[at_type[i]].zetas[shell][1];
      double coeff1 = modpar.PT[at_type[i]].coeffs[shell][0];
      double coeff2 = modpar.PT[at_type[i]].coeffs[shell][1];

                            // Here (shell.substr(1)) we turn 1s to s, 2p to p, etc.
      add_basis_ao(at_type[i], R[i], shell.substr(1), Nzeta, Nquant, IP, exp1, exp2, coeff1, coeff2, basis_ao);

    }// for iter - all shells of the atom of given elements

    if(verb){
      cout<<"Atom "<<i<<" ("<<at_type[i]<<") contributes "<<basis_ao.size() - Norb<<" orbitals\n";
      cout<<"Atom "<<i<<" ("<<at_type[i]<<") contributes "<<modpar.PT[at_type[i]].Nval<<" electrons\n";
    }

    // Set atom to AOs mapping
    vector<int> orbs_i; // orbitals of atom i
    for(int o=Norb; o<basis_ao.size(); o++){   orbs_i.push_back(o);    }
    atom_to_ao_map.push_back(orbs_i);

    // Update total number of AO created so far
    Norb = basis_ao.size();

  }// for i - all atoms in the system

 
  // Set AOs to atom mapping
  ao_to_atom_map = vector<int>(Norb);

  for(int a=0; a<Nat; a++){
    int norbs_a = atom_to_ao_map[a].size(); //  how many orbitals on the atoms a
    for(i=0; i<norbs_a; i++){
      int I = atom_to_ao_map[a][i]; // index of orbital i of atom a 
      ao_to_atom_map[I] = a;
    }// for i
  }// for a


  cout<<"Total number of orbitals added is = "<< Norb<<endl;
  cout<<"Total number of electrons is = "<< Nelec<<endl;

  cout<<"Atom to AO mapping:\n";
  show_mapping(atom_to_ao_map);

  cout<<"AO to Atom mapping:\n";
  for(i=0; i<Norb; i++){
    cout<<"orbital "<< i << " is sitting on atom "<< ao_to_atom_map[i]<<endl;
  }

  // Now print all AOs:
  if(verb){
    cout<<"Printing atomic basis: \n";
    for(i=0;i<Norb;i++){   cout<<"orbital"<<i<<endl;   basis_ao[i].show_info();  }
  }


}//int set_basis_STO_3G_DZ(...)


boost::python::list set_basis_STO_3G_DZ(vector<std::string>& at_type, vector<VECTOR>& R,  Model_Parameters& modpar, int verb){

  vector<AO> basis_ao;
  int Nelec, Norb;
  vector<vector<int> > atom_to_ao_map;
  vector<int> ao_to_atom_map;

  set_basis_STO_3G_DZ(at_type, R, modpar, verb, basis_ao, Nelec, Norb, atom_to_ao_map, ao_to_atom_map );

  boost::python::list res;
  res.append(basis_ao);
  res.append(Nelec);
  res.append(Norb);
  res.append(atom_to_ao_map);
  res.append(ao_to_atom_map);

  return res;

}




}// namespace libbasis_setups
}// namespace liblibra


