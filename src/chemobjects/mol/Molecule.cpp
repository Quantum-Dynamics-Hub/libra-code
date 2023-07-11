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

#include "Molecule.h"

/// liblibra namespace
namespace liblibra{




namespace libchemobjects{
namespace libmol{


void Molecule::init_variables(){
  Molecule_Number_of_bonds = 0;
  Molecule_Number_of_angles = 0;
  Molecule_Number_of_dihedrals = 0;
  Molecule_Number_of_impropers = 0;

  is_Molecule_name = 0;
  is_Molecule_id   = 0;
  is_Molecule_RB   = 0;

}

void Molecule::copy_content(const Molecule& mol){

  globMolecule_Size  = mol.globMolecule_Size;
  locMolecule_Size   = mol.locMolecule_Size;
  Molecule_Size      = mol.Molecule_Size;
//  globGroup_Index    = mol.globGroup_Index;
//  locGroup_Index     = mol.locGroup_Index;
  globAtom_Index    = mol.globAtom_Index;
  locAtom_Index     = mol.locAtom_Index;
  globMolecule_Index = mol.globMolecule_Index;
  locMolecule_Index  = mol.locMolecule_Index;

  Molecule_Number_of_bonds  = mol.Molecule_Number_of_bonds;
  Molecule_Number_of_angles = mol.Molecule_Number_of_angles;
  Molecule_Number_of_dihedrals = mol.Molecule_Number_of_dihedrals;
  Molecule_Number_of_impropers = mol.Molecule_Number_of_impropers;

  if(mol.is_Molecule_name){  Molecule_name = mol.Molecule_name; is_Molecule_name = 1; }
  if(mol.is_Molecule_id)  {  Molecule_id   = mol.Molecule_id;   is_Molecule_id   = 1; }
  if(mol.is_Molecule_RB)  {  Molecule_RB   = mol.Molecule_RB;   is_Molecule_RB   = 1; }

}

Molecule::Molecule(){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
}

Molecule::Molecule(const Molecule& mol){
  /********************
    Copy constructor
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of at object which is defined
  copy_content(mol);
}

Molecule& Molecule::operator=(const Molecule& mol){
  /********************
    Assignment operator
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of at object which is defined
  copy_content(mol);
  return *this;
}

Molecule::~Molecule(){ }



void Molecule::set(object at){
  set_value(is_Molecule_name,   Molecule_name,   at,"Molecule_name");
  set_value(is_Molecule_id,     Molecule_id,     at,"Molecule_id");
}

void Molecule::show_info(){
  std::cout<<"Molecule "<<globMolecule_Index<<" properties: "<<std::endl;
//  std::cout<<"Molecule fragment indexes (global) : ";
//  for(int i=0;i<globMolecule_Size;i++){  std::cout<<globGroup_Index[i]<<"  ";  }
  std::cout<<"Molecule atom indexes (global) : ";
  for(int i=0;i<globMolecule_Size;i++){  std::cout<<globAtom_Index[i]<<"  ";  }

  std::cout<<std::endl;

  if(is_Molecule_name)   {std::cout<<"Molecule_name = "<<Molecule_name<<std::endl;   } 
  if(is_Molecule_id)     {std::cout<<"Molecule_id = "<<Molecule_id<<std::endl;       }
}


}// namespace libmol
}// namespace libchemobjects
}// liblibra
