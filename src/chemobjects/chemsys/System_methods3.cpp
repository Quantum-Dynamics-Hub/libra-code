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

#include "System.h"


namespace libchemobjects{
namespace libchemsys{


void System::move_atom_by_index(VECTOR& displ,int indx){
  if(Atoms[indx].is_Atom_RB){
     Atoms[indx].Atom_RB.shift_position(displ); 
  }
}

void System::move_fragment_by_index(VECTOR& displ,int indx){
  if(Fragments[indx].is_Group_RB){
    Fragments[indx].Group_RB.shift_position(displ);
  }
}

void System::move_molecule_by_index(VECTOR& displ,int indx){
  if(Molecules[indx].is_Molecule_RB){
    Molecules[indx].Molecule_RB.shift_position(displ);
  }
}





void System::update_atoms_for_fragment(int indx){
  RigidBody& ftop = Fragments[indx].Group_RB;

  for(int i=0;i<Fragments[indx].Group_Size;i++){
    int at_indx = Fragments[indx].globAtom_Index[i];
    Atoms[at_indx].Atom_RB.rb_cm = ftop.get_center_in_global_frame(i); 
//    displ -= Atoms[at_indx].Atom_RB.rb_cm;
/* TEST!!!
    VECTOR displ = (Atoms[at_indx].Atom_RB.rb_cm - Atoms[at_indx].Atom_RB_old.rb_cm); 
    Atoms[at_indx].Atom_displ2 += displ * displ; 
    Atoms[at_indx].is_Atom_displ2 = 1;
*/
  }
}

void System::update_fragments_for_molecule(int indx){
/* AAAAA
  RigidBody& mtop = Molecules[indx].Molecule_RB;

  for(int i=0;i<Molecules[indx].Molecule_Size;i++){
    int fr_indx = Molecules[indx].globGroup_Index[i];
    Fragments[fr_indx].Group_RB.rb_cm = mtop.get_center_in_global_frame(i);
  }
*/
}

void System::update_atoms_for_molecule(int indx){
  RigidBody& mtop = Molecules[indx].Molecule_RB;

  for(int i=0;i<Molecules[indx].Molecule_Size;i++){
    int at_indx = Molecules[indx].globAtom_Index[i];
    Atoms[at_indx].Atom_RB.rb_cm = mtop.get_center_in_global_frame(i);
//    int fr_indx = Molecules[indx].globGroup_Index[i];
//    Fragments[fr_indx].Group_RB.rb_cm = mtop.get_center_in_global_frame(i);
//    update_atoms_for_fragment(fr_indx);
  }
}

void System::rotate_atoms_of_fragment(int indx,MATRIX3x3& R){
  for(int i=0;i<Fragments[indx].Group_Size;i++){
    int at_indx = Fragments[indx].globAtom_Index[i];
    Atoms[at_indx].Atom_RB.Rotate(R);
  }
}

void System::rotate_fragments_of_molecule(int indx,MATRIX3x3& R){
/*
  for(int i=0;i<Molecules[indx].Molecule_Size;i++){
    int fr_indx = Molecules[indx].globGroup_Index[i];
    Fragments[fr_indx].Group_RB.Rotate(R);
  }
*/
}

void System::rotate_atoms_of_molecule(int indx,MATRIX3x3& R){
/*
  for(int i=0;i<Molecules[indx].Molecule_Size;i++){
    int fr_indx = Molecules[indx].globGroup_Index[i];
    Fragments[fr_indx].Group_RB.Rotate(R);
    rotate_atoms_of_fragment(fr_indx,R);
  }
*/
}


void System::TRANSLATE_ATOM(double amount,VECTOR direction,int At){
/*******************************************************************
 Simplest manipulation
 Translates atom with atom id "int At" on amount of "double amount"
 in direction of "VECTOR direction"
********************************************************************/
  int v;
  VECTOR displ = amount * direction.unit();
  v = get_atom_index_by_atom_id(At);
  if(v!=-1){
    int grp_indx = Atoms[v].globGroup_Index;
    int mol_indx = Atoms[v].globMolecule_Index;
    RigidBody& atop = Atoms[v].Atom_RB;
    RigidBody& gtop = Fragments[grp_indx].Group_RB;
    RigidBody& mtop = Molecules[mol_indx].Molecule_RB;

    atop.shift_position(displ);
    gtop.shift_position((atop.rb_mass/gtop.rb_mass)*displ);
    mtop.shift_position((atop.rb_mass/mtop.rb_mass)*displ);
    
  }// if v!=-1
}

void System::TRANSLATE_FRAGMENT(double amount,VECTOR direction,int Fr){
/*******************************************************************
 Simplest manipulation
 Translates Fragment with fragment id "int Fr" on amount of "double amount"
 in direction of "VECTOR direction"
********************************************************************/
  int v;
  VECTOR displ = amount * direction.unit();
  v = get_fragment_index_by_fragment_id(Fr);
  if(v!=-1){
    int mol_indx = Fragments[v].globMolecule_Index;
    RigidBody& gtop = Fragments[v].Group_RB;
    RigidBody& mtop = Molecules[mol_indx].Molecule_RB;

    gtop.shift_position(displ);
    mtop.shift_position((gtop.rb_mass/mtop.rb_mass)*displ);
    update_atoms_for_fragment(v);

  }// if v!=-1
}

void System::TRANSLATE_MOLECULE(double amount,VECTOR direction,int Mol){
/*******************************************************************
 Simplest manipulation
 Translates Molecule with molecule id "int Mol" on amount of "double amount"
 in direction of "VECTOR direction"
********************************************************************/
  int v;
  VECTOR displ = amount * direction.unit();
  v = get_molecule_index_by_molecule_id(Mol);
  if(v!=-1){
    int mol_indx = Fragments[v].globMolecule_Index;
    RigidBody& mtop = Molecules[mol_indx].Molecule_RB;

    mtop.shift_position(displ);
    update_atoms_for_molecule(v);

  }// if v!=-1
}

void System::ROTATE_FRAGMENT(double degree_amount, VECTOR direction,int Gr){

  MATRIX3x3 R;
  double phi = M_PI*degree_amount/180.0; // Convert to the radians
  double cs = cos(0.5*phi);
  double si = sin(0.5*phi);
  VECTOR u = direction.unit();
  QUATERNION quat(cs,si*u.x,si*u.y,si*u.z);
  QUATERNION_TO_MATRIX(quat,R);

  int v = get_fragment_index_by_fragment_id(Gr);
  if(v!=-1){
    Fragments[v].Group_RB.Rotate(R);
    rotate_atoms_of_fragment(v,R);
  }
  // Molecule orientation does not change because the center of mass of the
  // fragment v does not change

}

void System::ROTATE_MOLECULE(double degree_amount, VECTOR direction,int Mol){

  MATRIX3x3 R;
  double phi = M_PI*degree_amount/180.0; // Convert to the radians
  double cs = cos(0.5*phi);
  double si = sin(0.5*phi);
  VECTOR u = direction.unit();
  QUATERNION quat(cs,si*u.x,si*u.y,si*u.z);
  QUATERNION_TO_MATRIX(quat,R);

  int v = get_molecule_index_by_molecule_id(Mol);
  if(v!=-1){
    Molecules[v].Molecule_RB.Rotate(R);
    rotate_atoms_of_molecule(v,R);
  }

}

}// namespace libchemsys
}// namespace libchemobjects

