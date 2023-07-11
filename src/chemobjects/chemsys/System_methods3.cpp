/*********************************************************************************
* Copyright (C) 2015-2021 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file System_methods3.cpp
  \brief The file implements functions for manipulating chemical objects: translations, rotations and updates
    
*/

#include "System.h"
#include "../../Units.h"

/// liblibra namespace
namespace liblibra{

/// libchemobjects namespace
namespace libchemobjects{

/// libchemsys namespace
namespace libchemsys{


void System::move_atom_by_index(VECTOR& displ,int indx){
/**
  \param[in] displ The amount and direction of the translation
  \param[in] indx The index of the atom to be translated

  Translate atom with index indx on the amount given by displ
*/

  if(Atoms[indx].is_Atom_RB){
     Atoms[indx].Atom_RB.shift_position(displ); 
  }
}

void System::move_fragment_by_index(VECTOR& displ,int indx){
/**
  \param[in] displ The amount and direction of the translation
  \param[in] indx The index of the fragment to be translated

  Translate fragment with index indx on the amount given by displ
*/

  if(Fragments[indx].is_Group_RB){
    Fragments[indx].Group_RB.shift_position(displ);
  }
}

void System::move_molecule_by_index(VECTOR& displ,int indx){
/**
  \param[in] displ The amount and direction of the translation
  \param[in] indx The index of the molecule to be translated

  Translate molecule with index indx on the amount given by displ
*/

  if(Molecules[indx].is_Molecule_RB){
    Molecules[indx].Molecule_RB.shift_position(displ);
  }
}



void System::update_atoms_for_fragment(int indx){
/**
  \param[in] indx The index of the fragment for which we update positions of atoms

  Recompute Cartesian coordinates of all atoms included in the fragment with index indx
  This is needed when only rigid-body (fragmental) variables are propagated or changed
*/

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
/// So far this function does nothing

/* AAAAA
  RigidBody& mtop = Molecules[indx].Molecule_RB;

  for(int i=0;i<Molecules[indx].Molecule_Size;i++){
    int fr_indx = Molecules[indx].globGroup_Index[i];
    Fragments[fr_indx].Group_RB.rb_cm = mtop.get_center_in_global_frame(i);
  }
*/
}

void System::update_atoms_for_molecule(int indx){
/**
  \param[in] indx The index of the molecule for which we update positions of atoms

  Recompute Cartesian coordinates of all atoms included in the molecule with index indx
  This is needed when only rigid-body (fragmental) variables of the molecule are propagated or changed
*/

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
/**
  \param[in] amount The magnitude of translation
  \param[in] direction The vector definining the direction of the translation. The magnitude of this vector does not matter.
  \param[in] At The ID (not index!) of the atom to be translated

  Simplest manipulation
  Translates atom with atom id "int At" on amount of "double amount"
  in direction of "VECTOR direction"
*/

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
/**
  \param[in] amount The magnitude of the translation
  \param[in] direction The vector definining the direction of the translation. The magnitude of this vector does not matter.
  \param[in] Fr The ID (not index!) of the fragment to be translated

  Simplest manipulation
  Translates Fragment with fragment id "int Fr" on amount of "double amount"
  in direction of "VECTOR direction"
*/

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
/**
  \param[in] amount The magnitude of the translation
  \param[in] direction The vector definining the direction of the translation. The magnitude of this vector does not matter.
  \param[in] Mol The ID (not index!) of the molecule to be translated

  Simplest manipulation
  Translates Molecule with molecule id "int Mol" on amount of "double amount"
  in direction of "VECTOR direction"
*/

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




void System::ROTATE_FRAGMENT(double degree_amount, const VECTOR& rot_direction, int fr_id, const VECTOR& center){
/**
  \param[in] degree_amount The magnitude of rotation, in degrees
  \param[in] direction The vector definining the axis of rotation in the external coordinate system
             (moving or lab frame).The magnitude of this vector does not matter.
  \param[in] fr_id The ID (not index!) of the group/fragment to be rotated
  \param[in] center The vector defining the center of the rotating coordinate system

  Simplest manipulation
  Rotates the fragment with the fragment ID "int Gr" on amount of "double amount"
  around the axis given by "VECTOR direction" around the center given by "center"
*/


  double phi = M_PI*degree_amount/180.0; // Convert to the radians
  VECTOR u = rot_direction.unit();

  int fr_indx = get_fragment_index_by_fragment_id(fr_id);
  if(fr_indx!=-1){

    Fragments[fr_indx].Group_RB.Rotate(phi, u, center);
    update_atoms_for_fragment(fr_indx);
  }

}


void System::ROTATE_FRAGMENT(double degree_amount, const VECTOR& rot_direction, int fr_id, int center_indx){
/**
  \param[in] degree_amount The magnitude of rotation, in degrees
  \param[in] rot_direction The vector definining the axis of rotation. The magnitude of this vector does not matter.
  \param[in] fr_id The ID (not index!) of the group/fragment to be rotated
  \param[in] center_indx The index of the atom around which the rotation occurs 

  Simplest manipulation
  Rotates the fragment with the fragment ID "int Gr" on amount of "double amount"
  around the axis given by "VECTOR direction" around the given center
*/

  VECTOR center = Atoms[center_indx].Atom_RB.rb_cm;

  ROTATE_FRAGMENT(degree_amount, rot_direction, fr_id, center);

}



void System::ROTATE_FRAGMENT(double degree_amount, const VECTOR& direction, int fr_id){
/**
  \param[in] degree_amount The magnitude of rotation, in degrees
  \param[in] direction The vector definining the axis of rotation. in the external coordinate system
             (moving or lab frame). The magnitude of this vector does not matter.
  \param[in] Gr The ID (not index!) of the group/fragment to be rotated

  Simplest manipulation
  Rotates the fragment with the fragment ID "int Gr" on amount of "double amount"
  around the axis given by "VECTOR direction" around the fragment's center of mass
*/

  int fr_indx = get_fragment_index_by_fragment_id(fr_id);
  if(fr_indx!=-1){

    VECTOR center(Fragments[fr_indx].Group_RB.rb_cm);
    ROTATE_FRAGMENT(degree_amount, direction, fr_id, center);
  }

}



void System::ROTATE_FRAGMENT(const MATRIX3x3& rot, int fr_id, const VECTOR& center){
/**
  \param[in] rot The rotation matrix to be applied
  \param[in] Gr The ID (not index!) of the group/fragment to be rotated
  \param[in] center The vector defining the center of the rotating coordinate system

  Simplest manipulation
  Rotates the fragment with the fragment ID "int Gr" on amount amount and direction defined by the 
  rotation matrix. Rotation is around the center given by "center"
*/

  int fr_indx = get_fragment_index_by_fragment_id(fr_id);
  if(fr_indx!=-1){
    Fragments[fr_indx].Group_RB.Rotate(rot, center);
    update_atoms_for_fragment(fr_indx);
  }

}



void System::ROTATE_FRAGMENT(const MATRIX3x3& rot_matrix, int fr_id, int center_indx){
/**
  \param[in] rot The rotation matrix to be applied
  \param[in] Gr The ID (not index!) of the group/fragment to be rotated
  \param[in] center_indx The index of the atom around which the rotation occurs 

  Simplest manipulation
  Rotates the fragment with the fragment ID "int Gr" on amount and direction defined
  by the rotation matrix "rot". Rotation is around the position of the atoms with 
  index "center_indx".
*/

  int fr_indx = get_fragment_index_by_fragment_id(fr_id);
  if(fr_indx!=-1){

    VECTOR center = Atoms[center_indx].Atom_RB.rb_cm;
    ROTATE_FRAGMENT(rot_matrix, fr_id, center);
  }


}


void System::ROTATE_FRAGMENT(const MATRIX3x3& rot, int fr_id){
/**
  \param[in] rot The rotation matrix to be applied
  \param[in] fr_id The ID (not index!) of the group/fragment to be rotated

  Simplest manipulation
  Rotates the fragment with the fragment ID "int Gr" on amount and direction defined by the 
  rotation matrix. Rotation is around the fragment's center of mass
*/

  int fr_indx = get_fragment_index_by_fragment_id(fr_id);
  if(fr_indx!=-1){

    VECTOR center = Fragments[fr_indx].Group_RB.rb_cm;
    ROTATE_FRAGMENT(rot, fr_id, center);
  }
}





void System::ROTATE_MOLECULE(double degree_amount, VECTOR direction,int Mol){
/**
  \param[in] degree_amount The magnitude of rotation, in degrees
  \param[in] direction The vector definining the axis of rotation. The magnitude of this vector does not matter.
  \param[in] Mol The ID (not index!) of the molecule to be rotated

  Simplest manipulation
  Rotates Molecule with molecule id "int Mol" on amount of "double amount"
  around the axis given by "VECTOR direction"
*/


  MATRIX3x3 R;
  double phi = M_PI*degree_amount/180.0; // Convert to the radians
  double cs = cos(0.5*phi);
  double si = sin(0.5*phi);
  VECTOR u = direction.unit();
  QUATERNION quat(cs,si*u.x,si*u.y,si*u.z);
  QUATERNION_TO_MATRIX(quat,R);

  int v = get_molecule_index_by_molecule_id(Mol);
  if(v!=-1){
    Molecules[v].Molecule_RB.Rotate_I(degree_amount, direction);

    cout<<"We also need to update the positions of the fragments' positions and orientations. \
           Since this is not yet implemented, we just exit for now. Don't use this function yet. \n";
    exit(0);

    rotate_atoms_of_molecule(v,R);
  }

}

}// namespace libchemsys
}// namespace libchemobjects
}// liblibra
