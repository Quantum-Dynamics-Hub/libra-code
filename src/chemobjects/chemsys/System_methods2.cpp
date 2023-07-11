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
  \file System_methods2.cpp
  \brief The file implements topology building functions as well as general creation (building) functions
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <map>
#endif

#include "System.h"

/// liblibra namespace
namespace liblibra{


using namespace std;


/// libchemobjects namespace
namespace libchemobjects{

/// libchemsys namespace
namespace libchemsys{


int is_in_vector(int indx,vector<int> vect){
/**
   This functions is not to be exposed to user.
   It searches for index "indx" of int type in a vector of ints
   Return 1 if "indx" has been found, 0 - otherwise
*/

  int res = 0;
  int sz = vect.size();
  for(int i=0;i<sz;i++){
    if(vect[i]==indx) { res = 1; break;}
  }
  return res;
}

void System::create_bond(int a1,int a2,int exclude12){

  Group bond;

  bond.globGroup_Size = 2;
  bond.locGroup_Size  = 2;
  bond.Group_Size     = 2;

  bond.locAtom_Index.push_back(0);
  bond.locAtom_Index.push_back(1);

  int is_valid = 1;
  // If any two atoms are the same - it is wrong angle
  if(a1==a2){ is_valid = 0; }

  if(is_valid){
    //------------------------- Check on existence -------------------
    int is_bond_exist = 0;
    is_bond_exist = (Find_Bond(a1,a2)>-1);
    //---------------------------------------------------------------

    // If angle is not yet exist
    if(!is_bond_exist){
      //----------------- Create normal bond ----------------------
      if(bond.globAtom_Index.size()>0){ bond.globAtom_Index.clear(); }
        bond.globAtom_Index.push_back(a1);
        bond.globAtom_Index.push_back(a2);
        bond.globGroup_Index = Number_of_bonds;
        Number_of_bonds++;

        bond.locGroup_Index = Molecules[Atoms[a1].globMolecule_Index].Molecule_Number_of_angles;
        Molecules[Atoms[a1].globMolecule_Index].Molecule_Number_of_bonds++;
        bond.globMolecule_Index = Atoms[a1].globMolecule_Index;

        Bonds.push_back(bond);

      //-------------- Create fragmental bond ------------------------
      int g1,g2;
      g1 = Atoms[a1].globGroup_Index;
      g2 = Atoms[a2].globGroup_Index;

      if(g1==g2) { }  // All atoms belong to the same group
      else{   // At least 1 atom belong to different group
        if(!is_in_vector(bond.globGroup_Index, Frag_bonds)){
          Frag_bonds.push_back(bond.globGroup_Index);
          Number_of_frag_bonds++;
        }
      }

/*
      //-------------- Erase 1,2-pairs if necessary ---------------------
      if(exclude12){
        int frag_pair_indx = Find_Frag_Pair(a1,a2);
        if(frag_pair_indx!=-1){
          //------- Effective erasing of the frag_pair_indx-th frag pair ------------
          int last = Number_of_frag_pairs-1;
          // Copy last pair to i-th position
          Frag_pairs[frag_pair_indx] = Frag_pairs[last];

          // However, update some information of this pair
          Pairs[frag_pair_indx].globGroup_Index = frag_pair_indx;
          Pairs[frag_pair_indx].locGroup_Index  = frag_pair_indx; // = globGroup_Index

          // Delete last pair (just copied - avoid pair repetition)
          Frag_pairs.pop_back();
          Number_of_frag_pairs--;
        }// if frag_pair_indx!=-1
      }//if exclude12
*/
    }// !is_bond_exist
  }// is_valid
}

void System::create_angle(int a1,int a2,int a3,int exclude13){

  Group angle;

  angle.globGroup_Size = 3;
  angle.locGroup_Size  = 3;
  angle.Group_Size     = 3;

  angle.locAtom_Index.push_back(0);
  angle.locAtom_Index.push_back(1);
  angle.locAtom_Index.push_back(2);

  int is_valid = 1;
  // If any two atoms are the same - it is wrong angle
  if((a1==a3)||(a1==a2)||(a2==a3)){ is_valid = 0; }

  if(is_valid){
    //------------------------- Check on existence -------------------
    int is_angle_exist = 0;
    is_angle_exist = (Find_Angle(a1,a2,a3)>-1);
    //---------------------------------------------------------------

    // If angle is not yet exist
    if(!is_angle_exist){
      //----------------- Create normal angle ----------------------
      if(angle.globAtom_Index.size()>0){ angle.globAtom_Index.clear(); }
        angle.globAtom_Index.push_back(a1);
        angle.globAtom_Index.push_back(a2);
        angle.globAtom_Index.push_back(a3);
        angle.globGroup_Index = Number_of_angles;
        Number_of_angles++;

        angle.locGroup_Index = Molecules[Atoms[a1].globMolecule_Index].Molecule_Number_of_angles;
        Molecules[Atoms[a1].globMolecule_Index].Molecule_Number_of_angles++;
        angle.globMolecule_Index = Atoms[a1].globMolecule_Index;

        Angles.push_back(angle);
    
      //-------------- Create fragmental angle ------------------------
      int g1,g2,g3;
      g1 = Atoms[a1].globGroup_Index;
      g2 = Atoms[a2].globGroup_Index;
      g3 = Atoms[a3].globGroup_Index;

      if( (g1==g2) && (g2==g3) ) { }  // All atoms belong to the same group
      else{   // At least 1 atom belong to different group
        if(!is_in_vector(angle.globGroup_Index, Frag_angles)){
          Frag_angles.push_back(angle.globGroup_Index);
          Number_of_frag_angles++;
        }
      }

/*
      //-------------- Erase 1,3-pairs is necessary ---------------------
      if(exclude13){
        int frag_pair_indx = Find_Frag_Pair(a1,a3);
        if(frag_pair_indx!=-1){
          //------- Effective erasing of the frag_pair_indx-th frag pair ------------
          int last = Number_of_frag_pairs-1;
          // Copy last pair to i-th position
          Frag_pairs[frag_pair_indx] = Frag_pairs[last];

          // However, update some information of this pair
          Pairs[frag_pair_indx].globGroup_Index = frag_pair_indx;
          Pairs[frag_pair_indx].locGroup_Index  = frag_pair_indx; // = globGroup_Index

          // Delete last pair (just copied - avoid pair repetition)
          Frag_pairs.pop_back();
          Number_of_frag_pairs--;
        }// if frag_pair_indx!=-1
      }//if exclude13
*/
    }// !is_angle_exist
  }// is_valid
}

void System::create_dihedral(int a1,int a2,int a3,int a4,int exclude14){

  Group dihedral;

  dihedral.globGroup_Size = 4;
  dihedral.locGroup_Size  = 4;
  dihedral.Group_Size     = 4;

  dihedral.locAtom_Index.push_back(0);
  dihedral.locAtom_Index.push_back(1);
  dihedral.locAtom_Index.push_back(2);
  dihedral.locAtom_Index.push_back(3);

  int is_valid = 1;
  // If any two atoms are the same - it is wrong dihedral
  if((a1==a4)||(a1==a3)||(a1==a2)||(a2==a4)||(a2==a3)||(a3==a4)){ is_valid = 0; }

  if(is_valid){
    //------------------------- Check on existence -------------------
    int is_dihedral_exist = 0;
    is_dihedral_exist = (Find_Dihedral(a1,a2,a3,a4)>-1);
    //---------------------------------------------------------------

    // If dihedral is not yet exist
    if(!is_dihedral_exist){
      //--------------- Create normal dihedral --------------------------
      if(dihedral.globAtom_Index.size()>0){ dihedral.globAtom_Index.clear(); }
        dihedral.globAtom_Index.push_back(a1);
        dihedral.globAtom_Index.push_back(a2);
        dihedral.globAtom_Index.push_back(a3);
        dihedral.globAtom_Index.push_back(a4);
        dihedral.globGroup_Index = Number_of_dihedrals;
        Number_of_dihedrals++;

        dihedral.locGroup_Index = Molecules[Atoms[a1].globMolecule_Index].Molecule_Number_of_dihedrals;
        Molecules[Atoms[a1].globMolecule_Index].Molecule_Number_of_dihedrals++;
        dihedral.globMolecule_Index = Atoms[a1].globMolecule_Index;

        Dihedrals.push_back(dihedral);

      //-------------- Create fragmental dihedral ------------------------
      int g1,g2,g3,g4;
      g1 = Atoms[a1].globGroup_Index;
      g2 = Atoms[a2].globGroup_Index;
      g3 = Atoms[a3].globGroup_Index;
      g4 = Atoms[a4].globGroup_Index;

      if( (g1==g2) && (g2==g3) && (g3==g4) ) { }  // All atoms belong to the same group
      else{   // At least 1 atom belong to different group
        if(!is_in_vector(dihedral.globGroup_Index, Frag_dihedrals)){
          Frag_dihedrals.push_back(dihedral.globGroup_Index);
          Number_of_frag_dihedrals++;
        }
      }

/*
      //-------------- Erase 1,4-pairs is necessary ---------------------
      if(exclude14){
        int frag_pair_indx = Find_Frag_Pair(a1,a4);
        if(frag_pair_indx!=-1){
          //------- Effective erasing of the frag_pair_indx-th frag pair ------------
          int last = Number_of_frag_pairs-1;
          // Copy last pair to i-th position
          Frag_pairs[frag_pair_indx] = Frag_pairs[last];

          // However, update some information of this pair
          Pairs[frag_pair_indx].globGroup_Index = frag_pair_indx;
          Pairs[frag_pair_indx].locGroup_Index  = frag_pair_indx; // = globGroup_Index

          // Delete last pair (just copied - avoid pair repetition)
          Frag_pairs.pop_back();
          Number_of_frag_pairs--;
        }// if frag_pair_indx!=-1
      }//if exclude14
*/
    }// !is_dihedral_exist
  }// is_valid
}

void System::create_improper(int a1,int a2,int a3,int a4){

  Group improper;

  improper.globGroup_Size = 4;
  improper.locGroup_Size  = 4;
  improper.Group_Size     = 4;

  improper.locAtom_Index.push_back(0);
  improper.locAtom_Index.push_back(1);
  improper.locAtom_Index.push_back(2);
  improper.locAtom_Index.push_back(3);


  // Create it if and only if atom 1 has 2 connections (so this will be third)
  if(Atoms[a1].globAtom_Adjacent_Atoms.size()==2){
    //-------------- Create normal improper -------------------------
    if(improper.globAtom_Index.size()>0){ improper.globAtom_Index.clear(); }
    improper.globAtom_Index.push_back(a1);
    improper.globAtom_Index.push_back(a2);
    improper.globAtom_Index.push_back(a3);
    improper.globAtom_Index.push_back(a4);

    improper.globGroup_Index = Number_of_impropers;
    Number_of_impropers++;

    improper.locGroup_Index = Molecules[Atoms[a1].globMolecule_Index].Molecule_Number_of_impropers;
    Molecules[Atoms[a1].globMolecule_Index].Molecule_Number_of_impropers++;
    improper.globMolecule_Index = Atoms[a1].globMolecule_Index;
    Impropers.push_back(improper); 

    //-------------- Create fragmental improper ------------------------
    int g1,g2,g3,g4;
    g1 = Atoms[a1].globGroup_Index;
    g2 = Atoms[a2].globGroup_Index;
    g3 = Atoms[a3].globGroup_Index;
    g4 = Atoms[a4].globGroup_Index;

    if( (g1==g2) && (g2==g3) && (g3==g4) ) { }  // All atoms belong to the same group
    else{   // At least 1 atom belong to different group
      if(!is_in_vector(improper.globGroup_Index, Frag_impropers)){
        Frag_impropers.push_back(improper.globGroup_Index);
        Number_of_frag_impropers++;
      }
    }
  }// if ==2
  // But if atom 1 has 3 connections already => this will be forth - so this
  // will not be an improper => erase existing improper with the central atom
  // being atom 1
  else if(Atoms[a1].globAtom_Adjacent_Atoms.size()==3){
    int impr_indx = Find_Improper(Atoms[a1].globAtom_Index);
    if(impr_indx>=0){
      // Check if this improper is in frag_impropers list
      if(is_in_vector(Impropers[impr_indx].globGroup_Index, Frag_impropers)){
        int impr_frag_indx;
        // Find an index of corresponding Frag_improper
        for(int s=0;s<Number_of_frag_impropers;s++){
          if(Impropers[impr_indx].globGroup_Index==Frag_impropers[s]) { impr_frag_indx = s; break; }
        }

        // Erase it from the Frag_impropers list:
        //------- Effective erasing of the impr_frag_indx-th (for short i-th) frag improper ------------
        int last = Number_of_frag_impropers-1;
        // Copy last pair to i-th position
        if(Frag_impropers[last]==Impropers[Number_of_impropers-1].globGroup_Index){
          Frag_impropers[impr_frag_indx] = impr_indx;
        }else{
          Frag_impropers[impr_frag_indx] = Frag_impropers[last]; //##### NEW
        }

        // Delete last pair (just copied - avoid pair repetition)
        Frag_impropers.pop_back();
        Number_of_frag_impropers--;
        //-------------------------------------------------------------------------------------------

      }// if included

      // Erase it from the Impropers list:
      //------- Effective erasing of the impr_indx-th (for short i-th) improper ------------
      int last = Number_of_impropers-1;

      // Copy last pair to i-th position
      Impropers[impr_indx] = Impropers[last];
      Impropers[impr_indx].globGroup_Index = impr_indx;
      Impropers[impr_indx].locGroup_Index = impr_indx;

      // Delete last pair (just copied - avoid pair repetition)
      Impropers.pop_back();
      Number_of_impropers--;
      //-------------------------------------------------------------------------------------------

    }// if impr_indx>=0
  }// if == 3
  else{
  // Do nothing: for connections <2 this will not create impropers
  //             for connections >3 all has been deleted already
  }

}



///==========================================================================


void System::update_max_id(){
/**
    This function updated maximal ids of Atom, Fragment and Molecule
    It is not to be exposed to user. Usually it should be placed at the end
    of other building functions.
    It is necessary because sometimes atoms, fragments or molecules are 
    deleted so we need to keep max ids updated.
    Note that max_id is not the actual number of object of corresponding type
    E.g. in a system of 2 atoms they may have ids 20 and 40, so the max_atom_id
    will be 40 not 2(of 1 as index)
*/

  int max_atom_id;
  int max_frag_id;
  int max_mol_id;
 
  // Atoms
  if(Number_of_atoms>0){
    if(Atoms[0].is_Atom_id){ max_atom_id = Atoms[0].Atom_id;}
    for(int i=1;i<Number_of_atoms;i++){
      if(Atoms[i].is_Atom_id){
        max_atom_id = (max_atom_id>=Atoms[i].Atom_id)?max_atom_id:Atoms[i].Atom_id;
      }        
    }
  }

  // Fragments
  if(Number_of_fragments>0){
    if(Fragments[0].is_Group_id){ max_frag_id = Fragments[0].Group_id;}
    for(int i=1;i<Number_of_fragments;i++){
      if(Fragments[i].is_Group_id){
        max_frag_id = (max_frag_id>=Fragments[i].Group_id)?max_frag_id:Fragments[i].Group_id;
      }
    }
  }

  // Molecules
  if(Number_of_molecules>0){
    if(Molecules[0].is_Molecule_id){ max_mol_id = Molecules[0].Molecule_id;}
    for(int i=1;i<Number_of_molecules;i++){
      if(Molecules[i].is_Molecule_id){
        max_mol_id = (max_mol_id>=Molecules[i].Molecule_id)?max_mol_id:Molecules[i].Molecule_id;
      }
    }
  }

  max_atom_id     = max_atom_id;  is_max_atom_id     = 1;
  max_fragment_id = max_frag_id;  is_max_fragment_id = 1;
  max_molecule_id = max_atom_id;  is_max_molecule_id = 1;

}

void System::CREATE_ATOM(Atom at1){
/**
  This function creates atom, fragment and molecule instances which all
  describe the same simplest object. Function adds all them to the object
  space spacified. All interrelations (cross-references) are also constructed.
*/

    Atom at = at1;

    if(is_max_atom_id){    
         at.Atom_id = max_atom_id + 1;       at.is_Atom_id = 1;         
         max_atom_id++;

    }else{
         at.Atom_id = 1;     at.is_Atom_id = 1;
         max_atom_id = 1;    is_max_atom_id = 1;
    }
      at.globAtom_Index = Number_of_atoms;
      Number_of_atoms++;


    // This atom automatically becomes a single group(fragment)
    Group fr;
    if(is_max_fragment_id){
         fr.Group_id  = max_fragment_id + 1; fr.is_Group_id = 1;
         max_fragment_id++;
    }else{
         fr.Group_id = 1;        fr.is_Group_id = 1;
         max_fragment_id = 1;    max_fragment_id = 1;
    }
      fr.globGroup_Index = Number_of_fragments;
      Number_of_fragments++;

    // This atom automatically becomes a single molecule
    Molecule mol;
    if(max_molecule_id){
         mol.Molecule_id  = max_molecule_id + 1; mol.is_Molecule_id = 1;
         max_molecule_id++;
    }else{
         mol.Molecule_id = 1;     mol.is_Molecule_id = 1;
         max_molecule_id = 1;     max_molecule_id = 1;
    }
      mol.globMolecule_Index = Number_of_molecules;
      Number_of_molecules++;
     

    // Make cross-references
    at.locAtom_Index      = 0;
    at.globGroup_Index    = fr.globGroup_Index;
    at.globMolecule_Index = mol.globMolecule_Index;

    fr.Group_Size = 1;
    fr.globGroup_Size = 1;
    fr.locGroup_Size = 1;
    fr.globAtom_Index.push_back(at.globAtom_Index);
    fr.locAtom_Index.push_back(at.locAtom_Index);
    fr.locGroup_Index = 0;
    fr.globMolecule_Index = mol.globMolecule_Index;

    mol.Molecule_Size = 1;
    mol.globMolecule_Size = 1;
    mol.locMolecule_Size = 1;
// AAAAA Due to cancellation of globGroup_Index member of molecule class
//    mol.globGroup_Index.push_back(fr.globGroup_Index);
//    mol.locGroup_Index.push_back(fr.locGroup_Index);
// AAAAA Instead:
    mol.globAtom_Index.push_back(at.globAtom_Index);
    mol.locAtom_Index.push_back(at.locAtom_Index);

    mol.locMolecule_Index = mol.globMolecule_Index;    // !!! 
    
    // Create interfragmental pairs
    // Run over all atoms that existed BEFORE (that is why there is -1) this atom
    //for(int i=0;i<(Number_of_atoms-1);i++){

    Atoms.push_back(at);

    // Run over ALL atoms that currently exist, so this will include self-self pairs!
    for(int i=0;i<Number_of_atoms;i++){

        Group pair;
 
        pair.globGroup_Size = 2;
        pair.locGroup_Size  = 2;
        pair.Group_Size     = 2;

        pair.globAtom_Index.push_back(Atoms[i].globAtom_Index);
        pair.globAtom_Index.push_back(at.globAtom_Index);

        pair.locAtom_Index.push_back(0);
        pair.locAtom_Index.push_back(1);

        pair.globGroup_Index = Number_of_frag_pairs;  // Index of the pair in array of frag pairs
        Number_of_frag_pairs++;  // Increment the total number of pairs in our object space
        Number_of_pairs++;       // Increment the total number of pairs in our object space

        //------------- Create Pair itself --------------------------
       
        // In contrast to bonds the atoms of the pair may belong to different molecules
        // Thus we can not define local index of the pair - they do not stick to some
        // particular molecules. Instead they stick to entire object space
        pair.locGroup_Index  = pair.globGroup_Index;

        // Add the pair to the object space
//        Frag_pairs.push_back(pair.globGroup_Index);   //
//        Pairs.push_back(pair);                        //
        
    }// for i



    // Add all objects to the Object Space

//    Atoms.push_back(at);

    Fragments.push_back(fr);
    Molecules.push_back(mol);

    VERTEX<Atom*> vrtx(&Atoms[at.globAtom_Index]);
    AtomGraph.ADD_VERTEX(vrtx); is_AtomGraph = 1;

    // Increment number of degrees of freedom
    Nf_t += 3;

    // Undefine translation vectors
//    is_TV1 = 0;
//    is_TV2 = 0;
//    is_TV3 = 0;
   is_Box = 0;

}

/*
void System::CREATE_ATOM(){
  Atom at; // Using the default constructor
  CREATE_ATOM(at);
}
*/

int System::is_in_vector(int indx,vector<int>& vect){
/**
  This functions is not to be exposed to user.
  It searches for index "indx" of int type in a vector of ints
  Return 1 if "indx" has been found, 0 - otherwise
*/

   int res = 0;
   int sz = vect.size();
   for(int i=0;i<sz;i++){
       if(vect[i]==indx) { res = 1; break;}
   } 
   return res;
}

void System::LINK_ATOMS(Atom& at1,Atom& at2){
/**
  \param[in] at1 One of the connected atoms
  \param[in] at2 One of the connected atoms

  This function performs one of the simplest operations on two atoms
  It connects them. Other connections also taken into account - thus
  we create not only bonds (pairs of atoms), but also angles, dihedrals, etc.
*/
  int i,adj,adj1,adj2;
  // Here are some control parameters which then should be made user-definable
  // They control if we want to exclude 1,2, 1,3 and 1,4 - pairs from the pair list
  // This is related to force field specifications
  // for example 1,4 atoms are not exluded in UFF, while 1,2 and 1,3 are excluded
  int exclude12 = 1;
  int exclude13 = 1;
  int exclude14 = 1;

  // Check if this bond already exists in the object space
  int is_exist = 0;
  int is_included;
  int n_bonds = Number_of_bonds;
  is_exist = (Find_Bond(at1.globAtom_Index,at2.globAtom_Index)>-1);

  if(!is_exist){ // Create only if this bond is not yet defined 

   if(at1.globMolecule_Index==at2.globMolecule_Index){  }
   else{ // Connect atoms of different molecules -> total number of molecules decreases 

     // Molecule with the smallest index will remain (and actually become a union of both molecules)
     // Molecule with the biggest index will be deleted. All other molecules (which had index bigger)
     // will be renamed (shifted)
    
     //##################### Step 0 ############################################
     // First - look for biggest, smallest and the last indexes
     int max  = (at1.globMolecule_Index>at2.globMolecule_Index)?at1.globMolecule_Index:at2.globMolecule_Index;
     int min  = (at1.globMolecule_Index<at2.globMolecule_Index)?at1.globMolecule_Index:at2.globMolecule_Index;
     int last_mol = Number_of_molecules-1;
 
     //###################### Step 1 ############################################
     // Merge max molecule to min one = Copy and Rename (Reindex)    
     Molecules[min].globMolecule_Size += Molecules[max].globMolecule_Size;
     Molecules[min].locMolecule_Size  += Molecules[max].locMolecule_Size;
     Molecules[min].Molecule_Size     += Molecules[max].Molecule_Size;

     Molecules[min].Molecule_Number_of_bonds     += Molecules[max].Molecule_Number_of_bonds;
     Molecules[min].Molecule_Number_of_angles    += Molecules[max].Molecule_Number_of_angles;
     Molecules[min].Molecule_Number_of_dihedrals += Molecules[max].Molecule_Number_of_dihedrals;

     // Update molecule IDs
     //Minimal maximal molecule inherits the ID of minimal one
     // Other ID decrease by 1:
     //os.ObjectSpace_Molecules[max].Molecule_id = 
     
     // Molecule
     for(i=0;i<Molecules[max].Molecule_Size;i++){
         /* AAAAA Due to cancellation of .globGroup_Index member in Molecule class
         // Copy and reindex Fragments
         int gr_indx = Molecules[max].globGroup_Index[i];        
         Molecules[min].globGroup_Index.push_back(gr_indx);
	 Molecules[min].locGroup_Index.push_back(Molecules[min].Molecule_Size + Molecules[max].locGroup_Index[i]);        
         */
         // Copy and reindex Atoms
         int at_indx = Molecules[max].globAtom_Index[i];        
         Molecules[min].globAtom_Index.push_back(at_indx);
         Molecules[min].locAtom_Index.push_back(Molecules[min].Molecule_Size + Molecules[max].locAtom_Index[i]);        


     }// for i
     
     // Fragments
     for(i=0;i<Number_of_fragments;i++){
       
        if(Fragments[i].globMolecule_Index==max){
           Fragments[i].globMolecule_Index=min;           
           Fragments[i].locGroup_Index    += Molecules[min].Molecule_Size;
        }
        else if(Fragments[i].globMolecule_Index>max){
           Fragments[i].globMolecule_Index--;
        }
     }// for i

     // Atoms
     for(i=0;i<Number_of_atoms;i++){

        if(Atoms[i].globMolecule_Index==max){
           Atoms[i].globMolecule_Index=min;                     
        }
        else if(Atoms[i].globMolecule_Index>max){
           Atoms[i].globMolecule_Index--;
        }
     }// for i

     // Bonds
     for(i=0;i<Number_of_bonds;i++){

        if(Bonds[i].globMolecule_Index==max){
           Bonds[i].globMolecule_Index=min;
           Bonds[i].locGroup_Index += Molecules[min].Molecule_Number_of_bonds;
        }
        else if(Bonds[i].globMolecule_Index>max){
           Bonds[i].globMolecule_Index--;
        }
     }// for i

     // Angles
     for(i=0;i<Number_of_angles;i++){
       if(Angles[i].globMolecule_Index==max){
  	  Angles[i].globMolecule_Index=min;
          Angles[i].locGroup_Index += Molecules[min].Molecule_Number_of_angles;
       }
       else if(Angles[i].globMolecule_Index>max){
 	  Angles[i].globMolecule_Index--;
       }
     }// for i

     // Dihedrals
     for(i=0;i<Number_of_dihedrals;i++){
       if(Dihedrals[i].globMolecule_Index==max){
  	  Dihedrals[i].globMolecule_Index=min;
          Dihedrals[i].locGroup_Index += Molecules[min].Molecule_Number_of_dihedrals;
       }
       else if(Dihedrals[i].globMolecule_Index>max){
 	  Dihedrals[i].globMolecule_Index--;
       }
     }// for i

     //****************************************************************
     // Decrement id's of those molecules for which globMolecule_Index is bigger than max
     for(i=max+1;i<Number_of_molecules;i++){
         Molecules[i].globMolecule_Index--;
         Molecules[i].locMolecule_Index = Molecules[i].globMolecule_Index;
         Molecules[i].Molecule_id--;
     }
     // Also modify the global properties of the object space
     max_molecule_id--;
     //*****************************************************************

     // Delete molecule with Molecule_Index = max and decrement number of molecules
     
     // Old version: Time-consuming 
     vector<Molecule>::iterator it;    
     it = Molecules.begin();
     Molecules.erase(it+max);
     Number_of_molecules--;
     
     // Finally update variables:
     at1.globMolecule_Index = min;
     at2.globMolecule_Index = min;

   }// else if(at1.globMolecule_Index!=at2.globMolecule_Index)

   //##################### Step 2 ########################################
   // Now we are ready to create all new entries - in fact this should be
   // a common procedure for two cases: at1.Mol_Indx==at2.Mol_Indx and
   // at1.Mol_Indx!=at2.Mol_Indx

   //----------------- Bonds ---------------------
   create_bond(at1.globAtom_Index,at2.globAtom_Index,exclude12);
   EDGE<Group*> eg(at1.globAtom_Index,at2.globAtom_Index,0);
   AtomGraph.ADD_EDGE(eg);

   //----------------- Angles --------------------
   // Angles of type at1.Atom_Index - at2.Atom_Index - atoms adjacent to at2
   for(adj=0;adj<at2.globAtom_Adjacent_Atoms.size();adj++){
     int adj_indx = at2.globAtom_Adjacent_Atoms[adj];
     create_angle(at1.globAtom_Index,at2.globAtom_Index,adj_indx,exclude13);
   }
   // Angles of type at2.Atom_Index - at1.Atom_Index - atoms adjacent to at1
   for(adj=0;adj<at1.globAtom_Adjacent_Atoms.size();adj++){
     int adj_indx = at1.globAtom_Adjacent_Atoms[adj];
     create_angle(adj_indx,at1.globAtom_Index,at2.globAtom_Index,exclude13);
   }
  
   //----------------- Dihedrals -------------------
   // Dihedrals of type adj1 - at1 - at2 - adj2
   for(adj1=0;adj1<at1.globAtom_Adjacent_Atoms.size();adj1++){
     for(adj2=0;adj2<at2.globAtom_Adjacent_Atoms.size();adj2++){
       int adj_indx1 = at1.globAtom_Adjacent_Atoms[adj1];
       int adj_indx2 = at2.globAtom_Adjacent_Atoms[adj2];
       create_dihedral(adj_indx1,at1.globAtom_Index,at2.globAtom_Index,adj_indx2,exclude14);
     }
   }
   // Dihedrals of type at1 - at2 - adj2 - adj22
   for(adj2=0;adj2<at2.globAtom_Adjacent_Atoms.size();adj2++){
     Atom& at22 = Atoms[at2.globAtom_Adjacent_Atoms[adj2]];
     for(int adj22=0;adj22<at22.globAtom_Adjacent_Atoms.size();adj22++){
       int adj_indx2 = at22.globAtom_Adjacent_Atoms[adj22];
       create_dihedral(at1.globAtom_Index,at2.globAtom_Index,at22.globAtom_Index,adj_indx2,exclude14);
     }
   }
   // Dihedrals of type adj11 - adj1 - at1 - at2  (on the contrary - this one - adds two)
   for(adj1=0;adj1<at1.globAtom_Adjacent_Atoms.size();adj1++){
     Atom& at11 = Atoms[at1.globAtom_Adjacent_Atoms[adj1]];
     for(int adj11=0;adj11<at11.globAtom_Adjacent_Atoms.size();adj11++){
       int adj_indx1 = at11.globAtom_Adjacent_Atoms[adj11];
       create_dihedral(adj_indx1,at11.globAtom_Index,at1.globAtom_Index,at2.globAtom_Index,exclude14);
     }
   } 

   //--------------- Impropers -------------------------------
   if(at1.globAtom_Adjacent_Atoms.size()>=2){ 
   create_improper(at1.globAtom_Index,at1.globAtom_Adjacent_Atoms[0],at1.globAtom_Adjacent_Atoms[1],at2.globAtom_Index);
   } 
   if(at2.globAtom_Adjacent_Atoms.size()>=2){
   create_improper(at2.globAtom_Index,at2.globAtom_Adjacent_Atoms[0],at2.globAtom_Adjacent_Atoms[1],at1.globAtom_Index);
   }

   // Add each other to indexes to corresponding variable of adjacent atom indexes(global)
   at1.globAtom_Adjacent_Atoms.push_back(at2.globAtom_Index);
   at2.globAtom_Adjacent_Atoms.push_back(at1.globAtom_Index);

   }// if not exist
 
}


void System::LINK_ATOMS(int At1_id,int At2_id){
/**
  \param[in] at1 The atom ID (not index!) of one of the connected atoms
  \param[in] at2 The atom ID (not index!) of one of the connected atoms

  This function is a more user-frendly version of LINK_ATOMS(Atom&,Atom&)
  function. It searches for two atoms by their id and then performens the same task on
  that atoms as the other variant
*/

  int v1,v2;
  v1 = get_atom_index_by_atom_id(At1_id);
  v2 = get_atom_index_by_atom_id(At2_id);

  cout<<"in LINK_ATOMS: v1 = "<<v1<<" v2 = "<<v2<<endl;

  if((v1!=-1) && (v2!=-1)){  LINK_ATOMS(Atoms[v1],Atoms[v2]);  }
  else{
    cout<<"Error: Can not link atoms with "<<At1_id<<" and "<<At2_id<<endl;
    cout<<"Can not find one or both of the atoms with such id(s)\n";
  }
}

void System::UPDATE_FRAG_TOPOLOGY(){
/**
  This function is not to be exposed to user (or is to do this - i've not
  decided yet). Noramlly it should be used at the end of all functions 
  dealing with fragments and related topological elements.
  It simply excludes frag_bonds, frag_angles and frag_dihedrals. 
  Note it does not add these topological elements. 
*/
  int i,j;
  int a1,a2,a3,a4;
  int g1,g2,g3,g4;
  int indx1,indx2;
  vector<int>::iterator it;
  //------------------ Update bonds ------------------------
  it = Frag_bonds.begin();
  for(i=0;i<Number_of_bonds;i++){
    a1 = Bonds[i].globAtom_Index[0];
    a2 = Bonds[i].globAtom_Index[1];
    g1 = Atoms[a1].globGroup_Index;
    g2 = Atoms[a2].globGroup_Index;

    if(g1==g2){
      indx1 = Bonds[i].globGroup_Index;
      // Check if this group exist
      for(j=0;j<Number_of_frag_bonds;j++){
        indx2 = Frag_bonds[j];
        if(indx1==indx2){ // this group still exist in frag topology
                          // we need to delete it
          Frag_bonds.erase(it+j);
          Number_of_frag_bonds--;
          j--;
          break;      
        }
      }// for j
    }// g1==g2
  }// for i

  //--------------------- Update angles ------------------------
  it = Frag_angles.begin();
  for(i=0;i<Number_of_angles;i++){
    a1 = Angles[i].globAtom_Index[0];
    a2 = Angles[i].globAtom_Index[1];
    a3 = Angles[i].globAtom_Index[2];
    g1 = Atoms[a1].globGroup_Index;
    g2 = Atoms[a2].globGroup_Index;
    g3 = Atoms[a3].globGroup_Index;

    if((g1==g2)&&(g2==g3)){
      indx1 = Angles[i].globGroup_Index;
      // Check if this group exist
      for(j=0;j<Number_of_frag_angles;j++){
        indx2 = Frag_angles[j];
        if(indx1==indx2){ // this group still exist in frag topology
                          // we need to delete it
          Frag_angles.erase(it+j);
          Number_of_frag_angles--;
          j--;
          break;
        }
      }// for j
    }// g1==g2
  }// for i

  //--------------------- Update dihedrals ------------------------
  it = Frag_dihedrals.begin();
  for(i=0;i<Number_of_dihedrals;i++){
    a1 = Dihedrals[i].globAtom_Index[0];
    a2 = Dihedrals[i].globAtom_Index[1];
    a3 = Dihedrals[i].globAtom_Index[2];
    a4 = Dihedrals[i].globAtom_Index[3];
    g1 = Atoms[a1].globGroup_Index;
    g2 = Atoms[a2].globGroup_Index;
    g3 = Atoms[a3].globGroup_Index;
    g4 = Atoms[a4].globGroup_Index;

    if((g1==g2)&&(g2==g3)&&(g3==g4)){
      indx1 = Dihedrals[i].globGroup_Index;
      // Check if this group exist
      for(int j=0;j<Number_of_frag_dihedrals;j++){
        indx2 = Frag_dihedrals[j];
        if(indx1==indx2){ // this group still exist in frag topology
                          // we need to delete it
          Frag_dihedrals.erase(it+j);
          Number_of_frag_dihedrals--;
          j--;
          break;
        }
      }// for j
    }// g1==g2
  }// for i

  //--------------------- Update impropers ------------------------
  it = Frag_impropers.begin();
  for(i=0;i<Number_of_impropers;i++){
    a1 = Impropers[i].globAtom_Index[0];
    a2 = Impropers[i].globAtom_Index[1];
    a3 = Impropers[i].globAtom_Index[2];
    a4 = Impropers[i].globAtom_Index[3];
    g1 = Atoms[a1].globGroup_Index;
    g2 = Atoms[a2].globGroup_Index;
    g3 = Atoms[a3].globGroup_Index;
    g4 = Atoms[a4].globGroup_Index;

    if((g1==g2)&&(g2==g3)&&(g3==g4)){
      indx1 = Impropers[i].globGroup_Index;
      // Check if this group exist
      for(j=0;j<Number_of_frag_impropers;j++){
        indx2 = Frag_impropers[j];
        if(indx1==indx2){ // this group still exist in frag topology
                          // we need to delete it
          Frag_impropers.erase(it+j);
          Number_of_frag_impropers--;
          j--;
          break;
        }
      }// for j
    }// g1==g2
  }// for i

  //------------------ Update frag pairs ------------------------
  vector<int>::iterator it1;
  it1 = Frag_pairs.begin();
  for(i=0;i<Number_of_frag_pairs;i++){
    int fp = Frag_pairs[i];
    a1 = Pairs[fp].globAtom_Index[0];
    a2 = Pairs[fp].globAtom_Index[1];
    g1 = Atoms[a1].globGroup_Index;
    g2 = Atoms[a2].globGroup_Index;

    if(g1==g2){
      //------- Effective erasing of the i-th frag pair ------------
      int last = Number_of_frag_pairs-1;
      // Copy last pair to i-th position
      Frag_pairs[i] = Frag_pairs[last];

      // However, update some information of this pair
      // WARNING!!! [i] must be changed to something elese!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Pairs[i].globGroup_Index = i;
      Pairs[i].locGroup_Index = i; // = globGroup_Index

      // Delete last pair (just copied - avoid pair repetition)
      Frag_pairs.pop_back();
      Number_of_frag_pairs--;

      // Make one step back so to include just copied pair into consideration
      i--;
    }
  }// for i

  //----------------- Update excluded pairs ------------------------
  for(i=0;i<Number_of_pairs;i++){
    a1 = Pairs[i].globAtom_Index[0];
    a2 = Pairs[i].globAtom_Index[1];
    g1 = Atoms[a1].globGroup_Index;
    g2 = Atoms[a2].globGroup_Index;

    // Atoms belong to the same group - excluded
    if(g1==g2){
//!!!!!!!!!!!!!!!!!!!!!! WARNING : Exclusions are turned off !!!!!!!!!!!!!!!11
//      Pairs[i].Group_pair_excluded = 1;
//      Pairs[i].is_Group_pair_excluded = 1;
    }else{
      //------- Check if these atoms form interfragmental bond -------
      int res = 0;
      int b, b_indx;
      for(b=0;b<Number_of_frag_bonds;b++){
        b_indx = Frag_bonds[b];
        if((  (Bonds[b_indx].globAtom_Index[0]==a1)
            &&(Bonds[b_indx].globAtom_Index[1]==a2)
            )||
           (  (Bonds[b_indx].globAtom_Index[0]==a2)
            &&(Bonds[b_indx].globAtom_Index[1]==a1)
           )
          ){  res = 1;break;   }
      }// for b

      if(res){  
//!!!!!!!!!!!!!!!!!!!!!! WARNING : Exclusions are turned off !!!!!!!!!!!!!!!11
//        Pairs[i].Group_pair_excluded = 1;
//        Pairs[i].is_Group_pair_excluded = 1;
      }else{

      //------ Check if these atoms form interfragmental angle ----------
      for(b=0;b<Number_of_frag_angles;b++){
        b_indx = Frag_angles[b];
        if((  (Angles[b_indx].globAtom_Index[0]==a1)
            &&(Angles[b_indx].globAtom_Index[2]==a2)
            )||
           (  (Angles[b_indx].globAtom_Index[0]==a2)
            &&(Angles[b_indx].globAtom_Index[2]==a1)
           )
          ){   res = 1;break;   }
      }// for b
      if(res){
//!!!!!!!!!!!!!!!!!!!!!! WARNING : Exclusions are turned off !!!!!!!!!!!!!!!11
//        Pairs[i].Group_pair_excluded = 1;
//        Pairs[i].is_Group_pair_excluded = 1;
      }else{

      //------- Check if these atoms form interfragmental dihedral -----------
      for(b=0;b<Number_of_frag_dihedrals;b++){
        b_indx = Frag_dihedrals[b];
        if((  (Dihedrals[b_indx].globAtom_Index[0]==a1)
            &&(Dihedrals[b_indx].globAtom_Index[3]==a2)
            )||
           (  (Dihedrals[b_indx].globAtom_Index[0]==a2)
            &&(Dihedrals[b_indx].globAtom_Index[3]==a1)
           )
          ){   res = 1;break;   }
      }// for b
      if(res){
//!!!!!!!!!!!!!!!!!!!!!! WARNING : Exclusions are turned off !!!!!!!!!!!!!!!11
//        Pairs[i].Group_pair_excluded = 1;
//        Pairs[i].is_Group_pair_excluded = 1;
      }

      }// else - check dihedrals
      }// else - check angles
      }// else - atoms belond to different groups
  }// for i

}

void System::ADD_ATOM_TO_FRAGMENT(int Atom_Id,int Group_Id){
/**
  This function does nothing so far
*/
}

/*
void System::ADD_ATOM_TO_FRAGMENT(int Atom_Id,int Group_Id){
**************************************************************************
   This function adds an atom with Atom_id = Atom_Id to group of atoms
   whith the Group_id = Group_Id
   If atom is in the same group already - nothing changes
   If atom is not in this and not in any other group whith group size>1 -
   we just extend existing group by amount of 1 atom
   If atom is in different group with size >1 - size of other group decreses
   (on amount of this atom) and the size of this group increases by 1
   Overall - there should not be any atoms belonging to several different 
   groups(Fragments) at the same time
***************************************************************************
  // First - find indexes
  int i,a,g; // Atom and group indexes correspondingly
  int is_found = 0;
  int is_actually_added = 0;

  a = get_atom_index_by_atom_id(Atom_Id);
  g = get_fragment_index_by_fragment_id(Group_Id);
 
  // Now we can add atom to a group (fragment)
  int gr_of_atom = Atoms[a].globGroup_Index;     // index of the group to which atom belongs
  int mol_of_atom = Atoms[a].globMolecule_Index; // index of the molecule to which atom belongs

  if(mol_of_atom==Fragments[g].globMolecule_Index){ // Atoms which belong to different molecules 
                                                    // can not be grouped thus atom should belong
                                                    // to the same molecule to which the group belongs
    if(gr_of_atom!=g){ // If atom has already been included in the group - do nothing
                       // otherwise...
      int gr_sz = Fragments[gr_of_atom].Group_Size;

      // In any case the number of molecules does not change (because we can only)
      // group atoms belonging to the same molecule
      is_actually_added = 1;
      if(gr_sz==1){   // If that atom represents entire group
                      // In this case the number of fragments decreases by 1
                    
        // Fragments and Atoms
        //------------------------------------------------------------------------
        // consider a special case: i==g
        i = g;
        Fragments[i].globAtom_Index.push_back(a);
        Fragments[i].locAtom_Index.push_back(Fragments[i].Group_Size + 1);
        Fragments[i].Group_Size++;
        Fragments[i].globGroup_Size++;
        Fragments[i].locGroup_Size++;
        // Change number of degrees of freedom depending on number of 
        // atoms in a group after adding atom to it
        if(Fragments[i].Group_Size==2){ Nf_t += 2;    }
        else if(Fragments[i].Group_Size==3){ Nf_t += 1; }
  
        Atoms[a].locAtom_Index   = Fragments[i].Group_Size;
        Atoms[a].globGroup_Index = g;

        // Update atoms data:
        for(i=0;i<Number_of_atoms;i++){
          if(i==a){ }//Everything has been done      
          else{
            if(Atoms[i].globGroup_Index>gr_of_atom){
              Atoms[i].globGroup_Index--;
            }
          }
        }// for i

        // Update fragments data:
        for(i=gr_of_atom+1;i<Number_of_fragments;i++){
          Fragments[i].globGroup_Index--;
          Fragments[i].locGroup_Index--;
          Fragments[i].Group_id--;
        }                        
        //-------------------------------------------------------------------------------

        vector<Group>::iterator it;
        it = Fragments.begin();
        Fragments.erase(it+gr_of_atom);
        Number_of_fragments--;
        // Important!
        max_fragment_id--;
        // Decrement number of degrees of freedom
        Nf_t -= 3;

        // Molecules                                     
        vector<int>::iterator it1, it2;
        it1 = Molecules[mol_of_atom].globGroup_Index.begin();
        it2 = Molecules[mol_of_atom].locGroup_Index.begin();
        for(i=0;i<Molecules[mol_of_atom].Molecule_Size;i++){
          if(Molecules[mol_of_atom].globGroup_Index[i]==gr_of_atom){
            Molecules[mol_of_atom].globGroup_Index.erase(it1+i); 
            Molecules[mol_of_atom].locGroup_Index.erase(it2+i);
            break;
          }
        }// for i
        Molecules[mol_of_atom].globMolecule_Size--;
        Molecules[mol_of_atom].locMolecule_Size--;
        Molecules[mol_of_atom].Molecule_Size--;

        for(int m=0;m<Number_of_molecules;m++){
          for(i=0;i<Molecules[m].Molecule_Size;i++){
            if(Molecules[m].globGroup_Index[i]>gr_of_atom){
              Molecules[m].globGroup_Index[i]--;
              Molecules[m].locGroup_Index[i]--;
            }
          }
        }// for m       
      }// if gr_sz==1
      else{ // If the size of the group to which atom belongs
            // is bigger than 1
            // In this case number of fragments does not change
            // Change the number of degrees of freedom according to the size of
            // two fragments
        if(gr_sz==2){
          if(Fragments[g].Group_Size==1){
            Nf_t += 0; // does not change
          }
          else if(Fragments[g].Group_Size==2){
            Nf_t -= 1; 
          }
          else if(Fragments[g].Group_Size>=3){
            Nf_t -= 2;
          }
        }
        else if(gr_sz==3){
          if(Fragments[g].Group_Size==1){
            Nf_t += 1; 
          }
          else if(Fragments[g].Group_Size==2){
            Nf_t += 0; // does not change
          }
          else if(Fragments[g].Group_Size>=3){
            Nf_t -= 1;
          }
        }
        else if(gr_sz>3){
          if(Fragments[g].Group_Size==1){
            Nf_t += 2;
          }
          else if(Fragments[g].Group_Size==2){
            Nf_t += 1;
          }
          else if(Fragments[g].Group_Size>=3){
            Nf_t += 0; // does not change
          }
        }

        // Fragments and Atoms
        vector<int>::iterator it1,it2;
        it1 = Fragments[gr_of_atom].globAtom_Index.begin();
        it2 = Fragments[gr_of_atom].locAtom_Index.begin();
                           
        Fragments[g].globAtom_Index.push_back(a);
        Fragments[g].locAtom_Index.push_back(Fragments[g].Group_Size + 1);
        Fragments[g].Group_Size++;
        Fragments[g].globGroup_Size++;
        Fragments[g].locGroup_Size++;

        Atoms[a].locAtom_Index   = Fragments[g].Group_Size;
        Atoms[a].globGroup_Index = g;

        int loc_a;                
        for(i=0;i<Fragments[gr_of_atom].Group_Size;i++){
          if(Fragments[gr_of_atom].globAtom_Index[i]==a){
            loc_a = Fragments[gr_of_atom].locAtom_Index[i];
            Fragments[gr_of_atom].globAtom_Index.erase(it1+i);
            Fragments[gr_of_atom].locAtom_Index.erase(it2+i);
            break;
          } 
        }// for i
        Fragments[gr_of_atom].Group_Size--;
        Fragments[gr_of_atom].globGroup_Size--;
        Fragments[gr_of_atom].locGroup_Size--;

        for(i=0;i<Fragments[gr_of_atom].Group_Size;i++){
          if(Fragments[gr_of_atom].locAtom_Index[i]>loc_a){
            Fragments[gr_of_atom].locAtom_Index[i]--;
            int glob_a = Fragments[gr_of_atom].globAtom_Index[i];
            Atoms[glob_a].locAtom_Index--;
          }
        }
      }// else: gr_sz >1
    } //if gr_of_atom!=g
  }// if mol_of_atom!=Fragments[g].globMolecule_Index

}
*/

/*
void System::GROUP_ATOMS(boost::python::list atoms_list,int Group_Id){
***************************************************************
  This function groups atoms from atoms_list and makes a group
  The group just created is assigned the Group_Id id.
*****************************************************************
  for(int i=0;i<len(atoms_list);i++){
    int Atom_Id = extract<int>(atoms_list[i]);
    ADD_ATOM_TO_FRAGMENT(Atom_Id,Group_Id);
  }
  UPDATE_FRAG_TOPOLOGY();
}
*/
void System::GROUP_ATOMS(boost::python::list atoms_list,int Group_Id){
/**
  \param[in] atoms_list The Python list containing IDs (not indices!) of the atoms to be grouped
  \param[in] Group_Id The ID of the group to be created
 
  This function groups atoms from atoms_list and makes a group
  The group just created is assigned the Group_Id id.
  This gonna be improved version with more intuitive functioning
*/
  int i;
  // Get atom indexs
  vector<int> at_indxs;  
  for(i=0;i<len(atoms_list);i++){
    int Atom_Id = extract<int>(atoms_list[i]);
    int a_indx = get_atom_index_by_atom_id(Atom_Id);
    if(a_indx!=-1){
      at_indxs.push_back(a_indx);
    }
  }// for i

  // Get group index
  int g_indx = get_fragment_index_by_fragment_id(Group_Id);
 
  if(g_indx==-1){
  // No such fragment - create a new fragment
  Group fr;
  fr.globGroup_Index = Number_of_fragments;
  Number_of_fragments++;
  Fragments.push_back(fr);
  g_indx = fr.globGroup_Index;
  }
 
  // Modify existing fragment
  Group& fr = Fragments[g_indx];
  fr.Group_id = Group_Id;  fr.is_Group_id = 1;
  max_fragment_id = (Group_Id>max_fragment_id)?Group_Id:max_fragment_id;

  for(i=0;i<at_indxs.size();i++){
    Atom& at = Atoms[at_indxs[i]];

    // Erase atoms to be grouped from the previous groups
    // Atoms from the group, to which all other atoms will be
    // added are not erased
    vector<int>::iterator it;
    if(at.globGroup_Index!=g_indx){
    for(it=Fragments[at.globGroup_Index].globAtom_Index.begin();it!=Fragments[at.globGroup_Index].globAtom_Index.end();it++){
      if(*it==at.globAtom_Index){
        Fragments[at.globGroup_Index].globAtom_Index.erase(it);
        Fragments[at.globGroup_Index].Group_Size--;
        it--; // This is very important!
      }
    }
    }// if at.globGroup_Index!=g_indx

    // Add indexes of the atoms to g_indx group
    at.globGroup_Index = fr.globGroup_Index;
    if(!is_in_vector(at.globAtom_Index,fr.globAtom_Index)){
      fr.globAtom_Index.push_back(at.globAtom_Index);
    }

  }// for i
  fr.Group_Size =  fr.globAtom_Index.size();

  // Now erase groups with no atoms, change the number of fragments accordingly
  // The algorithm: Check each fragment - if the size is 0 then check
  // the rest of the fragments, starting from the end, toward the current fragment
  // If the end fragment is of zero size - pop it
  // If not - swap it and the starting fragment (it1) of the zero size. After that -
  // delete the last element of array
  // Finally, update relations of the atoms and the fragments
  i = 0;
  vector<Group>::iterator it1,it2;
  for(it1=Fragments.begin();it1!=Fragments.end();it1++,i++){
    if(it1->Group_Size==0){
      if(it1==Fragments.end()-1){ Fragments.pop_back(); Number_of_fragments--;it1--; }
      else{
        for(it2=Fragments.end()-1;it2!=it1;it2--){
          if(it2->Group_Size==0){ Fragments.pop_back(); Number_of_fragments--; }
          else{
            *it1 = *it2; 
            //-------- Update atom-fragment relations ----------
            it1->globGroup_Index = i;
            for(int j=0;j<it1->Group_Size;j++){
              Atoms[it1->globAtom_Index[j]].globGroup_Index = i;
            }
            //--------------------------------------------------
            Fragments.pop_back(); Number_of_fragments--;
            it2=it1+1;
          }// else
        }// for it2
      }// else 
    }// if Group_Size == 0
  }// for it1

  // Clean last empty places
  for(it1=Fragments.end()-1;it1!=Fragments.begin();it1--){
    if(it1->Group_Size==0){ Fragments.pop_back(); it1==Fragments.end()-1; Number_of_fragments--; }
    else{ break; }
  }

}



void System::CREATE_BONDS(boost::python::list atoms_list,boost::python::dict valence_by_element){
/**
  \param[in] atoms_list The list of the IDs of the atoms that are considered in automatic bonding assignment
  \param[in] valence_by_element Defines the valences (the maximal number of connections a given element may have) of elements

  This function will search for neighbour atoms. If the atoms
  are located on the distance smaller then the sum of their 
  radii (+/- some tolerance) and if the valence is satisfied then 
  the bond between atoms will be created. If necessary the 
  order of the bond will be assigned to value bigger then 1.
  Only atoms from the atom list will be considered.
  hyper_atoms dictionalry specifies some special cases of hypervalent atoms
  this dictionalry contains pairs: [ atom element : valence ]

*/

  vector<int> at_indexes; // Indexes of the atoms in the list, not their IDs
  vector<int> at_valences;
  // Sorted
  vector<int> at_indexes_sorted;
  vector<int> at_valences_sorted;
  int i,j;
  int id,indx;
  int grp,per;
  int sz;
  std::string blk;
  std::string elt;
  int val,max_val;  

//------- Zero step: Determine the maximal possible valence, determined by user
  max_val = 8;
  for(i=0;i<len(valence_by_element.values());i++){
    val = extract<int>(valence_by_element.values()[i]);
    if(val>max_val) {max_val = val;}
  }// for i

//-------- First step: Create array of atom indexes --------------
  sz = len(atoms_list);
  for(i=0;i<sz;i++){
    id = extract<int>(atoms_list[i]);      
    int is_found = 0;
    j = get_atom_index_by_atom_id(id);
    if(j!=-1){ at_indexes.push_back(j); }
    else{ cout<<"Error: Can not find atom with id = "<<id<<endl; }
  }

//--------- Second step: Calculate valence of the atoms -------------
  sz = at_indexes.size();
  for(i=0;i<sz;i++){
    indx = at_indexes[i];
       // Get the position of the atom in Periodic Table

//       if(os.ObjectSpace_Atoms[indx].is_Atom_group){ grp = os.ObjectSpace_Atoms[indx].Atom_group; }
//       else{ /* Alternative handling */    }

//       if(os.ObjectSpace_Atoms[indx].is_Atom_period){ per = os.ObjectSpace_Atoms[indx].Atom_period; }
//       else{ /* Alternative handling */    }
     
//       if(os.ObjectSpace_Atoms[indx].is_Atom_block){  blk = os.ObjectSpace_Atoms[indx].Atom_block;  }
//       else{ /* Alternative handling */    }

//       if(os.ObjectSpace_Atoms[indx].is_Atom_element){elt = os.ObjectSpace_Atoms[indx].Atom_element;  }
//       else{ /* Alternative handling */    }



       if(!valence_by_element.has_key(elt)){
       // Now calculate maximal valence
       if((blk == "S" )||(blk == "P")){
          val = ((grp<=4)?grp:(8-grp));
       }
       else if(blk == "D"){
          val = grp;
       }
       else if(blk == "F"){
          val = 3;
       }

       // Also take into account some other connections
//       val = val - os.ObjectSpace_Atoms[indx].globAtom_Adjacent_Atoms.size();

       }else{
       val = extract<int>(valence_by_element.get(elt));
       } 


       at_valences.push_back(val);

   }// for i

//--------- Third step: Sort valences and reindex at_indexes array ---------
   for(i=max_val;i>=1;i--){// The maximal valence is 8     
  
      for(j=0;j<at_valences.size();j++){

          if(at_valences[j]==i){
             at_valences_sorted.push_back(at_valences[j]); 
             at_indexes_sorted.push_back(at_indexes[j]);

             at_valences.erase(at_valences.begin()+j);
             at_indexes.erase(at_indexes.begin()+j);

             j--;
          }
      }

   }// for i
 
//--------- Fourth step: make pairs of atoms based on the relation between distance
//--------- and the sum of the atomic (covalent) radii of corresponding atoms

   vector<connect> all_pairs;

   double tol = 1.0;
   int indx1,indx2;
   double dist_real, dist_exp;

   sz = at_indexes_sorted.size();

   for(i=0;i<(sz-1);i++){

      connect pr;
      pr.self = at_indexes_sorted[i];

      for(j=i+1;j<sz;j++){

          indx1 = at_indexes_sorted[i];
          indx2 = at_indexes_sorted[j];

          if((at_valences_sorted[i]>=1)&&(at_valences_sorted[j]>=1)){ // Valence criterion

          dist_real = 100.0; // Default value - Very far away
//          if(os.ObjectSpace_Atoms[indx1].is_Atom_coords_cart&&os.ObjectSpace_Atoms[indx2].is_Atom_coords_cart){
//          dist_real = (os.ObjectSpace_Atoms[indx1].Atom_coords_cart - os.ObjectSpace_Atoms[indx2].Atom_coords_cart).length();
//          }

          dist_exp = 2.0;   // Default value - atoms which are separated by this much assumed to be connected
//          if(os.ObjectSpace_Atoms[indx1].is_Atom_atomic_radius&&os.ObjectSpace_Atoms[indx2].is_Atom_atomic_radius){
//          dist_exp = (os.ObjectSpace_Atoms[indx1].Atom_atomic_radius + os.ObjectSpace_Atoms[indx2].Atom_atomic_radius);
//          }

          if(dist_real<=(dist_exp + tol)){ // Distance criterion            

             pr.others.push_back(indx2);
             at_valences_sorted[i]--;
             at_valences_sorted[j]--;

          }// distance criterion

          }// valence criterion                 

      }// for j

      if(pr.others.size()>0){
         all_pairs.push_back(pr);
      }
      
   }// for i

//------- Fifth step: calculate bond orders (optional) ------------------


//-------Sixth step: Link atoms -----------------
   sz = all_pairs.size();
   int num_new_bonds = 0;
   std::cout<<"In CREATE_BONDS function: "<<std::endl;
   for(i=0;i<sz;i++){
 
     indx1 = all_pairs[i].self;

     for(j=0;j<all_pairs[i].others.size();j++){

         indx2 = all_pairs[i].others[j];

         std::cout<<"LINK atoms with indexes "<<indx1<<" and "<<indx2<<std::endl;

//         LINK_ATOMS(os,os.ObjectSpace_Atoms[indx1],os.ObjectSpace_Atoms[indx2]);
         num_new_bonds++;

     }// for j

   }// for i

//------------------------------------------------------
  
//   return num_new_bonds;

}

void System::CLONE_MOLECULE(int mol_id){
/**
  \param[in] mol_id The ID (not index!) of the molecule to be copied

  This function copies the molecule with id number = mol_id (not index
  but ID!) and all internal structure. Also it makes necessary changes to
  the structure of the object space
*/

  // First - find indexes
  int m; // Index of the molecule with mol_id
  int is_found = 0;
  int i,j,sz,indx;

  m = get_molecule_index_by_molecule_id(mol_id);
  if(m==-1){
    std::cout<<"Error: Can not duplicate (clone) molecule with mol_id = "<<mol_id<<" No such molecule\n";
    exit(97);
  }
  //================= Create corresponding atoms =====================
  int init_size = Number_of_atoms;
  int init_size1= max_fragment_id + 1;
  std::map<int,int> at_map;
  std::map<int,int> gr_map;
  std::map<int,int>::iterator it;
  int k=0;
  int k1=0;
  // Create new atoms
  for(i=0;i<init_size;i++){
    if(Atoms[i].globMolecule_Index==m){       
      CREATE_ATOM(Atoms[i]);
      at_map[Atoms[i].globAtom_Index] = Atoms[init_size + k].Atom_id;         
      it = gr_map.find(Atoms[i].globGroup_Index);
      if(it==gr_map.end()){  gr_map[Atoms[i].globGroup_Index] = init_size1 + k1;  k1++;  }
      k++;
    }
  }// for i

  // Now connenct the newly created atoms such that it resmbles the
  // connectivity in original molecule
  int at_id1,at_id2;
  for(i=0;i<init_size;i++){
    if(Atoms[i].globMolecule_Index==m){
      at_id1 = at_map[Atoms[i].globAtom_Index];
      for(j=0;j<Atoms[i].globAtom_Adjacent_Atoms.size();j++){
        indx = Atoms[i].globAtom_Adjacent_Atoms[j];
        at_id2 = at_map[Atoms[indx].globAtom_Index];
        LINK_ATOMS(at_id1,at_id2);
      }// for j
    }
  }// for i

  // Create new groups
  int gr_id;
  for(i=0;i<init_size;i++){
    if(Atoms[i].globMolecule_Index==m){
      at_id1 = at_map[Atoms[i].globAtom_Index];       
      gr_id  = gr_map[Atoms[i].globGroup_Index];
      ADD_ATOM_TO_FRAGMENT(at_id1,gr_id);
    }
  }// for i
  
}


}// namespace libchemsys
}// namespace libchemobjects
}// liblibra
