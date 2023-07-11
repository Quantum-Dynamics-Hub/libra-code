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
/**
  \file System.cpp
  \brief The file implements the basic methods of the System class
    
*/

#include "System.h"

/// liblibra namespace
namespace liblibra{

/// libchemobjects namespace
namespace libchemobjects{

/// libchemsys namespace
namespace libchemsys{


void System::init_variables(){

  max_atom_id = 0;     is_max_atom_id = 1;
  max_fragment_id = 0; is_max_fragment_id = 1;
  max_molecule_id = 0; is_max_molecule_id = 1;
  is_Box = 0;
  is_Boxold = 0;
  Box_origin = 0.0;    is_Box_origin = 1;
  dT_2 = 0.0;          is_dT_2 = 1;
//  is_TV1 = 0;
//  is_TV2 = 0;
//  is_TV3 = 0;
  
  GroupConnMatrix = NULL;  is_GroupConnMatrix = 0;
  is_AtomGraph = 0;
  
  Number_of_atoms = 0;
  Number_of_bonds = 0;
  Number_of_angles = 0;
  Number_of_dihedrals = 0; 
  Number_of_impropers = 0;
  Number_of_pairs = 0;
  Number_of_fragments = 0;
  Number_of_rings = 0;
  Number_of_molecules = 0;

  Number_of_frag_bonds = 0;
  Number_of_frag_angles = 0;
  Number_of_frag_dihedrals = 0;
  Number_of_frag_impropers = 0;
  Number_of_frag_pairs = 0;
  Number_of_surface_atoms = 0;

  Nf_t = 0;     is_Nf_t = 1;
  Nf_r = 0;     is_Nf_r = 1;
/*
  stress_opt = "fr"; is_stress_opt = 1;
  stress_at = 0.0;   is_stress_at = 0;
  stress_fr = 0.0;   is_stress_fr = 0;
  stress_ml = 0.0;   is_stress_ml = 0;
*/  
}

void System::copy_content(const System& sys){

  if(sys.is_max_atom_id){  max_atom_id = sys.max_atom_id;  is_max_atom_id = 1; }
  if(sys.is_max_fragment_id){  max_fragment_id = sys.max_fragment_id;  is_max_fragment_id = 1; }
  if(sys.is_max_molecule_id){  max_molecule_id = sys.max_molecule_id;  is_max_molecule_id = 1; }

  if(sys.is_GroupConnMatrix){
    //GroupConnMatrix = new MATRIX(sys.GroupConnMatrix->num_of_rows,sys.GroupConnMatrix->num_of_cols);  // CHECK MEMORY LEAK
    *GroupConnMatrix = *sys.GroupConnMatrix; 
    is_GroupConnMatrix = 1;
  }
  if(sys.is_AtomGraph){  AtomGraph = sys.AtomGraph;  is_AtomGraph = 1; }
  if(sys.is_Box){  Box = sys.Box;  is_Box = 1; }  
  if(sys.is_Boxold){  Boxold = sys.Boxold;  is_Boxold = 1; }
  if(sys.is_Box_origin){ Box_origin = sys.Box_origin; is_Box_origin = 1; }
  if(sys.is_dT_2){ dT_2 = sys.dT_2; is_dT_2 = 1; }

  Number_of_atoms = sys.Number_of_atoms; 
  Number_of_bonds = sys.Number_of_bonds;
  Number_of_angles = sys.Number_of_angles; 
  Number_of_dihedrals = sys.Number_of_dihedrals;
  Number_of_impropers = sys.Number_of_impropers;
  Number_of_pairs = sys.Number_of_pairs;
  Number_of_fragments = sys.Number_of_fragments;
  Number_of_rings = sys.Number_of_rings;
  Number_of_molecules = sys.Number_of_molecules;

  Atoms = sys.Atoms;
  Bonds = sys.Bonds;
  Angles = sys.Angles;
  Dihedrals = sys.Dihedrals;
  Impropers = sys.Impropers;
  Pairs = sys.Pairs;
  Fragments = sys.Fragments;
  Rings = sys.Rings;
  Molecules = sys.Molecules;

  Number_of_frag_bonds = sys.Number_of_frag_bonds;
  Number_of_frag_angles = sys.Number_of_frag_angles;
  Number_of_frag_dihedrals = sys.Number_of_frag_dihedrals;
  Number_of_frag_impropers = sys.Number_of_frag_impropers;
  Number_of_frag_pairs = sys.Number_of_frag_pairs;
  Number_of_surface_atoms = sys.Number_of_surface_atoms;

  Frag_bonds = sys.Frag_bonds;
  Frag_angles = sys.Frag_angles;
  Frag_dihedrals = sys.Frag_dihedrals;
  Frag_impropers = sys.Frag_impropers;
  Frag_pairs = sys.Frag_pairs;
  Surface_atoms = sys.Surface_atoms;

  if(sys.is_name){  name = sys.name;  is_name = 1; }
  if(sys.is_id){  id = sys.id;  is_id = 1; }
  if(sys.is_mass){ mass = sys.mass; is_mass = 1; }
  if(sys.is_Nf_t){ Nf_t = sys.Nf_t; is_Nf_t = 1; }
  if(sys.is_Nf_r){ Nf_r = sys.Nf_r; is_Nf_r = 1; }
/* 
  if(sys.is_stress_opt){ stress_opt = sys.stress_opt; is_stress_opt = 1; }
  if(sys.is_stress_at){ stress_at = sys.stress_at; is_stress_at = 1; }
  if(sys.is_stress_fr){ stress_fr = sys.stress_fr; is_stress_fr = 1; }
  if(sys.is_stress_ml){ stress_ml = sys.stress_ml; is_stress_ml = 1; }
*/
}


System::System(){
/**
  \brief Constructor
  Initialize variables to default values
*/
  // Initialize variables to default values
  init_variables();
}

System::System(const System& sys){
/**
  \brief Copy constructor
  \param[in] sys The input object
  Initialize variables to default values
  Copy content of the "sys" object which is defined
*/
  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(sys);
} 
  
System& System::operator=(const System& sys){
/**
  \brief Assignment operator
  Initialize variables to default values
  Copy content of the "sys" object which is defined
*/
  // Initialize variables to default values
  init_variables(); 
  // Copy content of th object which is defined
  copy_content(sys); 
  return *this;
} 
  
  
System::~System(){ 
/**
  \brief Destructor
  Does nothing at all
*/
}



void System::set(object at){
/** 
  \brief Set properties of the System object from an arbitrary Python object.
  \param[in] at The input object - must contain the members with the names that match the names of the internal variables.

*/

  set_value(is_name,    name,    at,"name");
  set_value(is_id,      id,      at,"id");
  set_value(is_mass,    mass,    at,"mass");
//  set_value(is_stress_opt, stress_opt, at, "stress_opt");
}

void System::show_info(){
/** 
  \brief Show System state and properties info
*/

  std::cout<<"========System information===========\n";
  std::cout<<"Number_of_atoms = "<<Number_of_atoms<<endl;
  std::cout<<"Number_of_bonds = "<<Number_of_bonds<<endl;
  std::cout<<"Number_of_angels = "<<Number_of_angles<<endl;
  std::cout<<"Number_of_dihedrals = "<<Number_of_dihedrals<<endl;
  std::cout<<"Number_of_impropers = "<<Number_of_impropers<<endl;
  std::cout<<"Number_of_pairs = "<<Number_of_pairs<<endl;
  std::cout<<"Number_of_fragments = "<<Number_of_fragments<<endl;
  std::cout<<"Number_of_rings = "<<Number_of_rings<<endl;
  std::cout<<"Number_of_molecules = "<<Number_of_molecules<<endl;

//  std::cout<<"Number_of_frag_bonds = "<<Number_of_frag_bonds<<endl;
//  std::cout<<"Number_of_frag_angles = "<<Number_of_frag_angles<<endl;
//  std::cout<<"Number_of_frag_dihedrals = "<<Number_of_frag_dihedrals<<endl;
//  std::cout<<"Number_of_frag_pairs = "<<Number_of_frag_pairs<<endl;


  if(is_Box)    {std::cout<<" Box = "<<Box<<std::endl; }
  if(is_Boxold) {std::cout<<" Boxold = "<<Boxold<<std::endl; }
  if(is_name)   {std::cout<<" name = "<<name<<std::endl;   }
  if(is_id)     {std::cout<<" id = "<<id<<std::endl;     }
  if(is_mass)   {std::cout<<" mass = "<<mass<<std::endl;   }
/*
  if(is_stress_opt){std::cout<<" stress_opt = "<<stress_opt<<std::endl; }
  if(is_stress_at) {std::cout<<" stress_at = "<<stress_at<<std::endl; }
  if(is_stress_fr) {std::cout<<" stress_fr = "<<stress_fr<<std::endl; }
  if(is_stress_ml) {std::cout<<" stress_ml = "<<stress_ml<<std::endl; }
*/
}


}// namespace libchemsys
}// namespace libchemobjects
}// liblibra
