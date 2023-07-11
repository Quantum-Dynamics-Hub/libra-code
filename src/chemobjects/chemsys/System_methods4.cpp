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
  \file System_methods4.cpp
  \brief The file implements the function for determining functional group name for each given atom.
    
*/

#include "System.h"

/// liblibra namespace
namespace liblibra{


/// libchemobjects namespace
namespace libchemobjects{

/// libchemsys namespace
namespace libchemsys{


int System::is(std::string X, int n, int m,int sz,Atom** at,vector<Atom*>& which){
/**
  This auxiliary function searches for an atom of element X with
  n connections in array of pointers to Atom objects "at", which is
  of size sz. 
  It returnes 1 if m such atoms was found in array, 0 otherwise
  The vector of matching atoms (addresses) is formed in the "which" variable
  Note: for the special cases like X = OH the adresses of next neighbors
  (in this case of H atom) will also be added to the "which" variable
*/

  if(which.size()>0){ which.clear(); }
  int res = 0;
  for(int i=0;i<sz;i++){
    if(X=="OH"){
      if(at[i]->Atom_element=="O" && at[i]->globAtom_Adjacent_Atoms.size()==2){
        Atom& at01 = Atoms[at[i]->globAtom_Adjacent_Atoms[0]];
        Atom& at02 = Atoms[at[i]->globAtom_Adjacent_Atoms[1]];
        if(at01.Atom_element=="H"){
          res += 1;
          which.push_back(at[i]);
          which.push_back(&at01);
        }
        else if(at02.Atom_element=="H"){
          res += 1;
          which.push_back(at[i]);
          which.push_back(&at02);
        }
      }
    }// X = OH
    else if(X=="NH"){
      if(at[i]->Atom_element=="N" && at[i]->globAtom_Adjacent_Atoms.size()==2){
        Atom& at01 = Atoms[at[i]->globAtom_Adjacent_Atoms[0]];
        Atom& at02 = Atoms[at[i]->globAtom_Adjacent_Atoms[1]];
        if(at01.Atom_element=="H"){
          res += 1;
          which.push_back(at[i]);
          which.push_back(&at01);
        }
        else if(at02.Atom_element=="H"){
          res += 1;
          which.push_back(at[i]);
          which.push_back(&at02);
        }
      }
    }// X = NH
    else{
      if(at[i]->Atom_element==X && at[i]->globAtom_Adjacent_Atoms.size()==n){ 
        res += 1;
        which.push_back(at[i]);
      }
    }// else X - any except for OH
  }
  res = (res==m)?1:0;  
  if(res==0){ which.clear(); }
  return res;
}

void System::store(vector<Atom*>& w,vector<Atom*>& W){
  if(W.size()>0){ W.clear(); }
  for(int i=0;i<w.size();i++){ W.push_back(w[i]); }
  if(w.size()>0){ w.clear();}
}

void System::set(vector<Atom*>& at,std::string grp_name){
/**
  This auxiliary function sets the Atom_functional_group property
  of the atoms at[which[i]] to the value of grp_name and sets the flag
  that it is defined
*/

  for(int i=0;i<at.size();i++){
    at[i]->Atom_functional_group = grp_name;
    at[i]->is_Atom_functional_group = 1;
  }
}

void System::determine_functional_groups(int assign_rings){
/**
  This function uses the information on atom connectivity 
  stored in AtomGraph variable (of type GRAPH<Atom*,Group*>)
  and determines to which functional group (from standard
  list of organic and inorganic functional groups) each atom
  belongs. This information is then stored in functional_groups
  variable of type (map<int,std::string>) where the int - is
  the index of the atom and std::string is the name of the 
  functional group
  Functional groups are defined in:
  http://en.wikipedia.org/wiki/Functional_group
 
  additional groups are from:
  http://www.thechemblog.com/?p=375

  http://en.wikipedia.org/wiki/Sulfonamide_%28chemistry%29   # sulfonamide

  order: most general - first
 
  Total list of functional groups:
  Alkyl
    Ketal
    Hemiketal
    Acetal
    Hemiacetal
    Orthoester
  Alkenyl
    Carbonyl
      Aldehyde 
      Haloformyl
    Carbonate_ester
    Carboxylate
    Ester
      Carboxyl
    Carboxamide
    Secondary_ketimine
    Primary_ketimine
      Secondary_aldimine
      Primary_aldimine
    Thione
      Thial
  Alkynyl
    Nitrile
    Cyanate
    Isocyanate
    Thiocyanate
    Isothiocyanate
  Isonitrile
  Ammoniun
  Tertiary_amine
  Secondary_amine
  Primary_amine
  Imide
  Nito
  Nitrate
  Azo
  Azide
  Nitroso
  Nitrosooxy
  Thioether
    Sulfhydril
    Disulphide
  Sulfinyl
  Sulfino
  Sulfonamide
  Phosphino
  Phosphono
  Phosphate
  Phosphodiester
  Ether
    Hydroxyl
    Water
    Peroxy
      Hydroperoxy
 
// Hydrogen functional groups!
  Sulfonamide_H_on_N
  Sulfonamide_H_on_S
  Hemiacetal_H_on_O
  Hemiacetal_H_on_C
  Primary_aldimine_H_on_N
  Primary_aldimine_H_on_C
  Imide_H_on_N
  Imide_H_on_C
  Carboxamide_H_on_N
  Carboxamide_H_on_C
*/
  int i,j;
  cout<<"In System::determine_functional_groups()\n";

  if(assign_rings==1){    Assign_Rings();   }

  // Set the Atom_functional_group property for all atoms
  // as undefined
  for(i=0;i<Number_of_atoms;i++){
    Atoms[i].is_Atom_functional_group = 0;
  }
      
  //111111111 First determine carbon-based groups 111111111111111111
      
  for(i=0;i<Number_of_atoms;i++){
    vector<Atom*> w1,w2,w3,W1,W2,W3,all;
    int sz1 = Atoms[i].globAtom_Adjacent_Atoms.size(); 
    std::string func_grp;
    Atom** at; // adjacent atoms
    at = new Atom*[sz1];
    for(j=0;j<sz1;j++){ at[j] = &Atoms[Atoms[i].globAtom_Adjacent_Atoms[j]]; }

    //--------------- Carbon --------------------------
    if(Atoms[i].Atom_element=="C" && !Atoms[i].is_Atom_functional_group){

      if(sz1==4){
        func_grp = "Alkyl";
        if(is("O",2,2,sz1,at,w1)){  func_grp = "Ketal"; store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",2,1,sz1,at,w1) && is("OH",2,1,sz1,at,w2)){  func_grp = "Hemiketal"; store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",2,2,sz1,at,w1) && is("H",1,1,sz1,at,w2)) {  func_grp = "Acetal"; store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",2,1,sz1,at,w1) && is("OH",2,1,sz1,at,w2) && is("H",1,1,sz1,at,w3)){ func_grp = "Hemiacetal"; store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",2,3,sz1,at,w1)){ func_grp = "Orthoester"; store(w1,W1);store(w2,W2);store(w3,W3); }
      }
 
      else if(sz1==3){
        func_grp = "Alkenyl";
        if(is("O",1,1,sz1,at,w1)){ func_grp = "Carbonyl"; store(w1,W1);store(w2,W2);store(w3,W3);}
          if(is("O",1,1,sz1,at,w1) && is("H",1,1,sz1,at,w2)){ func_grp = "Aldehyde"; store(w1,W1);store(w2,W2);store(w3,W3);}
          if(is("O",1,1,sz1,at,w1) && is("F",1,1,sz1,at,w2)){ func_grp = "Haloformyl"; store(w1,W1);store(w2,W2);store(w3,W3);}
          if(is("O",1,1,sz1,at,w1) && is("Cl",1,1,sz1,at,w2)){ func_grp = "Haloformyl"; store(w1,W1);store(w2,W2);store(w3,W3);}
          if(is("O",1,1,sz1,at,w1) && is("Br",1,1,sz1,at,w2)){ func_grp = "Haloformyl"; store(w1,W1);store(w2,W2);store(w3,W3);}
          if(is("O",1,1,sz1,at,w1) && is("I",1,1,sz1,at,w2)){ func_grp = "Haloformyl"; store(w1,W1);store(w2,W2);store(w3,W3);}

        if(is("O",1,1,sz1,at,w1) && is("O",2,2,sz1,at,w2)){ func_grp = "Carbonate_ester"; store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",1,2,sz1,at,w1)){ func_grp = "Carboxylate"; store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",1,1,sz1,at,w1) && is("O",2,1,sz1,at,w2)){ func_grp = "Ester"; store(w1,W1);store(w2,W2);store(w3,W3);}
          if(is("O",1,1,sz1,at,w1) && is("OH",2,1,sz1,at,w2)){ func_grp = "Carboxyl"; store(w1,W1);store(w2,W2);store(w3,W3);}

        if(is("O",1,1,sz1,at,w1) && is("N",3,1,sz1,at,w2)){ func_grp = "Carboxamide"; store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("N",2,1,sz1,at,w1)){ func_grp = "Secondary_ketimine"; store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("NH",2,1,sz1,at,w1)){ func_grp = "Primary_ketimine"; store(w1,W1);store(w2,W2);store(w3,W3);}
          if(is("N",2,1,sz1,at,w1) && is("H",1,1,sz1,at,w2)){ func_grp = "Secondary_aldimine"; store(w1,W1);store(w2,W2);store(w3,W3);}
          if(is("NH",2,1,sz1,at,w1)&& is("H",1,1,sz1,at,w2)){ func_grp = "Primary_aldimine"; store(w1,W1);store(w2,W2);store(w3,W3);}

        if(is("S",1,1,sz1,at,w1)){ func_grp = "Thione"; store(w1,W1);store(w2,W2);store(w3,W3);}
          if(is("S",1,1,sz1,at,w1) && is("H",1,1,sz1,at,w2)){ func_grp = "Thial"; store(w1,W1);store(w2,W2);store(w3,W3);}

      }// C with 3 connections

      else if(sz1==2){
        func_grp = "Alkynyl";
        if(is("N",1,1,sz1,at,w1)){ func_grp = "Nitrile";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("N",1,1,sz1,at,w1) && is("O",2,1,sz1,at,w2)){ func_grp = "Cyanate";   store(w1,W1);store(w2,W2);store(w3,W3); }
        if(is("N",2,1,sz1,at,w1) && is("O",1,1,sz1,at,w2)){ func_grp = "Isocyanate";store(w1,W1);store(w2,W2);store(w3,W3); }       
        if(is("N",1,1,sz1,at,w1) && is("S",2,1,sz1,at,w2)){ func_grp = "Thiocyanate";   store(w1,W1);store(w2,W2);store(w3,W3); }
        if(is("N",2,1,sz1,at,w1) && is("S",1,1,sz1,at,w2)){ func_grp = "Isothiocyanate";store(w1,W1);store(w2,W2);store(w3,W3); }


      }// C with 2 connections

      else if(sz1==1){
        if(is("N",2,1,sz1,at,w1)){ func_grp = "Isonitrile";   }

      }// C with 1 connection

    all.push_back(&Atoms[i]);
    for(j=0;j<W1.size();j++) { all.push_back(W1[j]); }
    for(j=0;j<W2.size();j++) { all.push_back(W2[j]); }
    for(j=0;j<W3.size();j++) { all.push_back(W3[j]); }

    set(all,func_grp);

    }// carbon
  
  }// for all atoms
  //1111111111111111111111111111111111111111111111111111111111111111

  //2222222 Second determine remaining nitrogen-based groups 2222222222

  for(i=0;i<Number_of_atoms;i++){
    vector<Atom*> w1,w2,w3,W1,W2,W3,all;
    int sz1 = Atoms[i].globAtom_Adjacent_Atoms.size();
    std::string func_grp;
    Atom** at; // adjacent atoms
    at = new Atom*[sz1];
    for(j=0;j<sz1;j++){ at[j] = &Atoms[Atoms[i].globAtom_Adjacent_Atoms[j]]; }

    //--------------- Remaining nitrogens --------------------------
    if(Atoms[i].Atom_element=="N" && !Atoms[i].is_Atom_functional_group){

      if(sz1==4){ func_grp="Ammonium";   }

      else if(sz1==3){
        func_grp = "Tertiary_amine";
        if(is("H",1,1,sz1,at,w1)){ func_grp = "Secondary_amine";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("H",1,2,sz1,at,w1)){ func_grp = "Primary_amine";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("C",1,2,sz1,at,w1)){ 
          // Check for imide explicitly
          if(w1[0]->is_Atom_functional_group && w1[1]->is_Atom_functional_group){          
            if(w1[0]->Atom_functional_group=="Carboxamide" && w1[1]->Atom_functional_group=="Carboxamide"){
              func_grp = "Imide";  store(w1,W1);store(w2,W2);store(w3,W3);
            }
          }
        }
        if(is("O",1,2,sz1,at,w1)){ func_grp = "Nitro";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",1,2,sz1,at,w1) && is("O",2,1,sz1,at,w2)){ func_grp = "Nitrate";  store(w1,W1);store(w2,W2);store(w3,W3);}

      }// N with 3 connections

      else if(sz1==2){

        if(is("N",2,1,sz1,at,w1) ){ func_grp = "Azo";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("N",1,1,sz1,at,w1) && is("N",2,1,sz1,at,w2)){ func_grp = "Azide";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",1,1,sz1,at,w1) ){ func_grp = "Nitroso";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",1,1,sz1,at,w1) && is("O",2,1,sz1,at,w2)){ func_grp = "Nitrosooxy";  store(w1,W1);store(w2,W2);store(w3,W3);}

      }
   
    all.push_back(&Atoms[i]);
    for(j=0;j<W1.size();j++) { all.push_back(W1[j]); }
    for(j=0;j<W2.size();j++) { all.push_back(W2[j]); }
    for(j=0;j<W3.size();j++) { all.push_back(W3[j]); }

    set(all,func_grp);

    }// remainig nitrogens

  }// for all atoms
  //22222222222222222222222222222222222222222222222222222222222222222222


  //3333333 Third determine remaining sulfur-based groups 33333333333

  for(i=0;i<Number_of_atoms;i++){
    vector<Atom*> w1,w2,w3,W1,W2,W3,all;
    int sz1 = Atoms[i].globAtom_Adjacent_Atoms.size();
    std::string func_grp;
    Atom** at; // adjacent atoms
    at = new Atom*[sz1];
    for(j=0;j<sz1;j++){ at[j] = &Atoms[Atoms[i].globAtom_Adjacent_Atoms[j]]; }

    //--------------- Remaining sulfurs --------------------------
    if(Atoms[i].Atom_element=="S" && !Atoms[i].is_Atom_functional_group){

      if(sz1==2){
        func_grp = "Thioether";
        if(is("H",1,1,sz1,at,w1) ){ func_grp = "Sulfhydryl";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("S",2,1,sz1,at,w1) ){ func_grp = "Disulfide";  store(w1,W1);store(w2,W2);store(w3,W3);}

      }// sulfur with 2 connections

      else if(sz1==3){
        if(is("O",1,1,sz1,at,w1) ){ func_grp = "Sulfinyl";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",1,1,sz1,at,w1) && is("OH",1,1,sz1,at,w2)){ func_grp = "Sulfino";  store(w1,W1);store(w2,W2);store(w3,W3);}

      }// sulfur with 3 connections

      else if(sz1==4){
        if(is("O",1,2,sz1,at,w1) ){ func_grp = "Sulfonyl";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",1,2,sz1,at,w1) && is("OH",1,1,sz1,at,w2)){ func_grp = "Sulfo";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",1,2,sz1,at,w1) && is("N",3,1,sz1,at,w2)){ func_grp = "Sulfonamide";  store(w1,W1);store(w2,W2);store(w3,W3);}

      }// sulfur with 4 connections

    all.push_back(&Atoms[i]);
    for(j=0;j<W1.size();j++) { all.push_back(W1[j]); }
    for(j=0;j<W2.size();j++) { all.push_back(W2[j]); }
    for(j=0;j<W3.size();j++) { all.push_back(W3[j]); }

    set(all,func_grp);

    }// remainig sulfurs

  }// for all atoms
  //33333333333333333333333333333333333333333333333333333333333333

  //444444444 Fourth determine remaining phosphorus-based groups 4444444444444

  for(i=0;i<Number_of_atoms;i++){
    vector<Atom*> w1,w2,w3,W1,W2,W3,all;
    int sz1 = Atoms[i].globAtom_Adjacent_Atoms.size();
    std::string func_grp;
    Atom** at; // adjacent atoms
    at = new Atom*[sz1];
    for(j=0;j<sz1;j++){ at[j] = &Atoms[Atoms[i].globAtom_Adjacent_Atoms[j]]; }

    //--------------- Remaining phosphorus --------------------------
    if(Atoms[i].Atom_element=="P" && !Atoms[i].is_Atom_functional_group){

      if(sz1==3){    func_grp = "Phosphino"; }
      if(sz1==4){
        if(is("O",1,1,sz1,at,w1)&& is("OH",1,2,sz1,at,w2) ){ func_grp = "Phosphono";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",1,1,sz1,at,w1)&& is("OH",1,2,sz1,at,w2) && is("O",2,1,sz1,at,w3)){ func_grp = "Phosphate";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",1,1,sz1,at,w1)&& is("OH",1,1,sz1,at,w2) && is("O",2,2,sz1,at,w3)){ func_grp = "Phosphodiester";  store(w1,W1);store(w2,W2);store(w3,W3);}

      }// phosphorus with 4 connections

    all.push_back(&Atoms[i]);
    for(j=0;j<W1.size();j++) { all.push_back(W1[j]); }
    for(j=0;j<W2.size();j++) { all.push_back(W2[j]); }
    for(j=0;j<W3.size();j++) { all.push_back(W3[j]); }

    set(all,func_grp);

    }// remainig phosphorus

  }// for all atoms
  //44444444444444444444444444444444444444444444444444444444444444444444444444


  //555555555 Fifth determine remaining oxygen-based groups 5555555555555555

  for(i=0;i<Number_of_atoms;i++){
    vector<Atom*> w1,w2,w3,W1,W2,W3,all;
    int sz1 = Atoms[i].globAtom_Adjacent_Atoms.size();
    std::string func_grp;
    Atom** at; // adjacent atoms
    at = new Atom*[sz1];
    for(j=0;j<sz1;j++){ at[j] = &Atoms[Atoms[i].globAtom_Adjacent_Atoms[j]]; }

    //--------------- Remaining oxygens --------------------------
    if(Atoms[i].Atom_element=="O" && !Atoms[i].is_Atom_functional_group){

      if(sz1==2){
        func_grp = "Ether";
        if(is("H",1,1,sz1,at,w1) ){ func_grp = "Hydroxyl";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("H",1,2,sz1,at,w1) ){ func_grp = "Water";  store(w1,W1);store(w2,W2);store(w3,W3);}
        if(is("O",2,1,sz1,at,w1) ){ func_grp = "Peroxy";  store(w1,W1);store(w2,W2);store(w3,W3);}
          if(is("OH",2,1,sz1,at,w1) ){ func_grp = "Hydroperoxy";  store(w1,W1);store(w2,W2);store(w3,W3);}

      }// oxygen with 2 connections

    all.push_back(&Atoms[i]);
    for(j=0;j<W1.size();j++) { all.push_back(W1[j]); }
    for(j=0;j<W2.size();j++) { all.push_back(W2[j]); }
    for(j=0;j<W3.size();j++) { all.push_back(W3[j]); }

    set(all,func_grp);

    }// remainig oxygens

  }// for all atoms
  //55555555555555555555555555555555555555555555555555555555555555555555555

  //6666666 Sixth asisgn the hydrogens with the group names 666666666666666

  for(i=0;i<Number_of_atoms;i++){
    vector<Atom*> w1,w2,w3,W1,W2,W3,all;
    int sz1 = Atoms[i].globAtom_Adjacent_Atoms.size();
    std::string func_grp;
    Atom** at; // adjacent atoms
    at = new Atom*[sz1];
    for(j=0;j<sz1;j++){ at[j] = &Atoms[Atoms[i].globAtom_Adjacent_Atoms[j]]; }

    //--------------- Hydrogen --------------------------
    if(Atoms[i].Atom_element=="H" /*&& !Atoms[i].is_Atom_functional_group*/){

      if(sz1==1){
        func_grp = Atoms[Atoms[i].globAtom_Adjacent_Atoms[0]].Atom_functional_group;
        if(is("H",1,1,sz1,at,w1) ){ func_grp = "Hydrogen";  store(w1,W1);store(w2,W2);store(w3,W3);}

        // Break ambiguity for hydrogens connected to different atoms in the same functional group
        if(func_grp=="Sulfonamide"){ 
          if(is("N",3,1,sz1,at,w1) ){ func_grp = "Sulfonamide_H_on_N"; }
          if(is("S",4,1,sz1,at,w1) ){ func_grp = "Sulfonamide_H_on_S"; }
        }
        else if(func_grp=="Hemiacetal"){
          if(is("O",2,1,sz1,at,w1) ){ func_grp = "Hemiacetal_H_on_O";  }
          if(is("C",4,1,sz1,at,w1) ){ func_grp = "Hemiacetal_H_on_C";  }
        }
        else if(func_grp=="Primary_aldimine"){
          if(is("N",2,1,sz1,at,w1) ){ func_grp = "Primary_aldimine_H_on_N";  }
          if(is("C",2,1,sz1,at,w1) ){ func_grp = "Primary_aldimine_H_on_C";  }
        }
        else if(func_grp=="Imide"){
          if(is("N",3,1,sz1,at,w1) ){ func_grp = "Imide_H_on_N";  }
          if(is("C",3,1,sz1,at,w1) ){ func_grp = "Imide_H_on_C";  }
        }
        else if(func_grp=="Carboxamide"){
          if(is("N",3,1,sz1,at,w1) ){ func_grp = "Carboxamide_H_on_N";  }
          if(is("C",3,1,sz1,at,w1) ){ func_grp = "Carboxamide_H_on_C";  }
        }


      }// hydrogen with 1 connection

    all.push_back(&Atoms[i]);
    for(j=0;j<W1.size();j++) { all.push_back(W1[j]); }
    for(j=0;j<W2.size();j++) { all.push_back(W2[j]); }
    for(j=0;j<W3.size();j++) { all.push_back(W3[j]); }

    set(all,func_grp);

    }// remainig oxygens

  }// for all atoms
  //55555555555555555555555555555555555555555555555555555555555555555555555

  for(i=0;i<Number_of_atoms;i++){    
    cout<<"Atom i="<<i;
    if(Atoms[i].is_Atom_functional_group){
      cout<<" belongs to the functional group= "<<Atoms[i].Atom_functional_group<<endl;
    }
    else{
      cout<<" does not have functional group assigned to it\n";
    }
  }// for i


}

}// namespace libchemsys
}// namespace libchemobjects
}// liblibra



