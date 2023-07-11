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
  \file Model_Parameters.cpp This file implements the functions and class methods for setting up
  parameters of atomistic Hamiltonians

*/

#include "Model_Parameters.h"

/// liblibra namespace
namespace liblibra{



namespace libmodel_parameters{



// HF_integrals class:
int HF_integrals::find_data(int a,int b,int c,int d){
  int i,j,sz; sz = data.size();
  j = -1;
  for(i=0;i<sz;i++){
    if(data[i].a==a){
      if(data[i].b==b){
        if(data[i].c==c){
          if(data[i].d==d){
            j = i; break;
          }
        }
      }
    }
  }// for i

  return j;
}

void HF_integrals::set_JK_values(int a,int b, int c, int d, double J, double K){
  int i = find_data(a,b,c,d);
  if(i>-1){ 
    data[i].J_abcd = J;
    data[i].K_adcb = K;
  }
  else{ 
    data_element x;
    x.a = a; x.b = b; x.c = c; x.d = d;
    x.J_abcd = J; x.K_adcb = K;
    data.push_back(x); 
  }
}

void HF_integrals::get_JK_values(int a,int b, int c, int d, double& J, double& K){

  int i = find_data(a,b,c,d);
  J = 0.0; K = 0.0;
  if(i>-1){ J = data[i].J_abcd; K = data[i].K_adcb; }

}







// EHT_K class:
int EHT_K::find_data(std::string elt1,std::string orb_type1,std::string elt2,std::string orb_type2){
  int i,j,sz; sz = data.size();
  j = -1;

  for(i=0;i<sz;i++){
/*
    if((data[i].elt1==elt1 && data[i].orb_type1==orb_type1 && data[i].elt2==elt2 && data[i].orb_type2==orb_type2)||
       (data[i].elt2==elt1 && data[i].orb_type2==orb_type1 && data[i].elt1==elt2 && data[i].orb_type1==orb_type2)){ j = i; break; }
*/
    int status = 0;

    if(data[i].elt1==elt1){
      if(data[i].orb_type1==orb_type1){
        if(data[i].elt2==elt2){
          if(data[i].orb_type2==orb_type2){
            status = 1;
          }
        }
      }
    }
    if(!status){  // try symmetric version
      if(data[i].elt2==elt1){
        if(data[i].orb_type2==orb_type1){
          if(data[i].elt1==elt2){
            if(data[i].orb_type1==orb_type2){
              status = 1;
            }
          }
        }
      }
    }
    if(status){ j = i; break;}


  } // for i
  return j;
}


int EHT_K::find_data(std::string elt1,std::string orb_type1){
  int i,j,sz; sz = pp_data.size();
  j = -1;
  for(i=0;i<sz;i++){
    int status = 0;

    if(pp_data[i].elt1==elt1){
      if(pp_data[i].orb_type1==orb_type1){

            status = 1;
      }
    }

    if(status){ j = i; break;}

  } // for i
  return j;
}

int EHT_K::find_data(int n1_,int n2_,int n3_,int n4_){
  int i,j,sz; sz = psps_data.size();
  j = -1;
  for(i=0;i<sz;i++){
    int status = 0;

    if(psps_data[i].n1==n1_){
      if(psps_data[i].n2==n2_){
        if(psps_data[i].n3==n3_){
          if(psps_data[i].n4==n4_){
            status = 1;
          }
        }
      }
    }
    if(status){ j = i; break;}

  } // for i
  return j;
}



void EHT_K::set_PSPS_value(int n1_,int n2_,int n3_,int n4_,double val){

  int i = find_data(n1_,n2_,n3_,n4_);
  if(i>-1){ psps_data[i].K = val; }
  else{  psps_data_element x; x.n1 = n1_; x.n2 = n2_; x.n3 = n3_; x.n4 = n4_; x.K = val;  psps_data.push_back(x);  }


}



void EHT_K::set_PPa_value(std::string elt1, std::string orb_type1, double val){

  int i = find_data(elt1,orb_type1);
  if(i>-1){ pp_data[i].PPa_value = val; }
  else{  pp_data_element x; x.elt1 = elt1; x.orb_type1 = orb_type1; x.PPa_value = val;  pp_data.push_back(x);  }
  
}

void EHT_K::set_PP0_value(std::string elt1, std::string orb_type1, double val){

  int i = find_data(elt1,orb_type1);
  if(i>-1){ pp_data[i].PP0_value = val; }
  else{  pp_data_element x; x.elt1 = elt1; x.orb_type1 = orb_type1; x.PP0_value = val;  pp_data.push_back(x);  }

}

void EHT_K::set_PP1_value(std::string elt1, std::string orb_type1, double val){

  int i = find_data(elt1,orb_type1);
  if(i>-1){ pp_data[i].PP1_value = val; }
  else{  pp_data_element x; x.elt1 = elt1; x.orb_type1 = orb_type1; x.PP1_value = val;  pp_data.push_back(x);  }

}

void EHT_K::set_PP2_value(std::string elt1, std::string orb_type1, double val){

  int i = find_data(elt1,orb_type1);
  if(i>-1){ pp_data[i].PP2_value = val; }
  else{  pp_data_element x; x.elt1 = elt1; x.orb_type1 = orb_type1; x.PP2_value = val;  pp_data.push_back(x);  }

}



double EHT_K::get_PPa_value(std::string elt1,std::string orb_type1){  
  int i = find_data(elt1,orb_type1);
  double res = 0.0;  if(i>-1){ res = pp_data[i].PPa_value; }
  return res;
}

double EHT_K::get_PP0_value(std::string elt1,std::string orb_type1){  
  int i = find_data(elt1,orb_type1);
  double res = 0.0;  if(i>-1){ res = pp_data[i].PP0_value; }
  return res;
}

double EHT_K::get_PP1_value(std::string elt1,std::string orb_type1){  
  int i = find_data(elt1,orb_type1);
  double res = 0.0;  if(i>-1){ res = pp_data[i].PP1_value; }
  return res;
}

double EHT_K::get_PP2_value(std::string elt1,std::string orb_type1){  
  int i = find_data(elt1,orb_type1);
  double res = 0.0;  if(i>-1){ res = pp_data[i].PP2_value; }
  return res;
}





 
void EHT_K::set_K_value(int k_indx,std::string elt1,std::string orb_type1,std::string elt2,std::string orb_type2, double K){
/**
  This function sets the K parameters of the (generalized) EHT models for different types of pair interactions
  If the pair of orbital is not found, the record will be created. If the pair exists, only the parameter will be changed

  \param[in] k_indx The index of the paramter. So far can be: 0, 1, 2, 3, 4 
  \param[in] elt1 The element name of the first atom in the interacting pair (e.g. H, Si, etc.)
  \param[in] orb_type1 The orbital type for the first atom (orbital) in the interacting pair (1s, 2p, etc.)
  \param[in] elt2 The element name of the second atom in the interacting pair (e.g. H, Si, etc.)
  \param[in] orb_type2 The orbital type for the second atom (orbital) in the interacting pair (1s, 2p, etc.)
  \param[in] K The value to which the corresponding parameter will be set

*/

  int i = find_data(elt1,orb_type1,elt2,orb_type2);
  if(i>-1){ data[i].K_value[k_indx] = K;  data[i].is_K_value[k_indx] = 1; }
  else{ 
    data_element x; 
    x.elt1 = elt1; x.orb_type1 = orb_type1; x.elt2 = elt2; x.orb_type2 = orb_type2; x.K_value[k_indx] = K; x.is_K_value[k_indx] = 1;
    data.push_back(x); 
  }
}

double EHT_K::get_K_value(int k_indx,std::string elt1,std::string orb_type1,std::string elt2,std::string orb_type2){  
/**
  This function returns the K parameters of the (generalized) EHT models for different types of pair interactions
  If the pair of orbital is not found, the defaul value will be returned. 

  \param[in] k_indx The index of the paramter. So far can be: 0, 1, 2, 3, 4 
  \param[in] elt1 The element name of the first atom in the interacting pair (e.g. H, Si, etc.)
  \param[in] orb_type1 The orbital type for the first atom (orbital) in the interacting pair (1s, 2p, etc.)
  \param[in] elt2 The element name of the second atom in the interacting pair (e.g. H, Si, etc.)
  \param[in] orb_type2 The orbital type for the second atom (orbital) in the interacting pair (1s, 2p, etc.)

*/

  int i = find_data(elt1,orb_type1,elt2,orb_type2);
  double res = K_default[k_indx]; 
  if(i>-1){ 
    if(data[i].is_K_value[k_indx] == 1){ res = data[i].K_value[k_indx];}
  }

  return res;
}

void EHT_K::set_C_value(int c_indx,std::string elt1,std::string orb_type1,std::string elt2,std::string orb_type2, double C){
/**
  This function sets the C parameters of the (generalized) EHT models for different types of pair interactions
  If the pair of orbital is not found, the record will be created. If the pair exists, only the parameter will be changed

  \param[in] c_indx The index of the paramter. So far can be: 0, 1, 2, 3, 4 
  \param[in] elt1 The element name of the first atom in the interacting pair (e.g. H, Si, etc.)
  \param[in] orb_type1 The orbital type for the first atom (orbital) in the interacting pair (1s, 2p, etc.)
  \param[in] elt2 The element name of the second atom in the interacting pair (e.g. H, Si, etc.)
  \param[in] orb_type2 The orbital type for the second atom (orbital) in the interacting pair (1s, 2p, etc.)
  \param[in] C The value to which the corresponding parameter will be set

*/

  int i = find_data(elt1,orb_type1,elt2,orb_type2);
  if(i>-1){ data[i].C_value[c_indx] = C;  data[i].is_C_value[c_indx] = 1; }
  else{ 
    data_element x; 
    x.elt1 = elt1; x.orb_type1 = orb_type1; x.elt2 = elt2; x.orb_type2 = orb_type2; x.C_value[c_indx] = C; x.is_C_value[c_indx] = 1;
    data.push_back(x); 
  }
}

double EHT_K::get_C_value(int c_indx,std::string elt1,std::string orb_type1,std::string elt2,std::string orb_type2){  
/**
  This function returns the C parameters of the (generalized) EHT models for different types of pair interactions
  If the pair of orbital is not found, the defaul value will be returned. 

  \param[in] c_indx The index of the paramter. So far can be: 0, 1, 2, 3, 4 
  \param[in] elt1 The element name of the first atom in the interacting pair (e.g. H, Si, etc.)
  \param[in] orb_type1 The orbital type for the first atom (orbital) in the interacting pair (1s, 2p, etc.)
  \param[in] elt2 The element name of the second atom in the interacting pair (e.g. H, Si, etc.)
  \param[in] orb_type2 The orbital type for the second atom (orbital) in the interacting pair (1s, 2p, etc.)

*/

  int i = find_data(elt1,orb_type1,elt2,orb_type2);
  double res = C_default[c_indx]; 
  if(i>-1){ 
    if(data[i].is_C_value[c_indx] == 1){ res = data[i].C_value[c_indx];}
  }

  return res;
}


void EHT_K::show(){
/**
  Prints all the orbital pair records existing in the EHT_K object
  Careful: this may be a very long listing
*/

  int i,sz;
  
  cout<<"Showing data...\n";
  sz = data.size();
  for(i=0;i<sz;i++){

    cout<<"i= "<<i<<" "<<data[i].elt1<<" "<<data[i].orb_type1<<" "<<data[i].elt2<<" "<<data[i].orb_type2<<"\n";   
    cout<<"K_values      C_values\n";
    for(int k=0;k<5;k++){
      cout<<data[i].K_value[k]<<"  "<<data[i].C_value[k]<<endl; 
    }
    cout<<endl;

  }

//  vector<pp_data_element> pp_data;  
//  vector<psps_data_element> psps_data;

}








void mEHT_K::set_mapping(EHT_K& eht_k, const vector<AO>& basis_ao){
/**
  This function maps the data of EHT_K type (which is used as an intermediate storage - when we read the parameters from the
  external files) onto the interanal variables of the mEHT_K type (which is much more efficient and, therefore, used in actual
  computations). The internal storage of the mEHT_K type is specific to given system (as represented by the AO basis provided in
  the input).

  \param[in] eht_k The object containing the EHT parameters for generic orbital pairs, but without knowing actual orbital pairs of the
                   system.
  \param[in] basis_ao The list of AOs of the particular system - the parameters for all pairs of these AOs will be precomputed and
                   stored internally in the mEHT_K object
*/
  int k;
  cout<<"In mEHT_K::set_mapping...\n";

  size = basis_ao.size();

  for(k=0;k<5;k++){
    eht_K[k].reserve(size*size);   eht_K[k].resize(size*size,  0.0);
    eht_C[k].reserve(size*size);   eht_C[k].resize(size*size,  0.0);
  }

  //========================= Part 1 =======================================
  map<pair<string,string>,int> at_types;
  map<pair<string,string>,int>::iterator it_type,it_type2;

  int ind = 0;
  for(int I=0;I<size;I++){ // create a map of distinct types of orbitals

    std::pair<std::string,std::string> pr(basis_ao[I].element,basis_ao[I].ao_shell);

    it_type = at_types.find(pr);

    // The pair pr is not found in the at_types (should actually be orb_types) array 
    if(it_type==at_types.end()){  // so, this combination of element and AO shell is new type
      at_types.insert(std::pair<pair<string,string>,int>(pr,ind)); ind++;  // so add it to the array
    }

  }// for I

  cout<<"Print distinct orbital types:\n";  
  cout<<"Atom_element  AO_shell  Index of the orbital\n";
  for(it_type=at_types.begin();it_type!=at_types.end();it_type++){
    cout<<it_type->first.first<<"  "<<it_type->first.second<<"  "<<it_type->second<<endl;
  }



  //========================= Part 2 =======================================
  int ntyp = at_types.size(); // number of distinct orbital types
  vector<vector<vector<double> > > k_vals(5,vector<vector<double> >(ntyp,vector<double>(ntyp,0.0)) );
  vector<vector<vector<double> > > c_vals(5,vector<vector<double> >(ntyp,vector<double>(ntyp,0.0)) );

  for(it_type=at_types.begin();it_type!=at_types.end();it_type++){

    std::string at1 = it_type->first.first;
    std::string sh1 = it_type->first.second;
    int i1 = it_type->second;   // index of this AO type

    for(it_type2=at_types.begin();it_type2!=at_types.end();it_type2++){

      std::string at2 = it_type2->first.first;
      std::string sh2 = it_type2->first.second;
      int i2 = it_type2->second; // index of this AO type

      for(k=0;k<5;k++){
        k_vals[k][i1][i2]  = eht_k.get_K_value(k,at1,sh1,at2,sh2);
        c_vals[k][i1][i2]  = eht_k.get_C_value(k,at1,sh1,at2,sh2);
      }

      //cout<<"i1= "<<i1<<" i2= "<<i2<<" at1= "<<at1<<" sh1= "<<sh1<<" at2= "<<at2<<" sh2= "<<sh2<<" k_vals[0][i1][i2]= "<<k_vals[0][i1][i2]<<endl;

    }// it_type2
  }// it_type


  //========================= Part 3 =======================================

  for(int I=0;I<size;I++){

    int i1 = at_types[std::pair<std::string,std::string>(basis_ao[I].element,basis_ao[I].ao_shell)];

    for(int J=0;J<size;J++){

      // New procedure
      int i2 = at_types[std::pair<std::string,std::string>(basis_ao[J].element,basis_ao[J].ao_shell)];

      //cout<<"I= "<<I<<" i1= "<<i1<<" J= "<<J<<" i2= "<<i2<<" k_vals[i1][i2]= "<<k_vals[0][i1][i2]<<endl;

      for(k=0;k<5;k++){
        eht_K[k][I*size+J]  = k_vals[k][i1][i2];
        eht_C[k][I*size+J]  = c_vals[k][i1][i2];
      }
    }// for j
  }// for i

}

//void 


void mEHT_K::set_mapping1(EHT_K& k, int nat, vector<std::string>& mol_at_types){
/**
  This is similar to set_mapping, but we now assume atomic (not orbital) arrays
*/

  //========================= Part 1 =======================================
  map<std::string,int> at_types;
  map<std::string,int>::iterator it_type,it_type2;

  int ind = 0;
  for(int n=0;n<nat;n++){ // create a map of distinct types of atoms
    
    it_type = at_types.find(mol_at_types[n]);

    if(it_type==at_types.end()){  // new type
      at_types.insert(std::pair<std::string,int>(mol_at_types[n],ind)); ind++;
    }

  }// for nat

 
  // Print distinct types
  int ntyp = at_types.size(); // number of distinct types
  cout<<"# of distinct atomic-resolved types = "<<ntyp<<endl;
  for(it_type=at_types.begin();it_type!=at_types.end();it_type++){

    cout<<it_type->first<<"  "<<it_type->second<<endl;
  }



  //========================= Part 2 =======================================

  vector<vector<std::string> > orb_types(ntyp);
  vector<vector<double> > ppa_vals(ntyp);
  vector<vector<double> > pp0_vals(ntyp);
  vector<vector<double> > pp1_vals(ntyp);
  vector<vector<double> > pp2_vals(ntyp);



  int cnt = 0;
  for(it_type=at_types.begin();it_type!=at_types.end();it_type++){

    std::string at1 = it_type->first;

    for(int i=0;i<k.pp_data.size();i++){

      if(k.pp_data[i].elt1==at1){

        orb_types[cnt].push_back(k.pp_data[i].orb_type1);
        ppa_vals[cnt].push_back(k.pp_data[i].PPa_value);
        pp0_vals[cnt].push_back(k.pp_data[i].PP0_value);
        pp1_vals[cnt].push_back(k.pp_data[i].PP1_value);
        pp2_vals[cnt].push_back(k.pp_data[i].PP2_value);

      }// if elt1==at1

    }// for i

    cout<<"Atoms of "<<cnt<<" type at the atoms of "<<at1<<endl;
    cout<<"--- the pseudopotential consists of "<<orb_types[cnt].size()<<" orbital components"<<endl;
    for(int i=0;i<orb_types[cnt].size();i++){
    cout<<"--- --- orbital component "<<i<<" is of type "<<orb_types[cnt][i]<<" with parameters: alp(a.u.)= "
                                      <<ppa_vals[cnt][i]<<" PP0(a.u.)= "
                                      <<pp0_vals[cnt][i]<<" PP1(a.u.)= "
                                      <<pp1_vals[cnt][i]<<" PP2(a.u.)= "
                                      <<pp2_vals[cnt][i]<<endl;

    }// for i
    

    cnt++;

  }// for it_type



  //========================= Part 3 =======================================
  eht_PPa.reserve(nat);    eht_PPa.resize(nat);
  eht_PP0.reserve(nat);    eht_PP0.resize(nat);
  eht_PP1.reserve(nat);    eht_PP1.resize(nat);
  eht_PP2.reserve(nat);    eht_PP2.resize(nat);


  for(int n=0;n<nat;n++){
     
    int i1 = at_types[mol_at_types[n]]; // index of given atomic type within the set of the atom types present in this system
    
    eht_PPa[n] = ppa_vals[i1];
    eht_PP0[n] = pp0_vals[i1];
    eht_PP1[n] = pp1_vals[i1];
    eht_PP2[n] = pp2_vals[i1];

//    cout<<"Atom # "<<n<<" is of type "<<mol.at_type[n]<<" internal classification index "<<i1<<endl;
//    cout<<"--- the pseudopotential consists of "<<orb_types[cnt].size()<<" orbital components"<<endl;

  }// for n



}// set_mapping1



void set_default_elements(map<std::string, pElement>& PT){

//----------------------------------------------------------------
// Initialize fundamental properties of the atoms - periodic table

  pElement elt;

  elt._set("H",  1);      PT["H"]   = elt;      PT["H"].set_mass(1.008);  
  elt._set("He", 2);      PT["He"]  = elt;      PT["He"].set_mass(4.0026);

  elt._set("Li", 3);      PT["Li"]  = elt;      PT["Li"].set_mass(6.94);  
  elt._set("Be", 4);      PT["Be"]  = elt;      PT["Be"].set_mass(9.0122);
  elt._set("B",  5);      PT["B"]   = elt;      PT["B"].set_mass(10.81);  
  elt._set("C",  6);      PT["C"]   = elt;      PT["C"].set_mass(12.011); 
  elt._set("N",  7);      PT["N"]   = elt;      PT["N"].set_mass(14.007); 
  elt._set("O",  8);      PT["O"]   = elt;      PT["O"].set_mass(15.999); 
  elt._set("F",  9);      PT["F"]   = elt;      PT["F"].set_mass(18.998); 
  elt._set("Ne",10);      PT["Ne"]  = elt;      PT["Ne"].set_mass(20.180);


  elt._set("Na",11);      PT["Na"]  = elt;      PT["Na"].set_mass(22.990);  
  elt._set("Mg",12);      PT["Mg"]  = elt;      PT["Mg"].set_mass(24.305);
  elt._set("Al",13);      PT["Al"]  = elt;      PT["Al"].set_mass(26.982);  
  elt._set("Si",14);      PT["Si"]  = elt;      PT["Si"].set_mass(28.085); 
  elt._set("P", 15);      PT["P"]   = elt;      PT["P"].set_mass(30.974); 
  elt._set("S", 16);      PT["S"]   = elt;      PT["S"].set_mass(32.060); 
  elt._set("Cl",17);      PT["Cl"]  = elt;      PT["Cl"].set_mass(35.450); 
  elt._set("Ar",18);      PT["Ar"]  = elt;      PT["Ar"].set_mass(39.948);


  elt._set("K", 19);      PT["K"]   = elt;      PT["K"].set_mass(39.098);  
  elt._set("Ca",20);      PT["Ca"]  = elt;      PT["Ca"].set_mass(40.078);
  elt._set("Sc",21);      PT["Sc"]  = elt;      PT["Sc"].set_mass(44.956);
  elt._set("Ti",22);      PT["Ti"]  = elt;      PT["Ti"].set_mass(47.867);
  elt._set("V" ,23);      PT["V"]   = elt;      PT["V"].set_mass(50.942);
  elt._set("Cr",24);      PT["Cr"]  = elt;      PT["Cr"].set_mass(51.996);
  elt._set("Mn",25);      PT["Mn"]  = elt;      PT["Mn"].set_mass(54.938);
  elt._set("Fe",26);      PT["Fe"]  = elt;      PT["Fe"].set_mass(55.845);
  elt._set("Co",27);      PT["Co"]  = elt;      PT["Co"].set_mass(58.933);
  elt._set("Ni",28);      PT["Ni"]  = elt;      PT["Ni"].set_mass(58.693);
  elt._set("Cu",29);      PT["Cu"]  = elt;      PT["Cu"].set_mass(63.546);
  elt._set("Zn",30);      PT["Zn"]  = elt;      PT["Zn"].set_mass(65.380);
  elt._set("Ga",31);      PT["Ga"]  = elt;      PT["Ga"].set_mass(69.723);  
  elt._set("Ge",32);      PT["Ge"]  = elt;      PT["Ge"].set_mass(72.630); 
  elt._set("As",33);      PT["As"]  = elt;      PT["As"].set_mass(74.922); 
  elt._set("Se",34);      PT["Se"]  = elt;      PT["Se"].set_mass(78.960); 
  elt._set("Br",35);      PT["Br"]  = elt;      PT["Br"].set_mass(79.904); 
  elt._set("Kr",36);      PT["Kr"]  = elt;      PT["Kr"].set_mass(83.798);

  elt._set("Ru",44);      PT["Ru"]  = elt;      PT["Ru"].set_mass(101.07);

  elt._set("Ta",73);      PT["Ta"]  = elt;      PT["Ta"].set_mass(180.95);


}


void Model_Parameters::set_PT_mapping(const vector<AO>& basis_ao){
// map PT to mPT
  int size = basis_ao.size();
 
  orb_params.reserve(size);  orb_params.resize(size);

  for(int i=0;i<size;i++){

    std::string elt = basis_ao[i].element;
    std::string sh  = basis_ao[i].ao_shell;

    pElement Elt = PT[elt];
    

    OrbParams op;

    op.IP = Elt.IP[sh]; 
    op.EA = Elt.EA[sh];
    op.Nquant = Elt.Nquant[sh];
    op.Nzeta = Elt.Nzeta[sh];
    op.zetas = Elt.zetas[sh];
    op.coeffs = Elt.coeffs[sh];

    op.J_param1 = Elt.J_param1[sh];
    op.J_param2 = Elt.J_param2[sh];
    op.J_param3 = Elt.J_param3[sh];
    op.J_param4 = Elt.J_param4[sh];

    op.G1 = Elt.G1[sh];
    op.F2 = Elt.F2[sh];
    op.beta0 = Elt.beta0[sh];


    orb_params[i] = op;

  }// for i


}


void set_parameters_eht_mapping(Model_Parameters& modprms,const vector<AO>& basis_ao){

  modprms.meht_k.set_mapping(modprms.eht_k, basis_ao);

  modprms.set_PT_mapping(basis_ao);


}

void set_parameters_eht_mapping1(Model_Parameters& modprms, int nat, vector<std::string>& mol_at_types){

  modprms.meht_k.set_mapping1(modprms.eht_k, nat, mol_at_types);
}

}// namespace libmodel_parameters
}// namespace liblibra




