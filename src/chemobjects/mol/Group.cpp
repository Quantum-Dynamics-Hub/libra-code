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

#include "Group.h"

/// liblibra namespace
namespace liblibra{

using namespace libio;

namespace libchemobjects{
namespace libmol{


void Group::init_variables(){

  is_Group_name    = 0;
  is_Group_id      = 0;
//  is_Group_mass    = 0;
  is_Group_radius  = 0;
  is_Group_RB      = 0;
  is_Group_ff_type = 0;

//  is_Group_equilibrium_geometry = 0;
//  is_Group_equilibrium_geometry1= 0;
//  is_Group_equilibrium_geometry2= 0;
//  is_Group_current_geometry     = 0;
//  is_Group_force_constant       = 0;
//  is_Group_force_constant1      = 0;
//  is_Group_force_constant2      = 0;
//  is_Group_force_constant_half  = 0;

  is_Group_bond_order           = 0;
  is_Group_bond_alpha           = 0;
//  is_Group_angle_Cijk           = 0;
//  is_Group_angle_C0             = 0;
//  is_Group_angle_C1             = 0;
//  is_Group_angle_C2             = 0;
//  is_Group_angle_coordination   = 0;
//  is_Group_dihedral_multiplicity= 0;
//  is_Group_dihedral_periodicity = 0;
//  is_Group_pair_A_vdw           = 0;
//  is_Group_pair_B_vdw           = 0;
//  is_Group_pair_sigma_vdw       = 0;
//  is_Group_pair_sigma7_vdw      = 0;
//  is_Group_pair_epsil_vdw       = 0;
//  is_Group_pair_alpha_vdw       = 0;
//  is_Group_pair_q1q2            = 0;
//  is_Group_pair_shift_elec      = 0;
//  is_Group_pair_dielectric      = 0;
//  is_Group_pair_wij             = 0;
//  is_Group_pair_wij_1           = 0;
//  is_Group_pair_wij_2           = 0;
//  is_Group_pair_alpij           = 0;
//  is_Group_pair_alpij_1         = 0;
//  is_Group_pair_alpij_2         = 0;
//  is_Group_pair_rij_1           = 0;
//  is_Group_pair_rij_2           = 0;

//  Group_pair_excluded = 1;        is_Group_pair_excluded = 1;

}

void Group::copy_content(const Group& at){

  // Basic information is copied without respect to if it has been defined
  // in destination object
  globGroup_Size      = at.globGroup_Size;
  locGroup_Size       = at.locGroup_Size;
  Group_Size          = at.Group_Size;
  globGroup_Index     = at.globGroup_Index;
  locGroup_Index      = at.locGroup_Index;
  globMolecule_Index  = at.globMolecule_Index;

  globAtom_Index = at.globAtom_Index;
  locAtom_Index  = at.locAtom_Index;

  // This operator copies data from the source if and only if the corresponding data
  // is defined in the source

  if(at.is_Group_name) { Group_name = at.Group_name;  is_Group_name = 1;}
  if(at.is_Group_id)     { Group_id   = at.Group_id;    is_Group_id   = 1;}
//  if(at.is_Group_mass) { Group_mass = at.Group_mass;  is_Group_mass = 1;}
  if(at.is_Group_radius){ Group_radius = at.Group_radius; is_Group_radius = 1;}
  if(at.is_Group_RB) { Group_RB = at.Group_RB; is_Group_RB = 1;}
  if(at.is_Group_ff_type) { Group_ff_type = at.Group_ff_type; is_Group_ff_type = 1;}

//  if(at.is_Group_equilibrium_geometry){ Group_equilibrium_geometry = at.Group_equilibrium_geometry; is_Group_equilibrium_geometry = 1;}
//  if(at.is_Group_equilibrium_geometry1){Group_equilibrium_geometry1= at.Group_equilibrium_geometry1;is_Group_equilibrium_geometry1= 1;}
//  if(at.is_Group_equilibrium_geometry2){Group_equilibrium_geometry2= at.Group_equilibrium_geometry2;is_Group_equilibrium_geometry2= 1;}

//  if(at.is_Group_current_geometry){ Group_current_geometry   = at.Group_current_geometry;    is_Group_current_geometry   = 1;}
//  if(at.is_Group_force_constant){ Group_force_constant  = at.Group_force_constant;  is_Group_force_constant  = 1;}
//  if(at.is_Group_force_constant1){Group_force_constant1 = at.Group_force_constant1; is_Group_force_constant1 = 1;}
//  if(at.is_Group_force_constant2){Group_force_constant2 = at.Group_force_constant2; is_Group_force_constant2 = 1;}

//  if(at.is_Group_force_constant_half){ Group_force_constant_half   = at.Group_force_constant_half;    is_Group_force_constant_half   = 1;}
  if(at.is_Group_bond_order)     { Group_bond_order   = at.Group_bond_order;    is_Group_bond_order   = 1;}
  if(at.is_Group_bond_alpha)     { Group_bond_alpha   = at.Group_bond_alpha;    is_Group_bond_alpha   = 1;}

//  if(at.is_Group_angle_Cijk)     { Group_angle_Cijk   = at.Group_angle_Cijk;    is_Group_angle_Cijk   = 1;}
//  if(at.is_Group_angle_C0)     { Group_angle_C0   = at.Group_angle_C0;    is_Group_angle_C0   = 1;}
//  if(at.is_Group_angle_C1)     { Group_angle_C1   = at.Group_angle_C1;    is_Group_angle_C1   = 1;}
//  if(at.is_Group_angle_C2)     { Group_angle_C2   = at.Group_angle_C2;    is_Group_angle_C2   = 1;}
//  if(at.is_Group_angle_coordination){ Group_angle_coordination   = at.Group_angle_coordination;    is_Group_angle_coordination   = 1;}
//  if(at.is_Group_dihedral_multiplicity){ Group_dihedral_multiplicity   = at.Group_dihedral_multiplicity;    is_Group_dihedral_multiplicity   = 1;}
//  if(at.is_Group_dihedral_periodicity){ Group_dihedral_periodicity   = at.Group_dihedral_periodicity;    is_Group_dihedral_periodicity   = 1;}
//  if(at.is_Group_pair_A_vdw){ Group_pair_A_vdw   = at.Group_pair_A_vdw;    is_Group_pair_A_vdw   = 1;}
//  if(at.is_Group_pair_B_vdw){ Group_pair_B_vdw   = at.Group_pair_B_vdw;    is_Group_pair_B_vdw   = 1;}
//  if(at.is_Group_pair_sigma_vdw){ Group_pair_sigma_vdw   = at.Group_pair_sigma_vdw;    is_Group_pair_sigma_vdw   = 1;}
//  if(at.is_Group_pair_sigma7_vdw){ Group_pair_sigma7_vdw   = at.Group_pair_sigma7_vdw;    is_Group_pair_sigma7_vdw   = 1;}
//  if(at.is_Group_pair_epsil_vdw){ Group_pair_epsil_vdw   = at.Group_pair_epsil_vdw;    is_Group_pair_epsil_vdw   = 1;}
//  if(at.is_Group_pair_alpha_vdw){ Group_pair_alpha_vdw   = at.Group_pair_alpha_vdw;    is_Group_pair_alpha_vdw   = 1;}

//  if(at.is_Group_pair_q1q2)     { Group_pair_q1q2   = at.Group_pair_q1q2;    is_Group_pair_q1q2   = 1;}
//  if(at.is_Group_pair_shift_elec)  { Group_pair_shift_elec   = at.Group_pair_shift_elec;    is_Group_pair_shift_elec   = 1;}
//  if(at.is_Group_pair_dielectric)  { Group_pair_dielectric   = at.Group_pair_dielectric;    is_Group_pair_dielectric   = 1;}

//  if(at.is_Group_pair_wij)      { Group_pair_wij     = at.Group_pair_wij;      is_Group_pair_wij     = 1;}
//  if(at.is_Group_pair_wij_1)    { Group_pair_wij_1   = at.Group_pair_wij_1;    is_Group_pair_wij_1   = 1;}
//  if(at.is_Group_pair_wij_2)    { Group_pair_wij_2   = at.Group_pair_wij_2;    is_Group_pair_wij_2   = 1;}
//  if(at.is_Group_pair_alpij)    { Group_pair_alpij   = at.Group_pair_alpij;    is_Group_pair_alpij   = 1;}
//  if(at.is_Group_pair_alpij_1)  { Group_pair_alpij_1 = at.Group_pair_alpij_1;  is_Group_pair_alpij_1   = 1;}
//  if(at.is_Group_pair_alpij_2)  { Group_pair_alpij_2 = at.Group_pair_alpij_2;  is_Group_pair_alpij_2   = 1;}
//  if(at.is_Group_pair_rij_1)    { Group_pair_rij_1   = at.Group_pair_rij_1;    is_Group_pair_rij_1   = 1;}
//  if(at.is_Group_pair_rij_2)    { Group_pair_rij_2   = at.Group_pair_rij_2;    is_Group_pair_rij_2   = 1;}

//  if(at.is_Group_pair_excluded)  { Group_pair_excluded = at.Group_pair_excluded;  is_Group_pair_excluded = 1; }

}

Group::Group(){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
}

Group::Group(const Group& gr){
  /********************
    Copy constructor
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of gr object which is defined
  copy_content(gr);
}

Group& Group::operator=(const Group& gr){
  /********************
    Assignment operator
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of gr object which is defined
  copy_content(gr);
  return *this;
}

void Group::set(object at){

  set_value(is_Group_name,   Group_name,   at,"Group_name");
  set_value(is_Group_id,     Group_id,     at,"Group_id");
//  set_value(is_Group_mass,   Group_mass,   at,"Group_mass");
  set_value(is_Group_ff_type, Group_ff_type,     at,"Group_ff_type");

//  set_value(is_Group_equilibrium_geometry, Group_equilibrium_geometry, at,"Group_equilibrium_geometry");
//  set_value(is_Group_equilibrium_geometry1,Group_equilibrium_geometry1,at,"Group_equilibrium_geometry1");
//  set_value(is_Group_equilibrium_geometry2,Group_equilibrium_geometry2,at,"Group_equilibrium_geometry2");

//  set_value(is_Group_current_geometry,     Group_current_geometry,     at,"Group_current_geometry");
//  set_value(is_Group_force_constant,       Group_force_constant,       at,"Group_force_constant");
//  set_value(is_Group_force_constant1,      Group_force_constant1,      at,"Group_force_constant1");
//  set_value(is_Group_force_constant2,      Group_force_constant2,      at,"Group_force_constant2");

//  set_value(is_Group_force_constant_half,  Group_force_constant_half,  at,"Group_force_constant_half");

  set_value(is_Group_bond_order,           Group_bond_order,           at,"Group_bond_order");

//  set_value(is_Group_angle_Cijk,           Group_angle_Cijk,           at,"Group_angle_Cijk");
//  set_value(is_Group_angle_C0,             Group_angle_C0,             at,"Group_angle_C0");
//  set_value(is_Group_angle_C1,             Group_angle_C1,             at,"Group_angle_C1");
//  set_value(is_Group_angle_C2,             Group_angle_C2,             at,"Group_angle_C2");
//  set_value(is_Group_angle_coordination,   Group_angle_coordination,   at,"Group_angle_coordination");

//  set_value(is_Group_dihedral_multiplicity,  Group_dihedral_multiplicity,  at,"Group_dihedral_multiplicity");
//  set_value(is_Group_dihedral_periodicity,   Group_dihedral_periodicity,   at,"Group_dihedral_periodicity");

//  set_value(is_Group_pair_A_vdw,           Group_pair_A_vdw,           at,"Group_pair_A_vdw");
//  set_value(is_Group_pair_B_vdw,           Group_pair_B_vdw,           at,"Group_pair_B_vdw");
//  set_value(is_Group_pair_sigma_vdw,       Group_pair_sigma_vdw,       at,"Group_pair_sigma_vdw");
//  set_value(is_Group_pair_epsil_vdw,       Group_pair_epsil_vdw,       at,"Group_pair_epsil_vdw");
//  set_value(is_Group_pair_q1q2,            Group_pair_q1q2,            at,"Group_pair_q1q2");
//  set_value(is_Group_pair_dielectric,      Group_pair_dielectric,      at,"Group_pair_dielectric");
//  set_value(is_Group_pair_excluded,        Group_pair_excluded,        at,"Group_pair_excluded");

//  set_value(is_Group_pair_wij,             Group_pair_wij,             at,"Group_pair_wij");
//  set_value(is_Group_pair_wij_1,           Group_pair_wij_1,           at,"Group_pair_wij_1");
//  set_value(is_Group_pair_wij_2,           Group_pair_wij_2,           at,"Group_pair_wij_2");
//  set_value(is_Group_pair_alpij,           Group_pair_alpij,           at,"Group_pair_alpij");
//  set_value(is_Group_pair_alpij_1,         Group_pair_alpij_1,         at,"Group_pair_alpij_1");
//  set_value(is_Group_pair_alpij_2,         Group_pair_alpij_2,         at,"Group_pair_alpij_2");
//  set_value(is_Group_pair_rij_1,           Group_pair_rij_1,           at,"Group_pair_rij_1");
//  set_value(is_Group_pair_rij_2,           Group_pair_rij_2,           at,"Group_pair_rij_2");

}

void Group::show_info(){
  std::cout<<"Group "<<globGroup_Index<<" properties: "<<std::endl;
  std::cout<<"Group atom indexes (global) : ";
  for(int i=0;i<globAtom_Index.size()/*Group_Size*/;i++){
    std::cout<<globAtom_Index[i]<<"  ";
  }
  std::cout<<std::endl;

  if(is_Group_name)   {std::cout<<"Group_name = "<<Group_name<<std::endl;   } 
  if(is_Group_id)     {std::cout<<"Group_id = "<<Group_id<<std::endl;     }
//  if(is_Group_mass)   {std::cout<<"Group_mass = "<<Group_mass<<std::endl;   }
  if(is_Group_RB)     { Group_RB.show_info(); }
  if(is_Group_ff_type)   {std::cout<<"Group_ff_type = "<<Group_ff_type<<std::endl;     }

//  if(is_Group_equilibrium_geometry)   {std::cout<<"Group_equilibrium_geometry  = "<<Group_equilibrium_geometry<<std::endl;     }
//  if(is_Group_equilibrium_geometry1)  {std::cout<<"Group_equilibrium_geometry1 = "<<Group_equilibrium_geometry1<<std::endl;     }
//  if(is_Group_equilibrium_geometry2)  {std::cout<<"Group_equilibrium_geometry2 = "<<Group_equilibrium_geometry2<<std::endl;     }
//  if(is_Group_current_geometry)       {std::cout<<"Group_current_geometry = "<<Group_current_geometry<<std::endl;     }
//  if(is_Group_force_constant)         {std::cout<<"Group_force_constant  = "<<Group_force_constant<<std::endl;     }
//  if(is_Group_force_constant1)        {std::cout<<"Group_force_constant1 = "<<Group_force_constant1<<std::endl;     }
//  if(is_Group_force_constant2)        {std::cout<<"Group_force_constant2 = "<<Group_force_constant2<<std::endl;     }
//  if(is_Group_force_constant_half)    {std::cout<<"Group_force_constant_half = "<<Group_force_constant_half<<std::endl;     }
  if(is_Group_bond_order)             {std::cout<<"Group_bond_order = "<<Group_bond_order<<std::endl;     }
//  if(is_Group_angle_Cijk)             {std::cout<<"Group_angle_Cijk = "<<Group_angle_Cijk<<std::endl;     }
//  if(is_Group_angle_C0)               {std::cout<<"Group_angle_C0 = "<<Group_angle_C0<<std::endl;     }
//  if(is_Group_angle_C1)               {std::cout<<"Group_angle_C1 = "<<Group_angle_C1<<std::endl;     }
//  if(is_Group_angle_C2)               {std::cout<<"Group_angle_C2 = "<<Group_angle_C2<<std::endl;     }
//  if(is_Group_angle_coordination)     {std::cout<<"Group_angle_coordination = "<<Group_angle_coordination<<std::endl;     }
//  if(is_Group_dihedral_multiplicity)  {std::cout<<"Group_dihedral_multiplicity = "<<Group_dihedral_multiplicity<<std::endl;     }
//  if(is_Group_dihedral_periodicity)   {std::cout<<"Group_dihedral_periodicity = "<<Group_dihedral_periodicity<<std::endl;     }
//  if(is_Group_pair_A_vdw)             {std::cout<<"Group_pair_A_vdw = "<<Group_pair_A_vdw<<std::endl;     }
//  if(is_Group_pair_B_vdw)             {std::cout<<"Group_pair_B_vdw = "<<Group_pair_B_vdw<<std::endl;     }
//  if(is_Group_pair_sigma_vdw)         {std::cout<<"Group_pair_sigma_vdw = "<<Group_pair_sigma_vdw<<std::endl;     }
//  if(is_Group_pair_sigma7_vdw)        {std::cout<<"Group_pair_sigma7_vdw = "<<Group_pair_sigma7_vdw<<std::endl;     }
//  if(is_Group_pair_epsil_vdw)         {std::cout<<"Group_pair_epsil_vdw = "<<Group_pair_epsil_vdw<<std::endl;     }
//  if(is_Group_pair_alpha_vdw)         {std::cout<<"Group_pair_alpha_vdw = "<<Group_pair_alpha_vdw<<std::endl;     }
//  if(is_Group_pair_q1q2)              {std::cout<<"Group_pair_q1q2 = "<<Group_pair_q1q2<<std::endl;     }
//  if(is_Group_pair_shift_elec)        {std::cout<<"Group_pair_shift_elec = "<<Group_pair_shift_elec<<std::endl;     }
//  if(is_Group_pair_dielectric)        {std::cout<<"Group_pair_dielectric = "<<Group_pair_dielectric<<std::endl;     }
//  if(is_Group_pair_wij)               {std::cout<<"Group_pair_wij = "<<Group_pair_wij<<std::endl;     }
//  if(is_Group_pair_wij_1)             {std::cout<<"Group_pair_wij_1 = "<<Group_pair_wij_1<<std::endl;     }
//  if(is_Group_pair_wij_2)             {std::cout<<"Group_pair_wij_2 = "<<Group_pair_wij_2<<std::endl;     }
//  if(is_Group_pair_alpij)             {std::cout<<"Group_pair_alpij = "<<Group_pair_alpij<<std::endl;     }
//  if(is_Group_pair_alpij_1)           {std::cout<<"Group_pair_alpij_1 = "<<Group_pair_alpij_1<<std::endl;     }
//  if(is_Group_pair_alpij_2)           {std::cout<<"Group_pair_alpij_2 = "<<Group_pair_alpij_2<<std::endl;     }
//  if(is_Group_pair_rij_1)             {std::cout<<"Group_pair_rij_1 = "<<Group_pair_rij_1<<std::endl;     }
//  if(is_Group_pair_rij_2)             {std::cout<<"Group_pair_rij_2 = "<<Group_pair_rij_2<<std::endl;     }
//  if(is_Group_pair_excluded)          {std::cout<<"Group_pair_excluded = "<<Group_pair_excluded<<std::endl; } 

  std::cout<<std::endl;

}

void Group::save(boost::property_tree::ptree& pt,std::string path){

  libio::save(pt,path+".globGroup_Size",globGroup_Size);
  libio::save(pt,path+".locGroup_Size",locGroup_Size);
  libio::save(pt,path+".Group_Size",Group_Size);
  libio::save(pt,path+".globAtom_Index",globAtom_Index);
  libio::save(pt,path+".locAtom_Index",locAtom_Index);
  libio::save(pt,path+".globGroup_Index",globGroup_Index);
  libio::save(pt,path+".locGroup_Index",locGroup_Index);
  libio::save(pt,path+".globMolecule_Index",globMolecule_Index);

  if(is_Group_name){  libio::save(pt,path+".Group_name",Group_name);    }
  if(is_Group_id){  libio::save(pt,path+".Group_id",Group_id);    }
  if(is_Group_radius){  libio::save(pt,path+".Group_radius",Group_radius);    }
  if(is_Group_RB){  Group_RB.save(pt,path+".Group_RB");    }
  if(is_Group_ff_type){  libio::save(pt,path+".Group_ff_type",Group_ff_type);    }
  if(is_Group_bond_order){  libio::save(pt,path+".Group_bond_order",Group_bond_order);    }
  if(is_Group_bond_alpha){  libio::save(pt,path+".Group_bond_alpha",Group_bond_alpha);    }

}
 
void save(boost::property_tree::ptree& pt,std::string path,vector<Group>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".Group"+rt);
  }
}

void Group::load(boost::property_tree::ptree& pt,std::string path,int& status){
  int st;
  status = 0;

  libio::load(pt,path+".globGroup_Size",globGroup_Size,st); if(st==1) { status=1;}
  libio::load(pt,path+".locGroup_Size",locGroup_Size,st); if(st==1) { status=1;}
  libio::load(pt,path+".Group_Size",Group_Size,st); if(st==1) { status=1;}
  libio::load(pt,path+".globAtom_Index",globAtom_Index,st); if(st==1) { status=1;}
  libio::load(pt,path+".locAtom_Index",locAtom_Index,st); if(st==1) { status=1;}
  libio::load(pt,path+".globGroup_Index",globGroup_Index,st); if(st==1) { status=1;}
  libio::load(pt,path+".locGroup_Index",locGroup_Index,st); if(st==1) { status=1;}
  libio::load(pt,path+".globMolecule_Index",globMolecule_Index,st); if(st==1) { status=1;}

  libio::load(pt,path+".Group_name",Group_name,is_Group_name); if(is_Group_name==1) { status=1;}
  libio::load(pt,path+".Group_id",Group_id,is_Group_id); if(is_Group_id==1) { status=1;}
  libio::load(pt,path+".Group_radius",Group_radius,is_Group_radius); if(is_Group_radius==1) { status=1;}
  Group_RB.load(pt,path+".Group_RB",is_Group_RB); if(is_Group_RB==1) { status=1;}
  libio::load(pt,path+".Group_ff_type",Group_ff_type,is_Group_ff_type); if(is_Group_ff_type==1) { status=1;}
  libio::load(pt,path+".Group_bond_order",Group_bond_order,is_Group_bond_order); if(is_Group_bond_order==1) { status=1;}
  libio::load(pt,path+".Group_bond_alpha",Group_bond_alpha,is_Group_bond_alpha); if(is_Group_bond_alpha==1) { status=1;}

}

void load(boost::property_tree::ptree& pt,std::string path,vector<Group>& vt,int& status){
  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      Group x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}


}// namespace libmol
}// namespace libchemobjects
}// liblibra
