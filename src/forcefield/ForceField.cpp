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

#include "ForceField.h"

/// liblibra namespace
namespace liblibra{

using namespace libio;

namespace libforcefield{



/*
void ForceField::set(object at){

   set_value(is_ForceField_Name,       ForceField_Name,      at,"ForceField_Name");
   set_value(is_sigma_comb_rule,       sigma_comb_rule,      at,"sigma_comb_rule");
   set_value(is_epsilon_comb_rule,     epsilon_comb_rule,    at,"epsilon_comb_rule");

   set_value(is_system_pbc,            system_pbc,           at,"system_pbc");
   set_value(is_pbc_degree_x,          pbc_degree_x,         at,"pbc_degree_x");
   set_value(is_pbc_degree_y,          pbc_degree_y,         at,"pbc_degree_y");
   set_value(is_pbc_degree_z,          pbc_degree_z,         at,"pbc_degree_z");
   set_value(is_reciprocal_degree_x,   reciprocal_degree_x,  at,"reciprocal_degree_x");
   set_value(is_reciprocal_degree_y,   reciprocal_degree_y,  at,"reciprocal_degree_y");
   set_value(is_reciprocal_degree_z,   reciprocal_degree_z,  at,"reciprocal_degree_z");
   set_value(is_R_vdw_on,              R_vdw_on,             at,"R_vdw_on");
   set_value(is_R_vdw_off,             R_vdw_off,            at,"R_vdw_off");
   set_value(is_R_elec_on,             R_elec_on,            at,"R_elec_on");
   set_value(is_R_elec_off,            R_elec_off,           at,"R_elec_off");
   set_value(is_R_vlist,               R_vlist,              at,"R_vlist");
   set_value(is_elec_etha,             R_elec_etha,          at,"R_elec_etha");

   set_value(is_bond_functional,       bond_functional,      at,"bond_functional");
   set_value(is_angle_functional,      angle_functional,     at,"angle_functional");
   set_value(is_dihedral_functional,   dihedral_functional,  at,"dihedral_functional");
   set_value(is_oop_functional,        oop_functional,       at,"oop_functional");
   set_value(is_vdw_functional,        vdw_functional,       at,"vdw_functional");
   set_value(is_elec_functional,       elec_functional,      at,"elec_functional");

   set_value(is_vdw_scale12,           vdw_scale12,          at,"vdw_scale12");
   set_value(is_vdw_scale13,           vdw_scale13,          at,"vdw_scale13");
   set_value(is_vdw_scale14,           vdw_scale14,          at,"vdw_scale14");
   set_value(is_elec_scale12,          elec_scale12,         at,"elec_scale12");
   set_value(is_elec_scale13,          elec_scale13,         at,"elec_scale13");
   set_value(is_elec_scale14,          elec_scale14,         at,"elec_scale14");

}
*/

void ForceField::set(boost::python::dict d){
  extract_dictionary(d);
}

void ForceField::extract_dictionary(boost::python::dict d){
  std::string key;
  for(int i=0;i<len(d.values());i++){
    key = extract<std::string>(d.keys()[i]);

    if(key=="ForceField_Name"){ ForceField_Name = extract<std::string>(d.values()[i]); is_ForceField_Name = 1; }
    else if(key=="bond" || key=="bond_functional"){ bond_functional = extract<std::string>(d.values()[i]); is_bond_functional = 1;}
    else if(key=="angle" || key=="angle_functional"){ angle_functional = extract<std::string>(d.values()[i]); is_angle_functional = 1;}
    else if(key=="dihedral" || key=="dihedral_functional"){ dihedral_functional = extract<std::string>(d.values()[i]); is_dihedral_functional = 1;}
    else if(key=="oop" || key=="oop_functional"){ oop_functional = extract<std::string>(d.values()[i]); is_oop_functional = 1;}
    else if(key=="vdw" || key=="vdw_functional"){ vdw_functional = extract<std::string>(d.values()[i]); is_vdw_functional = 1;}
    else if(key=="elec" ||  key=="elec_functional"){ elec_functional = extract<std::string>(d.values()[i]); is_elec_functional = 1;}
    else if(key=="mb" || key=="mb_functional"){ mb_functional = extract<std::string>(d.values()[i]); is_mb_functional = 1;}
    else if(key=="cg" || key=="cg_functional"){ cg_functional = extract<std::string>(d.values()[i]); is_cg_functional = 1;}
    else if(key=="mb_excl" || key=="mb_excl_functional"){ mb_excl_functional = extract<std::string>(d.values()[i]); is_mb_excl_functional = 1;}

//    else if(key=="stress_opt"){ stress_opt = extract<std::string>(d.values()[i]); is_stress_opt = 1;}

    else if(key=="system_pbc"){ system_pbc = extract<std::string>(d.values()[i]); is_system_pbc = 1; }
    else if(key=="pbc_degree_x"){ pbc_degree_x = extract<int>(d.values()[i]); is_pbc_degree_x = 1; }
    else if(key=="pbc_degree_y"){ pbc_degree_y = extract<int>(d.values()[i]); is_pbc_degree_y = 1; }
    else if(key=="pbc_degree_z"){ pbc_degree_z = extract<int>(d.values()[i]); is_pbc_degree_z = 1; }
    else if(key=="reciprocal_degree_x"){ reciprocal_degree_x = extract<int>(d.values()[i]); is_reciprocal_degree_x = 1; }
    else if(key=="reciprocal_degree_y"){ reciprocal_degree_y = extract<int>(d.values()[i]); is_reciprocal_degree_y = 1; }
    else if(key=="reciprocal_degree_z"){ reciprocal_degree_z = extract<int>(d.values()[i]); is_reciprocal_degree_z = 1; }
    else if(key=="R_vdw_on"){ R_vdw_on = extract<double>(d.values()[i]); is_R_vdw_on = 1; }
    else if(key=="R_vdw_off"){ R_vdw_off = extract<double>(d.values()[i]); is_R_vdw_off = 1; }
    else if(key=="R_elec_on"){ R_elec_on = extract<double>(d.values()[i]); is_R_elec_on = 1; }
    else if(key=="R_elec_off"){ R_elec_off = extract<double>(d.values()[i]); is_R_elec_off = 1; }
    else if(key=="R_vlist"){ R_vlist = extract<double>(d.values()[i]); is_R_vlist = 1; }
    else if(key=="elec_etha"){ elec_etha = extract<double>(d.values()[i]); is_elec_etha = 1; }

    else if(key=="sigma_comb_rule"){ sigma_comb_rule = extract<std::string>(d.values()[i]); is_sigma_comb_rule = 1; }
    else if(key=="epsilon_comb_rule"){ epsilon_comb_rule = extract<std::string>(d.values()[i]); is_epsilon_comb_rule = 1; }
    else if(key=="vdw_scale12"){ vdw_scale12 = extract<double>(d.values()[i]); is_vdw_scale12 = 1; }
    else if(key=="vdw_scale13"){ vdw_scale13 = extract<double>(d.values()[i]); is_vdw_scale13 = 1; }
    else if(key=="vdw_scale14"){ vdw_scale14 = extract<double>(d.values()[i]); is_vdw_scale14 = 1; }
    else if(key=="elec_scale12"){ elec_scale12 = extract<double>(d.values()[i]); is_elec_scale12 = 1; }
    else if(key=="elec_scale13"){ elec_scale13 = extract<double>(d.values()[i]); is_elec_scale13 = 1; }
    else if(key=="elec_scale14"){ elec_scale14 = extract<double>(d.values()[i]); is_elec_scale14 = 1; }


  }// for i

}

void ForceField::init_variables(){
  is_ForceField_Name               = 0;

  is_system_pbc                    = 0;
  is_pbc_degree_x                  = 0;
  is_pbc_degree_y                  = 0;
  is_pbc_degree_z                  = 0;
  is_reciprocal_degree_x           = 0;
  is_reciprocal_degree_y           = 0;
  is_reciprocal_degree_z           = 0;
  is_R_vdw_on                      = 0;
  is_R_vdw_off                     = 0;
  is_R_elec_on                     = 0;
  is_R_elec_off                    = 0;
  is_R_vlist                       = 0;
  is_elec_etha                     = 0;

  is_bond_functional               = 0;
  is_angle_functional              = 0;
  is_dihedral_functional           = 0;
  is_oop_functional                = 0;
  is_vdw_functional                = 0;
  is_elec_functional               = 0;
  is_mb_functional                 = 0;
  is_cg_functional                 = 0;
  is_mb_excl_functional            = 0;


//  stress_opt = "fr";               is_stress_opt = 1;

  is_sigma_comb_rule    = 0;
  is_epsilon_comb_rule  = 0;
  is_vdw_scale12        = 0;
  is_vdw_scale13        = 0;
  is_vdw_scale14        = 0;
  is_elec_scale12       = 0;
  is_elec_scale13       = 0;
  is_elec_scale14       = 0;

}

void ForceField::copy_content(const ForceField& ff){

  if(ff.is_ForceField_Name){ ForceField_Name = ff.ForceField_Name; is_ForceField_Name = 1; }

  if(ff.is_bond_functional){ bond_functional = ff.bond_functional; is_bond_functional = 1; }
  if(ff.is_angle_functional){ angle_functional = ff.angle_functional; is_angle_functional = 1; }
  if(ff.is_dihedral_functional){ dihedral_functional = ff.dihedral_functional; is_dihedral_functional = 1; }
  if(ff.is_oop_functional){ oop_functional = ff.oop_functional; is_oop_functional = 1; }
  if(ff.is_vdw_functional){ vdw_functional = ff.vdw_functional; is_vdw_functional = 1; }
  if(ff.is_elec_functional){ elec_functional = ff.elec_functional; is_elec_functional = 1; }
  if(ff.is_mb_functional){ mb_functional = ff.mb_functional; is_mb_functional = 1; }
  if(ff.is_cg_functional){ cg_functional = ff.cg_functional; is_cg_functional = 1; }
  if(ff.is_mb_excl_functional){ mb_excl_functional = ff.mb_excl_functional; is_mb_excl_functional = 1; }


//  if(ff.is_stress_opt){ stress_opt = ff.stress_opt; is_stress_opt = 1; }

  if(ff.is_system_pbc){ system_pbc = ff.system_pbc; is_system_pbc = 1; }
  if(ff.is_pbc_degree_x){ pbc_degree_x = ff.pbc_degree_x; is_pbc_degree_x = 1; }
  if(ff.is_pbc_degree_y){ pbc_degree_y = ff.pbc_degree_y; is_pbc_degree_y = 1; }
  if(ff.is_pbc_degree_z){ pbc_degree_z = ff.pbc_degree_z; is_pbc_degree_z = 1; }
  if(ff.is_reciprocal_degree_x){ reciprocal_degree_x = ff.reciprocal_degree_x; is_reciprocal_degree_x = 1; }
  if(ff.is_reciprocal_degree_y){ reciprocal_degree_y = ff.reciprocal_degree_y; is_reciprocal_degree_y = 1; }
  if(ff.is_reciprocal_degree_z){ reciprocal_degree_z = ff.reciprocal_degree_z; is_reciprocal_degree_z = 1; }
  if(ff.is_R_vdw_on){ R_vdw_on = ff.R_vdw_on; is_R_vdw_on = 1; }
  if(ff.is_R_vdw_off){ R_vdw_off = ff.R_vdw_off; is_R_vdw_off = 1; }
  if(ff.is_R_elec_on){ R_elec_on = ff.R_elec_on; is_R_elec_on = 1; }
  if(ff.is_R_elec_off){ R_elec_off = ff.R_elec_off; is_R_elec_off = 1; }
  if(ff.is_R_vlist){ R_vlist = ff.R_vlist; is_R_vlist = 1; }
  if(ff.is_elec_etha){ elec_etha = ff.elec_etha; is_elec_etha = 1; }

  if(ff.is_sigma_comb_rule){ sigma_comb_rule = ff.sigma_comb_rule; is_sigma_comb_rule = 1; }
  if(ff.is_epsilon_comb_rule){ epsilon_comb_rule = ff.epsilon_comb_rule; is_epsilon_comb_rule = 1; }
  if(ff.is_vdw_scale12){ vdw_scale12 = ff.vdw_scale12; is_vdw_scale12 = 1; }
  if(ff.is_vdw_scale13){ vdw_scale13 = ff.vdw_scale13; is_vdw_scale13 = 1; }
  if(ff.is_vdw_scale14){ vdw_scale14 = ff.vdw_scale14; is_vdw_scale14 = 1; }
  if(ff.is_elec_scale12){ elec_scale12 = ff.elec_scale12; is_elec_scale12 = 1; }
  if(ff.is_elec_scale13){ elec_scale13 = ff.elec_scale13; is_elec_scale13 = 1; }
  if(ff.is_elec_scale14){ elec_scale14 = ff.elec_scale14; is_elec_scale14 = 1; }

  Atom_Records = ff.Atom_Records;
  Bond_Records = ff.Bond_Records;
  Angle_Records = ff.Angle_Records;
  Dihedral_Records = ff.Dihedral_Records;
  Improper_Records = ff.Improper_Records;
  Fragment_Records = ff.Fragment_Records;
}

ForceField::ForceField(){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
}

ForceField::ForceField(boost::python::dict d){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
  extract_dictionary(d);
}

ForceField::ForceField(const ForceField& ff){
  /********************
    Copy constructor
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(ff);
}

ForceField& ForceField::operator=(const ForceField& ff){
  /********************
    Assignment operator
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(ff);
  return *this;
}

ForceField::~ForceField(){
  if(Atom_Records.size()>0) { Atom_Records.clear(); }
  if(Bond_Records.size()>0) { Bond_Records.clear(); }
  if(Angle_Records.size()>0) { Angle_Records.clear(); }
  if(Dihedral_Records.size()>0) { Dihedral_Records.clear(); }
  if(Improper_Records.size()>0) { Improper_Records.clear(); }
  if(Fragment_Records.size()>0) { Fragment_Records.clear(); }
}

void ForceField::show_info() const {

  std::cout<<"ForceField properties:"<<std::endl;
  if(is_ForceField_Name)      {std::cout<<"ForceField_Name = "<<ForceField_Name<<std::endl;   }
  if(is_sigma_comb_rule)      {std::cout<<"sigma_comb_rule = "<<sigma_comb_rule<<std::endl;   }
  if(is_epsilon_comb_rule)    {std::cout<<"epsilon_comb_rule = "<<epsilon_comb_rule<<std::endl;   }

  if(is_system_pbc)           {std::cout<<"system_pbc = "<<system_pbc<<std::endl; }
  if(is_pbc_degree_x)         {std::cout<<"pbc_degree_x = "<<pbc_degree_x<<std::endl; }
  if(is_pbc_degree_y)         {std::cout<<"pbc_degree_y = "<<pbc_degree_y<<std::endl; }
  if(is_pbc_degree_z)         {std::cout<<"pbc_degree_z = "<<pbc_degree_z<<std::endl; }
  if(is_reciprocal_degree_x)  {std::cout<<"reciprocal_degree_x = "<<reciprocal_degree_x<<std::endl; }
  if(is_reciprocal_degree_y)  {std::cout<<"reciprocal_degree_y = "<<reciprocal_degree_y<<std::endl; }
  if(is_reciprocal_degree_z)  {std::cout<<"reciprocal_degree_z = "<<reciprocal_degree_z<<std::endl; }
  if(is_R_vdw_on)             {std::cout<<"R_vdw_on = "<<R_vdw_on<<std::endl; }
  if(is_R_vdw_off)            {std::cout<<"R_vdw_off = "<<R_vdw_off<<std::endl; }
  if(is_R_elec_on)            {std::cout<<"R_elec_on = "<<R_elec_on<<std::endl; }
  if(is_R_elec_off)           {std::cout<<"R_elec_off = "<<R_elec_off<<std::endl; }
  if(is_R_vlist)              {std::cout<<"R_vlist = "<<R_vlist<<std::endl; }
  if(is_elec_etha)            {std::cout<<"elec_etha = "<<elec_etha<<std::endl; }

  if(is_bond_functional)      {std::cout<<"bond_functional = "<<bond_functional<<std::endl; }
  if(is_angle_functional)     {std::cout<<"angle_functional = "<<angle_functional<<std::endl; }
  if(is_dihedral_functional)  {std::cout<<"dihedral_functional = "<<dihedral_functional<<std::endl; }
  if(is_oop_functional)       {std::cout<<"oop_functional = "<<oop_functional<<std::endl; }
  if(is_vdw_functional)       {std::cout<<"vdw_functional = "<<vdw_functional<<std::endl; }
  if(is_elec_functional)      {std::cout<<"elec_functional = "<<elec_functional<<std::endl; }
  if(is_mb_functional)        {std::cout<<"mb_functional = "<<mb_functional<<std::endl; }
  if(is_cg_functional)        {std::cout<<"cg_functional = "<<cg_functional<<std::endl; }
  if(is_mb_excl_functional)   {std::cout<<"mb_excl_functional = "<<mb_excl_functional<<std::endl; }

//  if(is_stress_opt)           {std::cout<<"stress_opt = "<<stress_opt<<std::endl; }

  if(is_vdw_scale12)          {std::cout<<"vdw_scale12 = "<<vdw_scale12<<std::endl; }
  if(is_vdw_scale13)          {std::cout<<"vdw_scale13 = "<<vdw_scale13<<std::endl; }
  if(is_vdw_scale14)          {std::cout<<"vdw_scale14 = "<<vdw_scale14<<std::endl; }
  if(is_elec_scale12)         {std::cout<<"elec_scale12 = "<<elec_scale12<<std::endl; }
  if(is_elec_scale13)         {std::cout<<"elec_scale13 = "<<elec_scale13<<std::endl; }
  if(is_elec_scale14)         {std::cout<<"elec_scale14 = "<<elec_scale14<<std::endl; }

}


void ForceField::save(boost::property_tree::ptree& pt,std::string path){

  if(is_ForceField_Name){  libio::save(pt,path+".ForceField_Name",ForceField_Name);    }
  if(is_bond_functional){  libio::save(pt,path+".bond_functional",bond_functional);    }
  if(is_angle_functional){  libio::save(pt,path+".angle_functional",angle_functional);    }
  if(is_dihedral_functional){  libio::save(pt,path+".dihedral_functional",dihedral_functional);    }
  if(is_oop_functional){  libio::save(pt,path+".oop_functional",oop_functional);    }
  if(is_vdw_functional){  libio::save(pt,path+".vdw_functional",vdw_functional);    }
  if(is_elec_functional){  libio::save(pt,path+".elec_functional",elec_functional);    }
  if(is_mb_functional){  libio::save(pt,path+".mb_functional",mb_functional);    }
  if(is_cg_functional){  libio::save(pt,path+".cg_functional",cg_functional);    }
  if(is_mb_excl_functional){  libio::save(pt,path+".mb_excl_functional",mb_excl_functional);    }

  if(is_system_pbc){  libio::save(pt,path+".system_pbc",system_pbc);    }
  if(is_pbc_degree_x){  libio::save(pt,path+".pbc_degree_x",pbc_degree_x);    }
  if(is_pbc_degree_y){  libio::save(pt,path+".pbc_degree_y",pbc_degree_y);    }
  if(is_pbc_degree_z){  libio::save(pt,path+".pbc_degree_z",pbc_degree_z);    }
  if(is_reciprocal_degree_x){  libio::save(pt,path+".reciprocal_degree_x",reciprocal_degree_x);    }
  if(is_reciprocal_degree_y){  libio::save(pt,path+".reciprocal_degree_y",reciprocal_degree_y);    }
  if(is_reciprocal_degree_z){  libio::save(pt,path+".reciprocal_degree_z",reciprocal_degree_z);    }
  if(is_R_vdw_on){  libio::save(pt,path+".R_vdw_on",R_vdw_on);    }
  if(is_R_vdw_off){  libio::save(pt,path+".R_vdw_off",R_vdw_off);    }
  if(is_R_elec_on){  libio::save(pt,path+".R_elec_on",R_elec_on);    }
  if(is_R_elec_off){  libio::save(pt,path+".R_elec_off",R_elec_off);    }
  if(is_R_vlist){  libio::save(pt,path+".R_vlist",R_vlist);    }
  if(is_elec_etha){  libio::save(pt,path+".elec_etha",elec_etha);    }

  if(is_sigma_comb_rule){  libio::save(pt,path+".sigma_comb_rule",sigma_comb_rule);    }
  if(is_epsilon_comb_rule){  libio::save(pt,path+".epsilon_comb_rule",epsilon_comb_rule);    }
  if(is_vdw_scale12){  libio::save(pt,path+".vdw_scale12",vdw_scale12);    }
  if(is_vdw_scale13){  libio::save(pt,path+".vdw_scale13",vdw_scale13);    }
  if(is_vdw_scale14){  libio::save(pt,path+".vdw_scale14",vdw_scale14);    }
  if(is_elec_scale12){  libio::save(pt,path+".elec_scale12",elec_scale12);    }
  if(is_elec_scale13){  libio::save(pt,path+".elec_scale13",elec_scale13);    }
  if(is_elec_scale14){  libio::save(pt,path+".elec_scale14",elec_scale14);    }

//  namespace here = libhamiltonian::libhamiltonian_atomistic::libhamiltonian_mm::libforcefield;


  if(Atom_Records.size()>0){  libforcefield::save(pt,path+".Atom_Records",Atom_Records);    }
  if(Bond_Records.size()>0){  libforcefield::save(pt,path+".Bond_Records",Bond_Records);    }
  if(Angle_Records.size()>0){  libforcefield::save(pt,path+".Angle_Records",Angle_Records);    }
  if(Dihedral_Records.size()>0){  libforcefield::save(pt,path+".Dihedral_Records",Dihedral_Records);    }
  if(Improper_Records.size()>0){  libforcefield::save(pt,path+".Improper_Records",Improper_Records);    }
  if(Fragment_Records.size()>0){  libforcefield::save(pt,path+".Fragment_Records",Fragment_Records);    }

}

void ForceField::load(boost::property_tree::ptree& pt,std::string path,int& status){
  int st;
  status = 0;

  libio::load(pt,path+".ForceField_Name",ForceField_Name,is_ForceField_Name); if(is_ForceField_Name==1) { status=1;}
  libio::load(pt,path+".bond_functional",bond_functional,is_bond_functional); if(is_bond_functional==1) { status=1;}
  libio::load(pt,path+".angle_functional",angle_functional,is_angle_functional); if(is_angle_functional==1) { status=1;}
  libio::load(pt,path+".dihedral_functional",dihedral_functional,is_dihedral_functional); if(is_dihedral_functional==1) { status=1;}
  libio::load(pt,path+".oop_functional",oop_functional,is_oop_functional); if(is_oop_functional==1) { status=1;}
  libio::load(pt,path+".vdw_functional",vdw_functional,is_vdw_functional); if(is_vdw_functional==1) { status=1;}
  libio::load(pt,path+".elec_functional",elec_functional,is_elec_functional); if(is_elec_functional==1) { status=1;}
  libio::load(pt,path+".mb_functional",mb_functional,is_mb_functional); if(is_mb_functional==1) { status=1;}
  libio::load(pt,path+".cg_functional",cg_functional,is_cg_functional); if(is_cg_functional==1) { status=1;}
  libio::load(pt,path+".mb_excl_functional",mb_excl_functional,is_mb_excl_functional); if(is_mb_excl_functional==1) { status=1;}


  libio::load(pt,path+".system_pbc",system_pbc,is_system_pbc); if(is_system_pbc==1) { status=1;}
  libio::load(pt,path+".pbc_degree_x",pbc_degree_x,is_pbc_degree_x); if(is_pbc_degree_x==1) { status=1;}
  libio::load(pt,path+".pbc_degree_y",pbc_degree_y,is_pbc_degree_y); if(is_pbc_degree_y==1) { status=1;}
  libio::load(pt,path+".pbc_degree_z",pbc_degree_z,is_pbc_degree_z); if(is_pbc_degree_z==1) { status=1;}
  libio::load(pt,path+".reciprocal_degree_x",reciprocal_degree_x,is_reciprocal_degree_x); if(is_reciprocal_degree_x==1) { status=1;}
  libio::load(pt,path+".reciprocal_degree_y",reciprocal_degree_y,is_reciprocal_degree_y); if(is_reciprocal_degree_y==1) { status=1;}
  libio::load(pt,path+".reciprocal_degree_z",reciprocal_degree_z,is_reciprocal_degree_z); if(is_reciprocal_degree_z==1) { status=1;}
  libio::load(pt,path+".R_vdw_on",R_vdw_on,is_R_vdw_on); if(is_R_vdw_on==1) { status=1;}
  libio::load(pt,path+".R_vdw_off",R_vdw_off,is_R_vdw_off); if(is_R_vdw_off==1) { status=1;}
  libio::load(pt,path+".R_elec_on",R_elec_on,is_R_elec_on); if(is_R_elec_on==1) { status=1;}
  libio::load(pt,path+".R_elec_off",R_elec_off,is_R_elec_off); if(is_R_elec_off==1) { status=1;}
  libio::load(pt,path+".R_vlist",R_vlist,is_R_vlist); if(is_R_vlist==1) { status=1;}
  libio::load(pt,path+".elec_etha",elec_etha,is_elec_etha); if(is_elec_etha==1) { status=1;}
  libio::load(pt,path+".sigma_comb_rule",sigma_comb_rule,is_sigma_comb_rule); if(is_sigma_comb_rule==1) { status=1;}
  libio::load(pt,path+".epsilon_comb_rule",epsilon_comb_rule,is_epsilon_comb_rule); if(is_epsilon_comb_rule==1) { status=1;}
  libio::load(pt,path+".vdw_scale12",vdw_scale12,is_vdw_scale12); if(is_vdw_scale12==1) { status=1;}
  libio::load(pt,path+".vdw_scale13",vdw_scale13,is_vdw_scale13); if(is_vdw_scale13==1) { status=1;}
  libio::load(pt,path+".vdw_scale14",vdw_scale14,is_vdw_scale14); if(is_vdw_scale14==1) { status=1;}
  libio::load(pt,path+".elec_scale12",elec_scale12,is_elec_scale12); if(is_elec_scale12==1) { status=1;}
  libio::load(pt,path+".elec_scale13",elec_scale13,is_elec_scale13); if(is_elec_scale13==1) { status=1;}
  libio::load(pt,path+".elec_scale14",elec_scale14,is_elec_scale14); if(is_elec_scale14==1) { status=1;}

//  namespace here = libhamiltonian::libhamiltonian_atomistic::libhamiltonian_mm::libforcefield;

  libforcefield::load(pt,path+".Atom_Records",Atom_Records,st); if(st==1) { status=1;}
  libforcefield::load(pt,path+".Bond_Records",Bond_Records,st); if(st==1) { status=1;}
  libforcefield::load(pt,path+".Angle_Records",Angle_Records,st); if(st==1) { status=1;}
  libforcefield::load(pt,path+".Dihedral_Records",Dihedral_Records,st); if(st==1) { status=1;}
  libforcefield::load(pt,path+".Improper_Records",Improper_Records,st); if(st==1) { status=1;}
  libforcefield::load(pt,path+".Fragment_Records",Fragment_Records,st); if(st==1) { status=1;}

}




int ForceField::Atom_Record_Index(int Atom_ff_int_type){
/****************************************************************
   This function searches for index of Atom_Records vector in which
   data about atom of force field type "Atom_ff_int_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Atom_Records.size();

   for(int i=0;i<sz;i++){
       if(Atom_Records[i].is_Atom_ff_int_type){
          if(Atom_Records[i].Atom_ff_int_type==Atom_ff_int_type){
             indx = i;
             break;
          }
       }
   }

   return indx;
}

int ForceField::Atom_Record_Index(std::string Atom_ff_type){
/****************************************************************
   This function searches for index of Atom_Records vector in which
   data about atom of force field type "Atom_ff_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Atom_Records.size();

   for(int i=0;i<sz;i++){
       if(Atom_Records[i].is_Atom_ff_type){
          if(Atom_Records[i].Atom_ff_type==Atom_ff_type){
             indx = i;
             break;
          }
       }
   }

   return indx;
}


int ForceField::Atom_Record_Index_by_Element(std::string Atom_element){
/****************************************************************
   This function searches for index of Atom_Records vector in which
   data about atom of element name "Atom_element" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Atom_Records.size();

   for(int i=0;i<sz;i++){
       if(Atom_Records[i].is_Atom_element){
          if(Atom_Records[i].Atom_element==Atom_element){
             indx = i;
             break;
          }
       }
   }

   return indx;
}


int ForceField::Atom_Record_Index_by_Element(std::string Atom_element,vector<int>& res){
/****************************************************************
   This function searches for index of Atom_Records vector in which
   data about atom of element name "Atom_element" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Atom_Records.size();
   if(res.size()>0) { res.clear(); }

   for(int i=0;i<sz;i++){
       if(Atom_Records[i].is_Atom_element){
          if(Atom_Records[i].Atom_element==Atom_element){
             indx = i;
             res.push_back(i);
          }
       }
   }

   return indx;
}

int ForceField::Atom_Record_Index_by_Element(int Atom_atomic_number,vector<int>& res){
/****************************************************************
   This function searches for index of Atom_Records vector in which
   data about atom of element number "Atom_atomic_number" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Atom_Records.size();
   if(res.size()>0) { res.clear(); }

   for(int i=0;i<sz;i++){
       if(Atom_Records[i].is_Atom_atomic_number){
          if(Atom_Records[i].Atom_atomic_number==Atom_atomic_number){
             indx = i;
             res.push_back(i);          
          }
       }
   }

   return indx;
}

int ForceField::Add_Atom_Record(Atom_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Atom_Record rec in array Atom_Records. If this atom type (atom record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the atom record rec is not valid - will not be added.
*********************************************************************/
   int res = 1;
   int sz = Atom_Records.size();

   int int_types = (rec.is_Atom_ff_int_type);
   int sym_types = (rec.is_Atom_ff_type);
   int elt_types = (rec.is_Atom_atomic_number);

   vector<int> at_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_type
   last_indx = Atom_Record_Index(rec.Atom_ff_int_type);
   at_indxs.push_back(last_indx);
   }
   else if(sym_types){
   // Search on the basis of symbolic type
   last_indx = Atom_Record_Index(rec.Atom_ff_type);
   at_indxs.push_back(last_indx);
   }
   else if(elt_types){
   // Search on the basis of atomic numbers
   last_indx = Atom_Record_Index_by_Element(rec.Atom_atomic_number,at_indxs);
   }
   else{
   // Record will not be added to the table of Angle_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2;
   }

   if(last_indx==-2){ // In this case the record was not added to Atom_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The atom record is valid but does not exist in Atom_Records array
                           // so it will be added completely
       Atom_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar atom record is already defined - try to add new properties

       for(int i=0;i<at_indxs.size();i++){
           Atom_Records[at_indxs[i]].merge(rec);
       }
       res = 0;
   }

   return res;

}

/*
int ForceField::Add_Atom_Record(Atom_Record rec){
*********************************************************************
   This is user-interface function. It checks for existance of the 
   Atom_Record rec in array Atom_Records. If this atom type (atom record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array (we did not add it this time). Returns -1
   if there is an error.
*********************************************************************
   int res = 1;   
   int sz = Atom_Records.size();
   int indx = sz;
   Atom_Record record; 
   record = rec;  

   if(rec.is_Atom_ff_type){
   
      for(int i=0;i<sz;i++){
       
         if(Atom_Records[i].is_Atom_ff_type){ 
            // This type already exist
            if(Atom_Records[i].Atom_ff_type==rec.Atom_ff_type){
                //res = 0;  - this is old version: means do not modify existing data
                res = 1; // new - add new data to existing record
                record = Atom_Records[i];
                indx = i; break; 
            }            
         }        
      }// for i

   }// is_Atom_ff_type == 1 (defined)
   else{
    std::cout<<"To add atom type record to the force field its Atom_ff_type should be defined"<<std::endl;
    res = -1;
   }

   if(res==1){
      // We will add this atom type to the force field 
      if(!rec.is_Atom_ff_int_type){
          // If Atom_ff_int_type not defined we define it
          // to be an index of this record in array
          rec.is_Atom_ff_int_type = 1;
          rec.Atom_ff_int_type = indx;
      }
      Atom_Records.push_back(record);
   }
   
   return res;

}
*/



int ForceField::Bond_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type){
/****************************************************************
   This function searches for index of Bond_Records vector in which
   data about bond formed by 2 atoms of of force field types
   "Atom1_ff_int_type" and "Atom2_ff_int_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Bond_Records.size();

   for(int i=0;i<sz;i++){
       if(Bond_Records[i].is_Atom1_ff_int_type && Bond_Records[i].is_Atom2_ff_int_type){
          if(((Bond_Records[i].Atom1_ff_int_type==Atom1_ff_int_type)&&(Bond_Records[i].Atom2_ff_int_type==Atom2_ff_int_type) )|| 
             ((Bond_Records[i].Atom1_ff_int_type==Atom2_ff_int_type)&&(Bond_Records[i].Atom2_ff_int_type==Atom1_ff_int_type) )
            )
          {
             indx = i;
             break;
          }
       }
   }

   return indx;
}


int ForceField::Bond_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type,int order){
/****************************************************************
   This function searches for index of Bond_Records vector in which
   data about bond formed by 2 atoms of of force field types
   "Atom1_ff_int_type" and "Atom2_ff_int_type" are stored
   Returns -1 if such index has not been found
   order parameter distinguish bonds 1-2 and 2-1
*****************************************************************/

   int indx = -1;
   int sz   = Bond_Records.size();

   int cmpr11,cmpr12,cmpr21,cmpr22,res;
   for(int i=0;i<sz;i++){
       if(Bond_Records[i].is_Atom1_ff_int_type && Bond_Records[i].is_Atom2_ff_int_type){
          cmpr11 = (Bond_Records[i].Atom1_ff_int_type==Atom1_ff_int_type);
          cmpr22 = (Bond_Records[i].Atom2_ff_int_type==Atom2_ff_int_type);
          cmpr12 = (Bond_Records[i].Atom1_ff_int_type==Atom2_ff_int_type);
          cmpr21 = (Bond_Records[i].Atom2_ff_int_type==Atom1_ff_int_type);

          res = 0; 
          if(order==1){  res = (cmpr11 && cmpr22);       }
          else{          res = ( (cmpr11 && cmpr22) || (cmpr12 && cmpr21) );   }


          if(res){ indx = i;  break;   }

       }// if is_Atom1_ff_int_type && is_Atom2_ff_int_type

   }// for i

   return indx;
}


int ForceField::Bond_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Bond_type_index, vector<int>& res){
/****************************************************************
   This function searches for indexes of Bond_Records vector in which
   data about bond formed by 2 atoms of of force field types
   "Atom1_ff_int_type" and "Atom2_ff_int_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of bond_records which match the search parameters
   The bond record will be chosen on the basis of additional comparison
   the Bond_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Bond_Records.size();
   int cmpr11,cmpr12,cmpr21,cmpr22,cmpr;

   for(int i=0;i<sz;i++){

       cmpr11 = cmpr12 = cmpr21 = cmpr22 = 0;
       if(Bond_Records[i].is_Atom1_ff_int_type && Bond_Records[i].is_Atom2_ff_int_type){

          cmpr11 = (Bond_Records[i].Atom1_ff_int_type == Atom1_ff_int_type);
          cmpr22 = (Bond_Records[i].Atom2_ff_int_type == Atom2_ff_int_type);
          cmpr12 = (Bond_Records[i].Atom1_ff_int_type == Atom2_ff_int_type);
          cmpr21 = (Bond_Records[i].Atom2_ff_int_type == Atom1_ff_int_type);

       }// if both Atom1_ff_int_type and Atom2_ff_int_type are defined    
       

       if(Bond_Records[i].is_Bond_type_index){
          
          cmpr = (Bond_Records[i].Bond_type_index == Bond_type_index);
          
       }// if Bond_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if((cmpr11&&cmpr22&&cmpr)||(cmpr12&&cmpr21&&cmpr)){
          indx = i; 
          res.push_back(i);
       }
   }// for i
   return indx;
}


int ForceField::Bond_Record_Index_by_Element(int Atom1_atomic_number, int Atom2_atomic_number, int Bond_type_index, vector<int>& res){
/****************************************************************
   This function searches for indexes of Bond_Records vector in which
   data about bond formed by 2 atoms of of atomic numbers
   "Atom1_atomic_number" and "Atom2_atomic_number" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of bond_records which match the search parameters
   The bond record will be chosen on the basis of additional comparison
   the Bond_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Bond_Records.size();
   int cmpr11,cmpr12,cmpr21,cmpr22,cmpr;

   for(int i=0;i<sz;i++){

       cmpr11 = cmpr12 = cmpr21 = cmpr22 = 0;
       if(Bond_Records[i].is_Atom1_atomic_number && Bond_Records[i].is_Atom2_atomic_number){

          cmpr11 = (Bond_Records[i].Atom1_atomic_number == Atom1_atomic_number);
          cmpr22 = (Bond_Records[i].Atom2_atomic_number == Atom2_atomic_number);
          cmpr12 = (Bond_Records[i].Atom1_atomic_number == Atom2_atomic_number);
          cmpr21 = (Bond_Records[i].Atom2_atomic_number == Atom1_atomic_number);

       }// if both Atom1_ff_int_type and Atom2_ff_int_type are defined


       if(Bond_Records[i].is_Bond_type_index){

          cmpr = (Bond_Records[i].Bond_type_index == Bond_type_index);

       }// if Bond_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if((cmpr11&&cmpr22&&cmpr)||(cmpr12&&cmpr21&&cmpr)){
          indx = i;
          res.push_back(i);
       }
   }// for i
   return indx;
}

int ForceField::Bond_Record_Index_by_Element(std::string Atom1_element, std::string Atom2_element, int Bond_type_index, vector<int>& res){
/****************************************************************
   This function searches for indexes of Bond_Records vector in which
   data about bond formed by 2 atoms of elements
   "Atom1_element" and "Atom2_element" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of bond_records which match the search parameters
   The bond record will be chosen on the basis of additional comparison
   the Bond_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Bond_Records.size();
   int cmpr11,cmpr12,cmpr21,cmpr22,cmpr;

   for(int i=0;i<sz;i++){

       cmpr11 = cmpr12 = cmpr21 = cmpr22 = 0;
       if(Bond_Records[i].is_Atom1_element && Bond_Records[i].is_Atom2_element){

          cmpr11 = (Bond_Records[i].Atom1_element == Atom1_element);
          cmpr22 = (Bond_Records[i].Atom2_element == Atom2_element);
          cmpr12 = (Bond_Records[i].Atom1_element == Atom2_element);
          cmpr21 = (Bond_Records[i].Atom2_element == Atom1_element);

       }// if both Atom1_ff_int_type and Atom2_ff_int_type are defined


       if(Bond_Records[i].is_Bond_type_index){

          cmpr = (Bond_Records[i].Bond_type_index == Bond_type_index);

       }// if Bond_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if((cmpr11&&cmpr22&&cmpr)||(cmpr12&&cmpr21&&cmpr)){
          indx = i;
          res.push_back(i);
       }
   }// for i
   return indx;
}

int ForceField::Bond_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type){
/****************************************************************
   This function searches for indexes of Bond_Records vector in which
   data about bond formed by 2 atoms of of force field types
   "Atom1_ff_type" and "Atom2_ff_type" are stored
   Returns -1 if such index has not been found, otherwise returns the index
   of last of bond_records which match the search parameters
   The bond record will be chosen on the basis of additional comparison
   the Bond_type_index properties if they are available
*****************************************************************/
  int indx = -1;
  int sz   = Bond_Records.size();
  int cmpr11,cmpr12,cmpr21,cmpr22,cmpr;

  for(int i=0;i<sz;i++){
    cmpr11 = cmpr12 = cmpr21 = cmpr22 = 0;
    if(Bond_Records[i].is_Atom1_ff_type && Bond_Records[i].is_Atom2_ff_type){
      cmpr11 = (Bond_Records[i].Atom1_ff_type == Atom1_ff_type);
      cmpr22 = (Bond_Records[i].Atom2_ff_type == Atom2_ff_type);
      cmpr12 = (Bond_Records[i].Atom1_ff_type == Atom2_ff_type);
      cmpr21 = (Bond_Records[i].Atom2_ff_type == Atom1_ff_type);
    }// if both Atom1_ff_int_type and Atom2_ff_int_type are defined
    // --------- Conclusions -------
    if((cmpr11&&cmpr22&&cmpr)||(cmpr12&&cmpr21&&cmpr)){ indx = i;  }
  }// for i
  return indx;
}


int ForceField::Bond_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type, int Bond_type_index, vector<int>& res){
/****************************************************************
   This function searches for indexes of Bond_Records vector in which
   data about bond formed by 2 atoms of of force field types
   "Atom1_ff_type" and "Atom2_ff_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of bond_records which match the search parameters
   The bond record will be chosen on the basis of additional comparison
   the Bond_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Bond_Records.size();
   int cmpr11,cmpr12,cmpr21,cmpr22,cmpr;

   for(int i=0;i<sz;i++){

       cmpr11 = cmpr12 = cmpr21 = cmpr22 = 0;
       if(Bond_Records[i].is_Atom1_ff_type && Bond_Records[i].is_Atom2_ff_type){

          cmpr11 = (Bond_Records[i].Atom1_ff_type == Atom1_ff_type);
          cmpr22 = (Bond_Records[i].Atom2_ff_type == Atom2_ff_type);
          cmpr12 = (Bond_Records[i].Atom1_ff_type == Atom2_ff_type);
          cmpr21 = (Bond_Records[i].Atom2_ff_type == Atom1_ff_type);

       }// if both Atom1_ff_int_type and Atom2_ff_int_type are defined


       if(Bond_Records[i].is_Bond_type_index){

          cmpr = (Bond_Records[i].Bond_type_index == Bond_type_index);

       }// if Bond_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if((cmpr11&&cmpr22&&cmpr)||(cmpr12&&cmpr21&&cmpr)){
          indx = i;
          res.push_back(i);
       }
   }// for i
   return indx;
}





int ForceField::Add_Bond_Record(Bond_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Bond_Record rec in array Bond_Records. If this bond type (bond record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the bond record rec is not valid - will not be added.
   If the Bond_type_index property is defined then it is taken into account
   during check on bond existence
*********************************************************************/
   int res = 1;
   int sz = Bond_Records.size();

   int int_types = (rec.is_Atom1_ff_int_type && rec.is_Atom2_ff_int_type );
   int sym_types = (rec.is_Atom1_ff_type && rec.is_Atom2_ff_type);
   int elt_types = (rec.is_Atom1_atomic_number && rec.is_Atom2_atomic_number);
   int is_bt_indx = rec.is_Bond_type_index;
   int bt_indx = 0;
   if(is_bt_indx){ bt_indx = rec.Bond_type_index; }

   vector<int> bnd_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_types
   last_indx = Bond_Record_Index(rec.Atom1_ff_int_type,rec.Atom2_ff_int_type,bt_indx,bnd_indxs);
   }
   else if(sym_types){
   // Search on the basis of symbolic types
   last_indx = Bond_Record_Index(rec.Atom1_ff_type,rec.Atom2_ff_type,bt_indx,bnd_indxs);
   }
   else if(elt_types){
   // Search on the basis of atomic numbers
   last_indx = Bond_Record_Index_by_Element(rec.Atom1_atomic_number,rec.Atom2_atomic_number,bt_indx,bnd_indxs);
   }
   else{ 
   // Record will not be added to the table of Bond_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2; 
   }

  
   if(last_indx==-2){ // In this case the record was not added to Bond_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The bond record is valid but does not exist in Bond_Records array
                           // so it will be added completely
       Bond_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar bond record is already defined - try to add new properties

       for(int i=0;i<bnd_indxs.size();i++){
           Bond_Records[bnd_indxs[i]].merge(rec);
       }
       res = 0;
   }


/* Old function body

   if(rec.is_Atom1_ff_type&&rec.is_Atom2_ff_type){

      for(int i=0;i<sz;i++){

         if(Bond_Records[i].is_Atom1_ff_type&&Bond_Records[i].is_Atom2_ff_type){
            // This type already exist
            if(((Bond_Records[i].Atom1_ff_type==rec.Atom1_ff_type)&&(Bond_Records[i].Atom2_ff_type==rec.Atom2_ff_type))||
               ((Bond_Records[i].Atom1_ff_type==rec.Atom2_ff_type)&&(Bond_Records[i].Atom2_ff_type==rec.Atom1_ff_type))
              )
              {
               if(Bond_Records[i].is_Bond_type_index){
               res = 0; break; 
               }// is_Bond_type_index
              }// Atom1_ff_type and Atom2_ff_type
         }// is_Atom1_ff_type and is_Atom2_ff_type
      }// for i

   }
   else{
    std::cout<<"To add bond type record to the force field its Atom1_ff_type and Atom2_ff_type should be defined"<<std::endl;
    res = -1;
   }

   if(res==1){
      // Before adding bond record to array - update Atom1(2)_ff_int_type
      sz = Atom_Records.size();
      for(int i=0;i<sz;i++){
          if(Atom_Records[i].is_Atom_ff_int_type){

             // Atom1
             if(Atom_Records[i].Atom_ff_type == rec.Atom1_ff_type){               
                   rec.Atom1_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom1_ff_int_type = 1;                
             }// if Atom1

             // Atom2  
             if(Atom_Records[i].Atom_ff_type == rec.Atom2_ff_type){               
                   rec.Atom2_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom2_ff_int_type = 1;                
             }// if Atom2

          }
      }// for i     
      Bond_Records.push_back(rec);
   }
*/

   return res;
}


int ForceField::Angle_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Atom3_ff_int_type){
/****************************************************************
   This function searches for index of Angle_Records vector in which
   data about angle formed by 3 atoms of of force field types
   "Atom1_ff_int_type", "Atom2_ff_int_type" and "Atom3_ff_int_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Angle_Records.size();

   for(int i=0;i<sz;i++){
       if(Angle_Records[i].is_Atom2_ff_int_type){
          if(Angle_Records[i].Atom2_ff_int_type==Atom2_ff_int_type){
             if(Angle_Records[i].is_Atom1_ff_int_type && Angle_Records[i].is_Atom3_ff_int_type){
                if(((Angle_Records[i].Atom1_ff_int_type==Atom1_ff_int_type)&&(Angle_Records[i].Atom3_ff_int_type==Atom3_ff_int_type) )||
                   ((Angle_Records[i].Atom1_ff_int_type==Atom3_ff_int_type)&&(Angle_Records[i].Atom3_ff_int_type==Atom1_ff_int_type) )
                  )
                  {
                    indx = i;
                    break;
                  } // if side atom types matches
             }// if side atom int_types are defined
          }// if center atom matches
       } // if Atom2_ff_int_type is defined
   }// for i

   return indx;
}

int ForceField::Angle_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Atom3_ff_int_type, int Angle_type_index, int order, vector<int>& res){
/****************************************************************
   This function searches for indexes of Atom_Records vector in which
   data about angle formed by 3 atoms of of force field types
   "Atom1_ff_int_type" "Atom2_ff_int_type" and "Atom3_ff_int_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of angle_records which match the search parameters
   The angle record will be chosen on the basis of additional comparison
   the Angle_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Angle_Records.size();
   int cmpr11,cmpr13,cmpr31,cmpr33,cmpr,cmpr2;

   for(int i=0;i<sz;i++){

       cmpr2 = 0;
       if(Angle_Records[i].is_Atom2_ff_int_type){
           cmpr2 = (Angle_Records[i].Atom2_ff_int_type == Atom2_ff_int_type);
       }

       if(cmpr2){

       cmpr11 = cmpr13 = cmpr31 = cmpr33 = 0;
       if(Angle_Records[i].is_Atom1_ff_int_type && Angle_Records[i].is_Atom3_ff_int_type){

          cmpr11 = (Angle_Records[i].Atom1_ff_int_type == Atom1_ff_int_type);
          cmpr33 = (Angle_Records[i].Atom3_ff_int_type == Atom3_ff_int_type);
          cmpr13 = (Angle_Records[i].Atom1_ff_int_type == Atom3_ff_int_type);
          cmpr31 = (Angle_Records[i].Atom3_ff_int_type == Atom1_ff_int_type);

       }// if both Atom1_ff_int_type and Atom3_ff_int_type are defined


       if(Angle_Records[i].is_Angle_type_index){

          cmpr = (Angle_Records[i].Angle_type_index == Angle_type_index);

       }// if Angle_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if(order==1){
           if(cmpr11&&cmpr33&&cmpr){
               indx = i;
               res.push_back(i);
           }
       }else{
           if((cmpr11&&cmpr33&&cmpr)||(cmpr13&&cmpr31&&cmpr)){
               indx = i;
               res.push_back(i);
           }
       }// else order==0

       }// if cmpr2==1

   }// for i
   return indx;
}


int ForceField::Angle_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type, std::string Atom3_ff_type, int Angle_type_index, int order, vector<int>& res){
/****************************************************************
   order = 0 - no direction 
   order = 1 - should be directed
   This function searches for indexes of Atom_Records vector in which
   data about angle formed by 3 atoms of of force field types
   "Atom1_ff_type" "Atom2_ff_type" and "Atom3_ff_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of angle_records which match the search parameters
   The angle record will be chosen on the basis of additional comparison
   the Angle_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Angle_Records.size();
   int cmpr11,cmpr13,cmpr31,cmpr33,cmpr,cmpr2;

   for(int i=0;i<sz;i++){

       cmpr2 = 0;
       if(Angle_Records[i].is_Atom2_ff_type){
           cmpr2 = (Angle_Records[i].Atom2_ff_type == Atom2_ff_type);
       }

       if(cmpr2){

       cmpr11 = cmpr13 = cmpr31 = cmpr33 = 0;
       if(Angle_Records[i].is_Atom1_ff_type && Angle_Records[i].is_Atom3_ff_type){

          cmpr11 = (Angle_Records[i].Atom1_ff_type == Atom1_ff_type);
          cmpr33 = (Angle_Records[i].Atom3_ff_type == Atom3_ff_type);
          cmpr13 = (Angle_Records[i].Atom1_ff_type == Atom3_ff_type);
          cmpr31 = (Angle_Records[i].Atom3_ff_type == Atom1_ff_type);

       }// if both Atom1_ff_int_type and Atom3_ff_int_type are defined


       if(Angle_Records[i].is_Angle_type_index){

          cmpr = (Angle_Records[i].Angle_type_index == Angle_type_index);

       }// if Angle_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if(order==1){
           if(cmpr11&&cmpr33&&cmpr){
               indx = i;
               res.push_back(i);
           }
       }else{
           if((cmpr11&&cmpr33&&cmpr)||(cmpr13&&cmpr31&&cmpr)){
               indx = i;
               res.push_back(i);
           }
       }// else

       }// if cmpr2==1

   }// for i
   return indx;
}


int ForceField::Angle_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type, std::string Atom3_ff_type){
/****************************************************************
   This function searches for indexes of Atom_Records vector in which
   data about angle formed by 3 atoms of of force field types
   "Atom1_ff_type" "Atom2_ff_type" and "Atom3_ff_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of angle_records which match the search parameters
   The angle record will be chosen on the basis of additional comparison
   the Angle_type_index properties if they are available
*****************************************************************/
  int indx = -1;
  cout<<"In ForceField::Angle_Record_Index\n";

  int sz   = Angle_Records.size();
  int cmpr11,cmpr13,cmpr31,cmpr33,cmpr,cmpr2;
  for(int i=0;i<sz;i++){
    cmpr2 = 0;
    if(Angle_Records[i].is_Atom2_ff_type){
      cmpr2 = (Angle_Records[i].Atom2_ff_type == Atom2_ff_type);
    }
    cout<<"i = "<<i<<" cmpr2 = "<<cmpr2<<endl;
    if(cmpr2){
      cmpr11 = cmpr13 = cmpr31 = cmpr33 = 0;
      if(Angle_Records[i].is_Atom1_ff_type && Angle_Records[i].is_Atom3_ff_type){
        cmpr11 = (Angle_Records[i].Atom1_ff_type == Atom1_ff_type);
        cmpr33 = (Angle_Records[i].Atom3_ff_type == Atom3_ff_type);
        cmpr13 = (Angle_Records[i].Atom1_ff_type == Atom3_ff_type);
        cmpr31 = (Angle_Records[i].Atom3_ff_type == Atom1_ff_type);
      }// if both Atom1_ff_int_type and Atom3_ff_int_type are defined

      cout<<"cmpr11 = "<<cmpr11<<" cmpr33 = "<<cmpr33<<" cmpr13 = "<<cmpr13<<" cmpr31 = "<<cmpr31<<endl;

      // --------- Conclusions -------
//      if((cmpr11&&cmpr33&&cmpr)||(cmpr13&&cmpr31&&cmpr)){ indx = i;  }
      if((cmpr11&&cmpr33)||(cmpr13&&cmpr31)){ indx = i;  }
    }// if cmpr2==1
  }// for i
  return indx;
}



int ForceField::Add_Angle_Record(Angle_Record rec,int order){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Angle_Record rec in array Angle_Records. If this angle type (angle record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the angle record rec is not valid - will not be added.
   If the Angle_type_index property is defined then it is taken into account
   during check on angle existence
   order parameter controls if we should consider order of indices
   (that is ijk is not the same as kji) or to treat them similarly
*********************************************************************/
   int res = 1;
   int sz = Angle_Records.size();

   int int_types = (rec.is_Atom1_ff_int_type && rec.is_Atom2_ff_int_type && rec.is_Atom3_ff_int_type);
   int sym_types = (rec.is_Atom1_ff_type && rec.is_Atom2_ff_type && rec.is_Atom3_ff_type);
//   int elt_types = (rec.is_Atom1_atomic_number && rec.is_Atom2_atomic_number && rec.is_Atom3_atomic_number);
   int is_at_indx = rec.is_Angle_type_index;
   int at_indx = 0;
   if(is_at_indx){ at_indx = rec.Angle_type_index; }

   vector<int> ang_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_types
   last_indx = Angle_Record_Index(rec.Atom1_ff_int_type,rec.Atom2_ff_int_type,rec.Atom3_ff_int_type,at_indx, order, ang_indxs);
   }
   else if(sym_types){
   // Search on the basis of symbolic types
   last_indx = Angle_Record_Index(rec.Atom1_ff_type,rec.Atom2_ff_type, rec.Atom3_ff_type, at_indx, order,ang_indxs);
   }
   else{
   // Record will not be added to the table of Angle_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2;
   }

   if(last_indx==-2){ // In this case the record was not added to Angle_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The angle record is valid but does not exist in Angle_Records array
                           // so it will be added completely
       Angle_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar bond record is already defined - try to add new properties

       for(int i=0;i<ang_indxs.size();i++){
           Angle_Records[ang_indxs[i]].merge(rec);
       }
       res = 0;
   }


   return res;
}



int ForceField::Add_Angle_Record(Angle_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Angle_Record rec in array Angle_Records. If this angle type (angle record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array (we did not add it this time). Returns -1
   if there is an error.
*********************************************************************/
   int res = 1;
   int sz = Angle_Records.size();

   if(rec.is_Atom1_ff_type&&rec.is_Atom2_ff_type&&rec.is_Atom3_ff_type){

      for(int i=0;i<sz;i++){

         if(Angle_Records[i].is_Atom1_ff_type&&Angle_Records[i].is_Atom2_ff_type&&Angle_Records[i].is_Atom3_ff_type){
            // This type already exist
	   if(Angle_Records[i].Atom2_ff_type==rec.Atom2_ff_type){
             if(((Angle_Records[i].Atom1_ff_type==rec.Atom1_ff_type)&&(Angle_Records[i].Atom2_ff_type==rec.Atom3_ff_type))||
                ((Angle_Records[i].Atom1_ff_type==rec.Atom3_ff_type)&&(Angle_Records[i].Atom2_ff_type==rec.Atom1_ff_type))
               ){ res = 0; break; }
	   }
         }
      }// for i

   }
   else{
    std::cout<<"To add angle type record to the force field its Atom1_ff_type, Atom2_ff_type  and Atom3_ff_type should be defined"<<std::endl;
    res = -1;
   }

   if(res==1){
      // Before adding angle record to array - update Atom1(2,3)_ff_int_type
      sz = Atom_Records.size();
      for(int i=0;i<sz;i++){
          if(Atom_Records[i].is_Atom_ff_int_type){

             // Atom1
             if(Atom_Records[i].Atom_ff_type == rec.Atom1_ff_type){               
                   rec.Atom1_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom1_ff_int_type = 1;                
             }// if Atom1

             // Atom2  
             if(Atom_Records[i].Atom_ff_type == rec.Atom2_ff_type){               
                   rec.Atom2_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom2_ff_int_type = 1;                
             }// if Atom2

             // Atom3  
             if(Atom_Records[i].Atom_ff_type == rec.Atom3_ff_type){               
                   rec.Atom3_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom3_ff_int_type = 1;                
             }// if Atom3

          }
      }// for i     
      Angle_Records.push_back(rec);
   }

   return res;
}


int ForceField::Dihedral_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Atom3_ff_int_type, int Atom4_ff_int_type){
/****************************************************************
   This function searches for index of Dihedral_Records vector in which
   data about dihedral formed by 4 atoms of of force field types
   "Atom1_ff_int_type", "Atom2_ff_int_type","Atom3_ff_int_type" and 
   "Atom4_ff_int_type" are stored
    Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Dihedral_Records.size();

   for(int i=0;i<sz;i++){
       if(Dihedral_Records[i].is_Atom2_ff_int_type&&Dihedral_Records[i].is_Atom3_ff_int_type){

          if((Dihedral_Records[i].Atom2_ff_int_type==Atom2_ff_int_type)&&(Dihedral_Records[i].Atom3_ff_int_type==Atom3_ff_int_type)){
             if(Dihedral_Records[i].is_Atom1_ff_int_type && Dihedral_Records[i].is_Atom4_ff_int_type){
                if((Dihedral_Records[i].Atom1_ff_int_type==Atom1_ff_int_type)&&(Dihedral_Records[i].Atom4_ff_int_type==Atom4_ff_int_type) )                  
                  {
                    indx = i;
                    break;
                  } // if side atom types matches: 1=1 and 4=4
             }// if side atoms int_types are defined
          }// if center 2 atoms matches : 2=2 and 3=3

          if((Dihedral_Records[i].Atom2_ff_int_type==Atom3_ff_int_type)&&(Dihedral_Records[i].Atom3_ff_int_type==Atom2_ff_int_type)){
             if(Dihedral_Records[i].is_Atom1_ff_int_type && Dihedral_Records[i].is_Atom4_ff_int_type){
                if((Dihedral_Records[i].Atom1_ff_int_type==Atom4_ff_int_type)&&(Dihedral_Records[i].Atom4_ff_int_type==Atom1_ff_int_type) )
                  {
                    indx = i;
                    break;
                  } // if side atom types matches: 1=4 and 4=1
             }// if side atoms int_types are defined
          }// if center 2 atoms matches : 2=3 and 3=2


       } // if Atom2_ff_int_type and Atom3_ff_int_type are defined
   }// for i

   return indx;
}

int ForceField::Improper_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Atom3_ff_int_type, int Atom4_ff_int_type){
/****************************************************************
   This function searches for index of Improper_Records vector in which
   data about improper formed by 4 atoms of of force field types
   "Atom1_ff_int_type", "Atom2_ff_int_type","Atom3_ff_int_type" and
   "Atom4_ff_int_type" are stored
    Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int res = 0;
   int sz   = Improper_Records.size();

   for(int i=0;i<sz;i++){
       if(Improper_Records[i].is_Atom1_ff_int_type &&
          Improper_Records[i].is_Atom2_ff_int_type &&
          Improper_Records[i].is_Atom3_ff_int_type &&
          Improper_Records[i].is_Atom4_ff_int_type)
       {

          res = 0;
          if(Improper_Records[i].Atom2_ff_int_type==Atom2_ff_int_type){

          // Now search match of any of 6 permutations of 3 periferal atoms
              if(Improper_Records[i].Atom1_ff_int_type==Atom1_ff_int_type){
                  if(Improper_Records[i].Atom3_ff_int_type==Atom3_ff_int_type && Improper_Records[i].Atom4_ff_int_type==Atom4_ff_int_type){
                      res = 1;
                  }
                  else if(Improper_Records[i].Atom3_ff_int_type==Atom4_ff_int_type && Improper_Records[i].Atom4_ff_int_type==Atom3_ff_int_type){
                      res = 1;
                  }
              }

              else if(Improper_Records[i].Atom3_ff_int_type==Atom3_ff_int_type){
                  if(Improper_Records[i].Atom1_ff_int_type==Atom1_ff_int_type && Improper_Records[i].Atom4_ff_int_type==Atom4_ff_int_type){
                      res = 1;
                  }
                  else if(Improper_Records[i].Atom1_ff_int_type==Atom4_ff_int_type && Improper_Records[i].Atom4_ff_int_type==Atom1_ff_int_type){
                      res = 1;
                  }
              }

              else if(Improper_Records[i].Atom4_ff_int_type==Atom4_ff_int_type){
                  if(Improper_Records[i].Atom3_ff_int_type==Atom3_ff_int_type && Improper_Records[i].Atom1_ff_int_type==Atom1_ff_int_type){
                      res = 1;
                  }
                  else if(Improper_Records[i].Atom3_ff_int_type==Atom1_ff_int_type && Improper_Records[i].Atom1_ff_int_type==Atom3_ff_int_type){
                      res = 1;
                  }
              }

          }// if 2 == 2 (central atoms are the same)

          if(res){ indx = i; break; }

       }// if all in types are defined

   }// for i

   return indx;
}



int ForceField::Dihedral_Record_Index(int Atom2_ff_int_type, int Atom3_ff_int_type){
/****************************************************************
   This function searches for index of Dihedral_Records vector based on
   only middle 2 atoms of force field types "Atom2_ff_int_type" and 
   "Atom3_ff_int_type". This is necessary because of in many cases
   the dihedral parmeter are not sensitive to 2 terminal atoms.
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Dihedral_Records.size();

   for(int i=0;i<sz;i++){
       if(Dihedral_Records[i].is_Atom2_ff_int_type&&Dihedral_Records[i].is_Atom3_ff_int_type){

          if(((Dihedral_Records[i].Atom2_ff_int_type==Atom2_ff_int_type)&&(Dihedral_Records[i].Atom3_ff_int_type==Atom3_ff_int_type))|| 
             ((Dihedral_Records[i].Atom2_ff_int_type==Atom3_ff_int_type)&&(Dihedral_Records[i].Atom3_ff_int_type==Atom2_ff_int_type))
            ){
                    indx = i;
                    break;              

             }// if center 2 atoms matches : 2=2, 3=3 or 2=3, 3=2


       } // if Atom2_ff_int_type and Atom3_ff_int_type are defined
   }// for i

   return indx;
}

int ForceField::Dihedral_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Atom3_ff_int_type, int Atom4_ff_int_type, int Dihedral_type_index, int order, vector<int>& res){
/****************************************************************
   This function searches for indexes of Dihedral_Records vector in which
   data about dihedral formed by 4 atoms of of force field types
   "Atom1_ff_int_type" "Atom2_ff_int_type" "Atom3_ff_int_type" and "
   Atom4_ff_int_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of dihedral_records which match the search parameters
   The dihedral record will be chosen on the basis of additional comparison
   the Dihedral_type_index properties if they are available
   order - is a parameter defining importance of indices order
   if order = 1 => ijkl is not the same as lkji
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Dihedral_Records.size();
   int cmpr11,cmpr14,cmpr41,cmpr23,cmpr32,cmpr22,cmpr33,cmpr44,cmpr;

   for(int i=0;i<sz;i++){

       if(Dihedral_Records[i].is_Dihedral_type_index){

          cmpr = (Dihedral_Records[i].Dihedral_type_index == Dihedral_type_index);

       }// if Dihedral_type_index is defined
       else{ cmpr = 1; }


       cmpr22 = cmpr33 = 0;
       if(Dihedral_Records[i].is_Atom1_ff_int_type &&
          Dihedral_Records[i].is_Atom2_ff_int_type &&
          Dihedral_Records[i].is_Atom3_ff_int_type &&
          Dihedral_Records[i].is_Atom4_ff_int_type
       ){
           cmpr11 = (Dihedral_Records[i].Atom1_ff_int_type == Atom1_ff_int_type);
           cmpr22 = (Dihedral_Records[i].Atom2_ff_int_type == Atom2_ff_int_type);
           cmpr33 = (Dihedral_Records[i].Atom3_ff_int_type == Atom3_ff_int_type);
           cmpr23 = (Dihedral_Records[i].Atom2_ff_int_type == Atom3_ff_int_type);
           cmpr32 = (Dihedral_Records[i].Atom3_ff_int_type == Atom2_ff_int_type);
           cmpr44 = (Dihedral_Records[i].Atom4_ff_int_type == Atom4_ff_int_type);
           cmpr14 = (Dihedral_Records[i].Atom1_ff_int_type == Atom4_ff_int_type);
           cmpr41 = (Dihedral_Records[i].Atom4_ff_int_type == Atom1_ff_int_type);

       }


       if(order==1){     

           if(cmpr&&cmpr11&&cmpr22&&cmpr33&&cmpr44){
               indx = i;
               res.push_back(i);
           }

       }else{

           if( (cmpr&&cmpr11&&cmpr22&&cmpr33&&cmpr44) ||
               (cmpr&&cmpr14&&cmpr23&&cmpr32&&cmpr41)
             ){
               indx = i;
               res.push_back(i);
           }

       }// order != 1


   }// for i
   return indx;
}

int ForceField::Dihedral_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type, std::string Atom3_ff_type, std::string Atom4_ff_type){
/****************************************************************
   This function searches for indexes of Dihedral_Records vector in which
   data about dihedral formed by 4 atoms of of force field types
   "Atom1_ff_int_type" "Atom2_ff_int_type" "Atom3_ff_int_type" and "
   Atom4_ff_int_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of dihedral_records which match the search parameters
   The dihedral record will be chosen on the basis of additional comparison
   the Dihedral_type_index properties if they are available
   order - is a parameter defining importance of indices order
   if order = 1 => ijkl is not the same as lkji
*****************************************************************/
  vector<int> tmp;
  return Dihedral_Record_Index(Atom1_ff_type,Atom2_ff_type,Atom3_ff_type,Atom4_ff_type,-1,0,tmp);
}

int ForceField::Dihedral_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type, std::string Atom3_ff_type, std::string Atom4_ff_type, int Dihedral_type_index, int order, vector<int>& res){
/****************************************************************
   This function searches for indexes of Dihedral_Records vector in which
   data about dihedral formed by 4 atoms of of force field types
   "Atom1_ff_int_type" "Atom2_ff_int_type" "Atom3_ff_int_type" and "
   Atom4_ff_int_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of dihedral_records which match the search parameters
   The dihedral record will be chosen on the basis of additional comparison
   the Dihedral_type_index properties if they are available
   order - is a parameter defining importance of indices order
   if order = 1 => ijkl is not the same as lkji
   If Dihedral_type_index == -1 it will not be taken into account
*****************************************************************/
  if(res.size()>0) { res.clear(); }
  int indx = -1;
  int sz   = Dihedral_Records.size();
  int cmpr11,cmpr14,cmpr41,cmpr23,cmpr32,cmpr22,cmpr33,cmpr44,cmpr;
  for(int i=0;i<sz;i++){
    //--------------------------------------------------
    if(Dihedral_type_index==-1){ cmpr = 1; }
    else {
      if(Dihedral_Records[i].is_Dihedral_type_index){
        cmpr = (Dihedral_Records[i].Dihedral_type_index == Dihedral_type_index);
      }// if Dihedral_type_index is defined
      else{ cmpr = 1; }
    }

    //----------------------------------------------------
    cmpr22 = cmpr33 = 0;
       if(Dihedral_Records[i].is_Atom1_ff_type &&
          Dihedral_Records[i].is_Atom2_ff_type &&
          Dihedral_Records[i].is_Atom3_ff_type &&
          Dihedral_Records[i].is_Atom4_ff_type
       ){
           cmpr11 = (Dihedral_Records[i].Atom1_ff_type == Atom1_ff_type);
           cmpr22 = (Dihedral_Records[i].Atom2_ff_type == Atom2_ff_type);
           cmpr33 = (Dihedral_Records[i].Atom3_ff_type == Atom3_ff_type);
           cmpr23 = (Dihedral_Records[i].Atom2_ff_type == Atom3_ff_type);
           cmpr32 = (Dihedral_Records[i].Atom3_ff_type == Atom2_ff_type);
           cmpr44 = (Dihedral_Records[i].Atom4_ff_type == Atom4_ff_type);
           cmpr14 = (Dihedral_Records[i].Atom1_ff_type == Atom4_ff_type);
           cmpr41 = (Dihedral_Records[i].Atom4_ff_type == Atom1_ff_type);
       }


       if(order==1){

           if(cmpr&&cmpr11&&cmpr22&&cmpr33&&cmpr44){
               indx = i;
               res.push_back(i);
           }

       }else{

           if( (cmpr&&cmpr11&&cmpr22&&cmpr33&&cmpr44) ||
               (cmpr&&cmpr14&&cmpr23&&cmpr32&&cmpr41)
             ){
               indx = i;
               res.push_back(i);
           }

       }// order != 1


   }// for i
   return indx;
}

int ForceField::Add_Dihedral_Record(Dihedral_Record rec,int order){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Dihedral_Record rec in array Dihedral_Records. If this dihedral type 
   (dihedral record) exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the dihedral record rec is not valid - will not be added.
   If the Dihedral_type_index property is defined then it is taken into account
   during check on dihedral existence
   order parameter controls if we should consider order of indices
   (that is ijkl is not the same as lkji) or to treat them similarly
*********************************************************************/
   int res = 1;
   int sz = Dihedral_Records.size();

   int int_types = (rec.is_Atom1_ff_int_type && rec.is_Atom2_ff_int_type && rec.is_Atom3_ff_int_type && rec.is_Atom4_ff_int_type);
   int sym_types = (rec.is_Atom1_ff_type && rec.is_Atom2_ff_type && rec.is_Atom3_ff_type && rec.is_Atom4_ff_type);
//   int elt_types = (rec.is_Atom1_atomic_number && rec.is_Atom2_atomic_number && rec.is_Atom3_atomic_number);
   int is_dt_indx = rec.is_Dihedral_type_index;
   int dt_indx = 0;
   if(is_dt_indx){ dt_indx = rec.Dihedral_type_index; }

   vector<int> dih_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_types
   last_indx = Dihedral_Record_Index(rec.Atom1_ff_int_type,rec.Atom2_ff_int_type,rec.Atom3_ff_int_type, rec.Atom4_ff_int_type,dt_indx, order, dih_indxs);
   }
   else if(sym_types){
   // Search on the basis of symbolic types
   last_indx = Dihedral_Record_Index(rec.Atom1_ff_type,rec.Atom2_ff_type, rec.Atom3_ff_type, rec.Atom4_ff_type,dt_indx, order,dih_indxs);
   }
   else{
   // Record will not be added to the table of Dihedral_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2;
   }

   if(last_indx==-2){ // In this case the record was not added to Dihedral_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The dihedral record is valid but does not exist in Dihedral_Records array
                           // so it will be added completely
       Dihedral_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar bond record is already defined - try to add new properties

       for(int i=0;i<dih_indxs.size();i++){
           Dihedral_Records[dih_indxs[i]].merge(rec);
       }
       res = 0;
   }


   return res;
}





int ForceField::Add_Dihedral_Record(Dihedral_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Dihedral_Record rec in array Dihedral_Records. If this dihedral type (dihedral record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array (we did not add it this time). Returns -1
   if there is an error.
*********************************************************************/
   int res = 1;
   int sz = Dihedral_Records.size();

   if(rec.is_Atom1_ff_type&&rec.is_Atom2_ff_type&&rec.is_Atom3_ff_type&&rec.is_Atom4_ff_type){

      for(int i=0;i<sz;i++){

         if(Dihedral_Records[i].is_Atom1_ff_type&&Dihedral_Records[i].is_Atom2_ff_type
          &&Dihedral_Records[i].is_Atom3_ff_type&&Dihedral_Records[i].is_Atom4_ff_type){
           // This type already exist
	   if((Dihedral_Records[i].Atom2_ff_type==rec.Atom2_ff_type)&&(Dihedral_Records[i].Atom3_ff_type==rec.Atom3_ff_type)){
             if((Dihedral_Records[i].Atom1_ff_type==rec.Atom1_ff_type)&&(Dihedral_Records[i].Atom4_ff_type==rec.Atom4_ff_type))               
               { res = 0; break; }
	   }
           if((Dihedral_Records[i].Atom2_ff_type==rec.Atom3_ff_type)&&(Dihedral_Records[i].Atom3_ff_type==rec.Atom2_ff_type)){
             if((Dihedral_Records[i].Atom1_ff_type==rec.Atom4_ff_type)&&(Dihedral_Records[i].Atom4_ff_type==rec.Atom1_ff_type))               
               { res = 0; break; }
	   }

         }
      }// for i

   }
   else{
    std::cout<<"To add dihedral type record to the force field its Atom1_ff_type, Atom2_ff_type, Atom3_ff_type and Atom4_ff_type should be defined"<<std::endl;
    res = -1;
   }

   if(res==1){
      // Before adding angle record to array - update Atom1(2,3,4)_ff_int_type
      sz = Atom_Records.size();
      for(int i=0;i<sz;i++){
          if(Atom_Records[i].is_Atom_ff_int_type){

             // Atom1
             if(Atom_Records[i].Atom_ff_type == rec.Atom1_ff_type){               
                   rec.Atom1_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom1_ff_int_type = 1;                
             }// if Atom1

             // Atom2  
             if(Atom_Records[i].Atom_ff_type == rec.Atom2_ff_type){               
                   rec.Atom2_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom2_ff_int_type = 1;                
             }// if Atom2

             // Atom3  
             if(Atom_Records[i].Atom_ff_type == rec.Atom3_ff_type){               
                   rec.Atom3_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom3_ff_int_type = 1;                
             }// if Atom3

             // Atom4  
             if(Atom_Records[i].Atom_ff_type == rec.Atom4_ff_type){               
                   rec.Atom4_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom4_ff_int_type = 1;                
             }// if Atom4

          }
      }// for i     
      Dihedral_Records.push_back(rec);
   }

   return res;
}

int ForceField::Add_Improper_Record(Dihedral_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Dihedral_Record rec in array Improper_Records. If this improper type
   (improper record) exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the improper record rec is not valid - will not be added.

*********************************************************************/
   int res = 1;
   int sz = Improper_Records.size();

   int int_types = (rec.is_Atom1_ff_int_type && rec.is_Atom2_ff_int_type && rec.is_Atom3_ff_int_type && rec.is_Atom4_ff_int_type);
   int sym_types = (rec.is_Atom1_ff_type && rec.is_Atom2_ff_type && rec.is_Atom3_ff_type && rec.is_Atom4_ff_type);

   int is_dt_indx = rec.is_Dihedral_type_index;
   int dt_indx = 0;
   if(is_dt_indx){ dt_indx = rec.Dihedral_type_index; }

   vector<int> dih_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_types
   last_indx = Improper_Record_Index(rec.Atom1_ff_int_type,rec.Atom2_ff_int_type,rec.Atom3_ff_int_type, rec.Atom4_ff_int_type);
   dih_indxs.push_back(last_indx);
   }
/*
   else if(sym_types){
   // Search on the basis of symbolic types
   last_indx = Dihedral_Record_Index(rec.Atom1_ff_type,rec.Atom2_ff_type, rec.Atom3_ff_type, rec.Atom4_ff_type,dt_indx, order,dih_indxs);
   }
*/
   else{
   // Record will not be added to the table of Improper_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2;
   }


   if(last_indx==-2){ // In this case the record was not added to Improper_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The dihedral record is valid but does not exist in Improper_Records array
                           // so it will be added completely
       Improper_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar improper record is already defined - try to add new properties

       for(int i=0;i<dih_indxs.size();i++){
           Improper_Records[dih_indxs[i]].merge(rec);
       }
       res = 0;
   }


   return res;
}


int ForceField::Fragment_Record_Index(int Fragment_ff_int_type){
/****************************************************************
   This function searches for index of Fragment_Records vector in which
   data about fragment of force field type "Fragment_ff_int_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Fragment_Records.size();

   for(int i=0;i<sz;i++){
       if(Fragment_Records[i].is_Fragment_ff_int_type){
          if(Fragment_Records[i].Fragment_ff_int_type==Fragment_ff_int_type){
             indx = i;
             break;
          }
       }
   }

   return indx;
}

int ForceField::Fragment_Record_Index(std::string Fragment_ff_type){
/****************************************************************
   This function searches for index of Fragment_Records vector in which
   data about fragment of force field type "Fragment_ff_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Fragment_Records.size();

   for(int i=0;i<sz;i++){
       if(Fragment_Records[i].is_Fragment_ff_type){
          if(Fragment_Records[i].Fragment_ff_type==Fragment_ff_type){
             indx = i;
             break;
          }
       }
   }

   return indx;
}

int ForceField::Add_Fragment_Record(Fragment_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Fragment_Record rec in array Fragment_Records. If this fragment type (fragment record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the atom record rec is not valid - will not be added.
*********************************************************************/
   int res = 1;
   int sz = Fragment_Records.size();

   int int_types = (rec.is_Fragment_ff_int_type);
   int sym_types = (rec.is_Fragment_ff_type);

   vector<int> at_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_type
   last_indx = Fragment_Record_Index(rec.Fragment_ff_int_type);
   at_indxs.push_back(last_indx);
   }
   else if(sym_types){
   // Search on the basis of symbolic type
   last_indx = Fragment_Record_Index(rec.Fragment_ff_type);
   at_indxs.push_back(last_indx);
   }
   else{
   // Record will not be added to the table of Fragment_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2;
   }

   if(last_indx==-2){ // In this case the record was not added to Fragment_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The atom record is valid but does not exist in Fragment_Records array
                           // so it will be added completely
       Fragment_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar atom record is already defined - try to add new properties

       for(int i=0;i<at_indxs.size();i++){
           Fragment_Records[at_indxs[i]].merge(rec);
       }
       res = 0;
   }

   return res;

}



int ForceField::show_atom_records(){

    for(int i=0;i<Atom_Records.size();i++){

        Atom_Records[i].show_info();
    }

    return 1;
}

int ForceField::show_bond_records(){

    for(int i=0;i<Bond_Records.size();i++){

        Bond_Records[i].show_info();
    }

    return 1;
}

int ForceField::show_angle_records(){

    for(int i=0;i<Angle_Records.size();i++){

        Angle_Records[i].show_info();
    }

    return 1;
}

int ForceField::show_dihedral_records(){

    for(int i=0;i<Dihedral_Records.size();i++){

        Dihedral_Records[i].show_info();
    }

    return 1;
}

int ForceField::show_improper_records(){

    for(int i=0;i<Improper_Records.size();i++){

        Improper_Records[i].show_info();
    }

    return 1;
}

int ForceField::show_fragment_records(){

    for(int i=0;i<Fragment_Records.size();i++){

        Fragment_Records[i].show_info();
    }

    return 1;
}


}// namespace libforcefield
}// namespace liblibra



