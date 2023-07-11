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
/**
  \file Element.cpp
  \brief The file describes the Element class and related functions    
*/

#include "Element.h"

/// liblibra namespace
namespace liblibra{

using namespace libio;

/// libchemobjects namespace
namespace libchemobjects{

/// libuniverse namespace
namespace libuniverse{


void Element::init_variables(){
  is_Elt_name = 0;
  is_Elt_mass = 0;
  is_Elt_number = 0;
  is_Elt_nucleus_charge = 0;
  is_Elt_period = 0;
  is_Elt_group = 0; 
  is_Elt_block = 0; 
  is_Elt_red = 0;
  is_Elt_green = 0;
  is_Elt_blue = 0;
  is_Elt_bond_length = 0;
  is_Elt_radius = 0;

}

void Element::copy_content(const Element& at){

  if(at.is_Elt_name){ Elt_name = at.Elt_name;  is_Elt_name = 1;}
  if(at.is_Elt_mass){ Elt_mass = at.Elt_mass;  is_Elt_mass = 1;}
  if(at.is_Elt_number){ Elt_number = at.Elt_number;  is_Elt_number = 1;}
  if(at.is_Elt_nucleus_charge){ Elt_nucleus_charge = at.Elt_nucleus_charge;  is_Elt_nucleus_charge = 1;}
  if(at.is_Elt_period){ Elt_period = at.Elt_period;  is_Elt_period = 1;}
  if(at.is_Elt_group){ Elt_group = at.Elt_group;  is_Elt_group = 1;}
  if(at.is_Elt_block){ Elt_block = at.Elt_block;  is_Elt_block = 1;}
  if(at.is_Elt_red)  { Elt_red   = at.Elt_red;   is_Elt_red = 1;}
  if(at.is_Elt_green){ Elt_green = at.Elt_green; is_Elt_green = 1;}
  if(at.is_Elt_blue) { Elt_blue  = at.Elt_blue;  is_Elt_blue = 1;}
  if(at.is_Elt_bond_length){ Elt_bond_length = at.Elt_bond_length; is_Elt_bond_length = 1; }
  if(at.is_Elt_radius){ Elt_radius = at.Elt_radius; is_Elt_radius = 1; }
 
}

Element::Element(){
/**
  \brief Constructor
  Initialize variables to default values
*/
  // Initialize variables to default values
  init_variables();
}

Element::Element(const Element& at){
/**
  \brief Copy constructor
  \param[in] at The input object

  Initialize variables to default values
  Copy content of the at object which is defined
*/
  // Initialize variables to default values
  init_variables();
  // Copy content of at object which is defined
  copy_content(at);
}

Element& Element::operator=(const Element& at){
/**
  \brief Assignment operator
  Initialize variables to default values
  Copy content of at object which is defined
*/
  // Initialize variables to default values
  init_variables();
  // Copy content of at object which is defined
  copy_content(at);
  return *this;
}

Element::~Element(){ }


void Element::set(object at){
/** 
  \brief Set properties of the Element object from an arbitrary Python object.

  \param[in] at The input object - must contain the members with the names that match the names of the internal variables.
*/

  set_value(is_Elt_name,    Elt_name,    at,"Elt_name");
  set_value(is_Elt_mass,    Elt_mass,    at,"Elt_mass");  
  set_value(is_Elt_number,  Elt_number,  at,"Elt_number");
  set_value(is_Elt_nucleus_charge,Elt_nucleus_charge,at,"Elt_nucleus_charge");
  set_value(is_Elt_period,  Elt_period,  at,"Elt_period");
  set_value(is_Elt_group,   Elt_group,   at,"Elt_group");
  set_value(is_Elt_block,   Elt_block,   at,"Elt_block");
  set_value(is_Elt_red,     Elt_red,     at,"Elt_red");
  set_value(is_Elt_green,   Elt_green,   at,"Elt_green");
  set_value(is_Elt_blue,    Elt_blue,    at,"Elt_blue");
  set_value(is_Elt_bond_length,Elt_bond_length, at, "Elt_bond_length");
  set_value(is_Elt_radius,  Elt_radius,  at,"Elt_radius");
}


void Element::show_info(){
/** 
  \brief Show info about Element properties
*/

  if(is_Elt_name) {std::cout<<"Elt_name = "<<Elt_name<<std::endl;   } 
  if(is_Elt_mass) {std::cout<<"Elt_mass = "<<Elt_mass<<std::endl;   }
  if(is_Elt_number) {std::cout<<"Elt_number = "<<Elt_number<<std::endl;}
  if(is_Elt_nucleus_charge) {std::cout<<"Elt_nucleus_charge = "<<Elt_nucleus_charge<<std::endl;}
  if(is_Elt_period){std::cout<<"Elt_period = "<<Elt_period<<std::endl;}
  if(is_Elt_group){std::cout<<"Elt_group = "<<Elt_group<<std::endl;}
  if(is_Elt_block){std::cout<<"Elt_block = "<<Elt_block<<std::endl;}
  if(is_Elt_red)  {std::cout<<"Elt_red = "<<Elt_red<<std::endl;}
  if(is_Elt_green){std::cout<<"Elt_green = "<<Elt_green<<std::endl;}
  if(is_Elt_blue) {std::cout<<"Elt_blue = "<<Elt_blue<<std::endl;}
  if(is_Elt_bond_length) {std::cout<<"Elt_bond_length = "<<Elt_bond_length<<std::endl;}
  if(is_Elt_radius){std::cout<<"Elt_radius = "<<Elt_radius<<std::endl; }
  std::cout<<std::endl;
}

void Element::save(boost::property_tree::ptree& pt,std::string path){
/**
  \brief Save the state of the Element object as a property tree

  Each defined data member is added as a node to the property tree. The nodes are added to 
  the level of the tree controlled by the path variable.
 
  \param[in,out] pt The property tree to which the properties of the Element are added
  \param[in] path The parameter controlling the level of the tree to which the Element members will be added.
*/

  if(is_Elt_name){  libio::save(pt,path+".Elt_name",Elt_name);    }
  if(is_Elt_id){  libio::save(pt,path+".Elt_id",Elt_id);    }
  if(is_Elt_mass){  libio::save(pt,path+".Elt_mass",Elt_mass);    }
  if(is_Elt_number){  libio::save(pt,path+".Elt_number",Elt_number);    }
  if(is_Elt_nucleus_charge){  libio::save(pt,path+".Elt_nucleus_charge",Elt_nucleus_charge);    }
  if(is_Elt_period){  libio::save(pt,path+".Elt_period",Elt_period);    }
  if(is_Elt_group){  libio::save(pt,path+".Elt_group",Elt_group);    }
  if(is_Elt_block){  libio::save(pt,path+".Elt_block",Elt_block);    }
  if(is_Elt_red){  libio::save(pt,path+".Elt_red",Elt_red);    }
  if(is_Elt_green){  libio::save(pt,path+".Elt_green",Elt_green);    }
  if(is_Elt_blue){  libio::save(pt,path+".Elt_blue",Elt_blue);    }
  if(is_Elt_bond_length){  libio::save(pt,path+".Elt_bond_length",Elt_bond_length);    }
  if(is_Elt_radius){  libio::save(pt,path+".Elt_radius",Elt_radius);    }

}

void save(boost::property_tree::ptree& pt,std::string path,vector<Element>& vt){
/**
  \brief Save the state of the vector of Element objects as a property tree

  Each Element object is added as a separate branch. 
 
  \param[in,out] pt The property tree to which the list of the Element objects will be added
  \param[in] path The parameter controlling the level of the tree to which the list of Element will be added.
  \param[in] vt The list of Element objects to be printed out into property tree
*/

  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".Element"+rt);
  }
}


void Element::load(boost::property_tree::ptree& pt,std::string path,int& status){
/**
  \brief Load the state of the Element object from a property tree

  Each data member found in the property tree is extracted as the member of the Barostat object. The
  status of each found data member is set to 1.
 
  \param[in] pt The property tree from which the properties of the Element will be extracted
  \param[in] path The parameter controlling from which level of the tree we try to extract the Element object
  \param[out] status Is the global status of the success of the operation. It is 1 is at least one Element member is found at
              given level of the property tree.
*/

  int st;
  status = 0;

  libio::load(pt,path+".Elt_name",Elt_name,is_Elt_name); if(is_Elt_name==1) { status=1;}
  libio::load(pt,path+".Elt_id",Elt_id,is_Elt_id); if(is_Elt_id==1) { status=1;}
  libio::load(pt,path+".Elt_mass",Elt_mass,is_Elt_mass); if(is_Elt_mass==1) { status=1;}
  libio::load(pt,path+".Elt_number",Elt_number,is_Elt_number); if(is_Elt_number==1) { status=1;}
  libio::load(pt,path+".Elt_nucleus_charge",Elt_nucleus_charge,is_Elt_nucleus_charge); if(is_Elt_nucleus_charge==1) { status=1;}
  libio::load(pt,path+".Elt_period",Elt_period,is_Elt_period); if(is_Elt_period==1) { status=1;}
  libio::load(pt,path+".Elt_group",Elt_group,is_Elt_group); if(is_Elt_group==1) { status=1;}
  libio::load(pt,path+".Elt_block",Elt_block,is_Elt_block); if(is_Elt_block==1) { status=1;}
  libio::load(pt,path+".Elt_red",Elt_red,is_Elt_red); if(is_Elt_red==1) { status=1;}
  libio::load(pt,path+".Elt_green",Elt_green,is_Elt_green); if(is_Elt_green==1) { status=1;}
  libio::load(pt,path+".Elt_blue",Elt_blue,is_Elt_blue); if(is_Elt_blue==1) { status=1;}
  libio::load(pt,path+".Elt_bond_length",Elt_bond_length,is_Elt_bond_length); if(is_Elt_bond_length==1) { status=1;}
  libio::load(pt,path+".Elt_radius",Elt_radius,is_Elt_radius); if(is_Elt_radius==1) { status=1;}

}

void load(boost::property_tree::ptree& pt,std::string path,vector<Element>& vt,int& status){
/**
  \brief Load the vector of Element objects from a property tree

  Each Element object is extracted from a separate branch. 
 
  \param[in] pt The property tree from which the vector of Element objects will be extracted
  \param[in] path The parameter controlling from which level of the property tree we will try to extract the vector of Element objects
  \param[out] status Is the global status of the success of the operation. It is 1 is at least one Element object is extracted
*/

  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      Element x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}


}// namespace libuniverse
}// namespace libchemobjects

}// liblibra
