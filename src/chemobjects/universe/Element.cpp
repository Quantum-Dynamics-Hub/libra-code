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

#include "Element.h"

namespace libchemobjects{
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
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
}

Element::Element(const Element& at){
  /********************
    Copy constructor
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of at object which is defined
  copy_content(at);
}

Element& Element::operator=(const Element& at){
  /********************
    Assignment operator
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of at object which is defined
  copy_content(at);
  return *this;
}

Element::~Element(){ }


void Element::set(object at){
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

  if(is_Elt_name){  ::save(pt,path+".Elt_name",Elt_name);    }
  if(is_Elt_id){  ::save(pt,path+".Elt_id",Elt_id);    }
  if(is_Elt_mass){  ::save(pt,path+".Elt_mass",Elt_mass);    }
  if(is_Elt_number){  ::save(pt,path+".Elt_number",Elt_number);    }
  if(is_Elt_nucleus_charge){  ::save(pt,path+".Elt_nucleus_charge",Elt_nucleus_charge);    }
  if(is_Elt_period){  ::save(pt,path+".Elt_period",Elt_period);    }
  if(is_Elt_group){  ::save(pt,path+".Elt_group",Elt_group);    }
  if(is_Elt_block){  ::save(pt,path+".Elt_block",Elt_block);    }
  if(is_Elt_red){  ::save(pt,path+".Elt_red",Elt_red);    }
  if(is_Elt_green){  ::save(pt,path+".Elt_green",Elt_green);    }
  if(is_Elt_blue){  ::save(pt,path+".Elt_blue",Elt_blue);    }
  if(is_Elt_bond_length){  ::save(pt,path+".Elt_bond_length",Elt_bond_length);    }
  if(is_Elt_radius){  ::save(pt,path+".Elt_radius",Elt_radius);    }

}

void save(boost::property_tree::ptree& pt,std::string path,vector<Element>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".Element"+rt);
  }
}


void Element::load(boost::property_tree::ptree& pt,std::string path,int& status){
  int st;
  status = 0;

  ::load(pt,path+".Elt_name",Elt_name,is_Elt_name); if(is_Elt_name==1) { status=1;}
  ::load(pt,path+".Elt_id",Elt_id,is_Elt_id); if(is_Elt_id==1) { status=1;}
  ::load(pt,path+".Elt_mass",Elt_mass,is_Elt_mass); if(is_Elt_mass==1) { status=1;}
  ::load(pt,path+".Elt_number",Elt_number,is_Elt_number); if(is_Elt_number==1) { status=1;}
  ::load(pt,path+".Elt_nucleus_charge",Elt_nucleus_charge,is_Elt_nucleus_charge); if(is_Elt_nucleus_charge==1) { status=1;}
  ::load(pt,path+".Elt_period",Elt_period,is_Elt_period); if(is_Elt_period==1) { status=1;}
  ::load(pt,path+".Elt_group",Elt_group,is_Elt_group); if(is_Elt_group==1) { status=1;}
  ::load(pt,path+".Elt_block",Elt_block,is_Elt_block); if(is_Elt_block==1) { status=1;}
  ::load(pt,path+".Elt_red",Elt_red,is_Elt_red); if(is_Elt_red==1) { status=1;}
  ::load(pt,path+".Elt_green",Elt_green,is_Elt_green); if(is_Elt_green==1) { status=1;}
  ::load(pt,path+".Elt_blue",Elt_blue,is_Elt_blue); if(is_Elt_blue==1) { status=1;}
  ::load(pt,path+".Elt_bond_length",Elt_bond_length,is_Elt_bond_length); if(is_Elt_bond_length==1) { status=1;}
  ::load(pt,path+".Elt_radius",Elt_radius,is_Elt_radius); if(is_Elt_radius==1) { status=1;}

}

void load(boost::property_tree::ptree& pt,std::string path,vector<Element>& vt,int& status){
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


