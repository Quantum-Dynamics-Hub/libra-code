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
  \file Thermostat.cpp
  \brief The file implements the Thermostat class for constant-temperature calculations
    
*/

#include "Thermostat.h"
#include "../../io/libio.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

/// libthermostat namespace
namespace libthermostat{


void Thermostat::set(object at){
/** 
  \briefSet properties of the Thermostat object from an arbitrary Python object.

  \param[in] at The input object - must contain the members with the names that match the names of the internal variables.
 
*/

 set_value(is_s_var,   s_var,    at,"s_var");
 set_value(is_Ps,      Ps,       at,"Ps");
 set_value(is_Q,       Q,        at,"Q");
 set_value(is_Nf_t,    Nf_t,     at,"Nf_t");
 set_value(is_Nf_r,    Nf_r,     at,"Nf_r");
 set_value(is_Nf_b,    Nf_b,     at,"Nf_b");
 set_value(is_nu_therm,nu_therm, at,"nu_therm");
 set_value(is_NHC_size,NHC_size, at,"NHC_size");
 set_value(is_Temperature,    Temperature,    at,"Temperature");
 set_value(is_thermostat_type, thermostat_type, at, "thermostat_type");
}

void Thermostat::show_info(){
/** 
  \brief Show info about Thermostat state and properties
*/

  std::cout<<"Thermostat properties:"<<std::endl;
  if(is_s_var)     {std::cout<<"s_var = "<<s_var<<" unitless"<<std::endl;   }
  if(is_Ps)   {std::cout<<"Ps = "<<Ps<<" some units"<<std::endl;   }
  if(is_Q)    {std::cout<<"Q  = "<<Q<<" some units"<<std::endl;   }
  if(is_Nf_t) {std::cout<<"Nf_t = "<<Nf_t<<std::endl; }
  if(is_Nf_r) {std::cout<<"Nf_r = "<<Nf_r<<std::endl; }
  if(is_Nf_b) {std::cout<<"Nf_b = "<<Nf_b<<std::endl; }
  if(is_nu_therm) {std::cout<<"nu_therm = "<<nu_therm<<" tau"<<std::endl; }
  if(is_NHC_size) {std::cout<<"NHC_size = "<<NHC_size<<" unitless"<<std::endl; }
  if(is_Temperature)   {std::cout<<"Temperature = "<<Temperature<<" K"<<std::endl;   }
  if(is_thermostat_type){ std::cout<<"thermostat_type = "<<thermostat_type<<std::endl; }
  std::cout<<std::endl;
}

void Thermostat::extract_dictionary(boost::python::dict d){
/** 
  \brief Set properties of the Thermostat object from the Python dictionary

  \param[in] d The input dictionary - must contain the keys that match the internal variables names
*/

  std::string key;
  for(int i=0;i<len(d.values());i++){
    key = extract<std::string>(d.keys()[i]);

    if(key=="s_var") { s_var = extract<double>(d.values()[i]);  is_s_var = 1; }
    else if(key=="Ps") { Ps = extract<double>(d.values()[i]);  is_Ps = 1; }
    else if(key=="Q") { Q = extract<double>(d.values()[i]);  is_Q = 1; }
    else if(key=="nu_therm") { nu_therm = extract<double>(d.values()[i]);  is_nu_therm = 1; }
    else if(key=="NHC_size") { NHC_size = extract<int>(d.values()[i]);  is_NHC_size = 1; }
    else if(key=="Temperature") { Temperature = extract<double>(d.values()[i]);  is_Temperature = 1; }
    else if(key=="thermostat_type") { thermostat_type = extract<std::string>(d.values()[i]);  is_thermostat_type = 1; }

  }
}

void Thermostat::init_variables(){
/** 
  \brief Initialize Thermostat variables to the default values
*/

  s_var = 1.0;                  is_s_var = 1;
  Ps = 0.0;                     is_Ps = 1;
  Q        = 100.0;             is_Q = 1;
  Nf_t = 1;                     is_Nf_t = 1;
  Nf_r = 0;                     is_Nf_r = 1;
  Nf_b = 0;                     is_Nf_b = 1;
  nu_therm = 1.0;               is_nu_therm = 1;
  NHC_size = 1;                 is_NHC_size = 1;
  Temperature = 300.0;          is_Temperature = 1;
  thermostat_type = "Nose-Poincare";     is_thermostat_type = 1;

  s_t_size = 0;
  s_r_size = 0;
  s_b_size = 0;
  ksi_t_size = 0;
  ksi_r_size = 0;
  ksi_b_size = 0;
  G_t_size = 0;
  G_r_size = 0;
  G_b_size = 0;
  Q_t_size = 0;
  Q_r_size = 0;
  Q_b_size = 0;
}


void Thermostat::copy_content(const Thermostat& th){
/** 
  \brief Copies one thermostat object into the other one. 

  Only the properties that are set in the source object are copied into the target project

  \param[in] th The input Thermostat object
*/

  if(th.is_s_var) { s_var = th.s_var;  is_s_var = 1;}
  if(th.is_Ps) { Ps = th.Ps;  is_Ps = 1;}
  if(th.is_Q)  { Q  = th.Q;   is_Q  = 1;}
  if(th.is_Nf_t){ Nf_t = th.Nf_t; is_Nf_t = 1; }
  if(th.is_Nf_r){ Nf_r = th.Nf_r; is_Nf_r = 1; }
  if(th.is_Nf_b){ Nf_b = th.Nf_b; is_Nf_b = 1; }
  if(th.is_nu_therm) { nu_therm = th.nu_therm;  is_nu_therm = 1; }
  if(th.is_NHC_size) { NHC_size = th.NHC_size;  is_NHC_size = 1; }
  if(th.is_Temperature) { Temperature = th.Temperature;  is_Temperature = 1;}
  if(th.is_thermostat_type){ thermostat_type = th.thermostat_type; is_thermostat_type = 1;}

  if(th.s_t_size>0)   {  s_t = th.s_t;  s_t_size = th.s_t_size; }
  if(th.s_r_size>0)   {  s_r = th.s_r;  s_r_size = th.s_r_size; }
  if(th.s_b_size>0)   {  s_b = th.s_b;  s_b_size = th.s_b_size; }
  if(th.ksi_t_size>0) {  ksi_t = th.ksi_t;  ksi_t_size = th.ksi_t_size; }
  if(th.ksi_r_size>0) {  ksi_r = th.ksi_r;  ksi_r_size = th.ksi_r_size; }
  if(th.ksi_b_size>0) {  ksi_b = th.ksi_b;  ksi_b_size = th.ksi_b_size; }
  if(th.G_t_size>0)   {  G_t = th.G_t;  G_t_size = th.G_t_size; }
  if(th.G_r_size>0)   {  G_r = th.G_r;  G_r_size = th.G_r_size; }
  if(th.G_b_size>0)   {  G_b = th.G_b;  G_b_size = th.G_b_size; }
  if(th.Q_t_size>0)   {  Q_t = th.Q_t;  Q_t_size = th.Q_t_size; }
  if(th.Q_r_size>0)   {  Q_r = th.Q_r;  Q_r_size = th.Q_r_size; }
  if(th.Q_b_size>0)   {  Q_b = th.Q_b;  Q_b_size = th.Q_b_size; }

}



Thermostat::Thermostat(){
/**
  \brief Default Constructor.

  Initializes variables to default values
*/

  init_variables();
}

Thermostat::Thermostat(boost::python::dict d){
/**
  \brief Constructor with Python dictionary argument

  Constructs the Thermostat object using the key:value pairs defined 
  in the input dictionary. The parameters not defined in the dictionary are set to the default values
 
  \param[in] d The Python dictionary. The key names must be consistent with the internal variable (class member) names
*/

  init_variables();
  extract_dictionary(d);
}


Thermostat::Thermostat(const Thermostat& th){
/**
  \brief Copy constructor

  Constructs the Thermostat object from another Thermostat object using only those members that are initialized in the
  input object. The parameters not defined in the input object are set to the default values
 
  \param[in] th The input Thermostat object.
*/

  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(th);
}

Thermostat& Thermostat::operator=(const Thermostat& th){
/**
  \brief Copy (assignment) operator

*/

  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(th);
  return *this;
}


Thermostat::~Thermostat(){
/**
  \brief Destructor

*/

  if(s_t_size>0)   {  s_t.clear(); }
  if(s_r_size>0)   {  s_r.clear(); }
  if(s_b_size>0)   {  s_b.clear(); }
  if(ksi_t_size>0) {  ksi_t.clear(); }
  if(ksi_r_size>0) {  ksi_r.clear(); }
  if(ksi_b_size>0) {  ksi_b.clear(); }
  if(G_t_size>0)   {  G_t.clear(); }
  if(G_r_size>0)   {  G_r.clear(); }
  if(G_b_size>0)   {  G_b.clear(); }
  if(Q_t_size>0)   {  Q_t.clear(); }
  if(Q_r_size>0)   {  Q_r.clear(); }
  if(Q_b_size>0)   {  Q_b.clear(); }

}

void Thermostat::save(boost::property_tree::ptree& pt,std::string path){
/**
  \brief Save the state of the Thermostat object as a property tree

  Each defined data member is added as a node to the property tree. The nodes are added to 
  the level of the tree controlled by the path variable.
 
  \param[in,out] pt The property tree to which the properties of the Thermostat are added
  \param[in] path The parameter controlling the level of the tree to which the Thermostat members will be added.
*/

  if(s_t_size>0){  libio::save(pt,path+".s_t",s_t);    }
  if(s_r_size>0){  libio::save(pt,path+".s_r",s_r);    }
  if(s_b_size>0){  libio::save(pt,path+".s_b",s_b);    }
  if(ksi_t_size>0){  libio::save(pt,path+".ksi_t",ksi_t);    }
  if(ksi_r_size>0){  libio::save(pt,path+".ksi_r",ksi_r);    }
  if(ksi_b_size>0){  libio::save(pt,path+".ksi_b",ksi_b);    }
  if(G_t_size>0){  libio::save(pt,path+".G_t",G_t);    }
  if(G_r_size>0){  libio::save(pt,path+".G_r",G_r);    }
  if(G_b_size>0){  libio::save(pt,path+".G_b",G_b);    }
  if(Q_t_size>0){  libio::save(pt,path+".Q_t",Q_t);    }
  if(Q_r_size>0){  libio::save(pt,path+".Q_r",Q_r);    }
  if(Q_b_size>0){  libio::save(pt,path+".Q_b",Q_b);    }


  if(is_Nf_t){  libio::save(pt,path+".Nf_t",Nf_t);    }
  if(is_Nf_r){  libio::save(pt,path+".Nf_r",Nf_r);    }
  if(is_Nf_b){  libio::save(pt,path+".Nf_b",Nf_b);    }
  if(is_s_var){  libio::save(pt,path+".s_var",s_var);    }
  if(is_Ps){  libio::save(pt,path+".Ps",Ps);    }
  if(is_Q){  libio::save(pt,path+".Q",Q);    }
  if(is_NHC_size){  libio::save(pt,path+".NHC_size",NHC_size);    }
  if(is_nu_therm){  libio::save(pt,path+".nu_therm",nu_therm);    }
  if(is_Temperature){  libio::save(pt,path+".Temperature",Temperature);    }
  if(is_thermostat_type){  libio::save(pt,path+".thermostat_type",thermostat_type);    }

}

void save(boost::property_tree::ptree& pt,std::string path,vector<Thermostat>& vt){
/**
  \brief Save the state of the vector of Thermostat objects as a property tree

  Each Thermostat object is added as a separate branch. 
 
  \param[in,out] pt The property tree to which the list of the Thermostat objects will be added
  \param[in] path The parameter controlling the level of the tree to which the list of Thermostats will be added.
  \param[in] vt The list of Thermostat objects to be printed out into property tree
*/

  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".Thermostat"+rt);
  }
}

void Thermostat::load(boost::property_tree::ptree& pt,std::string path,int& status){
/**
  \brief Load the state of the Thermostat object from a property tree

  Each data member found in the property tree is extracted as the member of the Thermostat object. The
  status of each found data member is set to 1.
 
  \param[in] pt The property tree from which the properties of the Thermostat will be extracted
  \param[in] path The parameter controlling from which level of the tree we try to extract the Thermostat object
  \param[out] status Is the global status of the success of the operation. It is 1 is at least one Thermostat member is found at
              given level of the property tree.
*/

  int st;
  status = 0;

  libio::load(pt,path+".s_t",s_t,st); if(st==1) { status=1; s_t_size = s_t.size(); }
  libio::load(pt,path+".s_r",s_r,st); if(st==1) { status=1; s_r_size = s_r.size(); }
  libio::load(pt,path+".s_b",s_b,st); if(st==1) { status=1; s_b_size = s_b.size(); }
  libio::load(pt,path+".ksi_t",ksi_t,st); if(st==1) { status=1; ksi_t_size = ksi_t.size(); }
  libio::load(pt,path+".ksi_r",ksi_r,st); if(st==1) { status=1; ksi_r_size = ksi_r.size(); }
  libio::load(pt,path+".ksi_b",ksi_b,st); if(st==1) { status=1; ksi_b_size = ksi_b.size(); }
  libio::load(pt,path+".G_t",G_t,st); if(st==1) { status=1; G_t_size = G_t.size(); }
  libio::load(pt,path+".G_r",G_r,st); if(st==1) { status=1; G_r_size = G_r.size(); }
  libio::load(pt,path+".G_b",G_b,st); if(st==1) { status=1; G_b_size = G_b.size(); }
  libio::load(pt,path+".Q_t",Q_t,st); if(st==1) { status=1; Q_t_size = Q_t.size(); }
  libio::load(pt,path+".Q_r",Q_r,st); if(st==1) { status=1; Q_r_size = Q_r.size(); }
  libio::load(pt,path+".Q_b",Q_b,st); if(st==1) { status=1; Q_b_size = Q_b.size(); }

  libio::load(pt,path+".Nf_t",Nf_t,is_Nf_t); if(is_Nf_t==1) { status=1;}
  libio::load(pt,path+".Nf_r",Nf_r,is_Nf_r); if(is_Nf_r==1) { status=1;}
  libio::load(pt,path+".Nf_b",Nf_b,is_Nf_b); if(is_Nf_b==1) { status=1;}
  libio::load(pt,path+".s_var",s_var,is_s_var); if(is_s_var==1) { status=1;}
  libio::load(pt,path+".Ps",Ps,is_Ps); if(is_Ps==1) { status=1;}
  libio::load(pt,path+".Q",Q,is_Q); if(is_Q==1) { status=1;}
  libio::load(pt,path+".NHC_size",NHC_size,is_NHC_size); if(is_NHC_size==1) { status=1;}
  libio::load(pt,path+".nu_therm",nu_therm,is_nu_therm); if(is_nu_therm==1) { status=1;}
  libio::load(pt,path+".Temperature",Temperature,is_Temperature); if(is_Temperature==1) { status=1;}
  libio::load(pt,path+".thermostat_type",thermostat_type,is_thermostat_type); if(is_thermostat_type==1) { status=1;}

}

void load(boost::property_tree::ptree& pt,std::string path,vector<Thermostat>& vt,int& status){
/**
  \brief Load the vector of Thermostat objects from a property tree

  Each Thermostat object is extracted from a separate branch. 
 
  \param[in] pt The property tree from which the vector of Thermostat objects will be extracted
  \param[in] path The parameter controlling from which level of the property tree we will try to extract the vector of Thermostat objects
  \param[out] vt The vector of created Thermostat objects
  \param[out] status Is the global status of the success of the operation. It is 1 is at least one Thermostat object is extracted
*/

  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      Thermostat x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}




}// namespace libthermostat
}// namespace libdyn
}// liblibra



