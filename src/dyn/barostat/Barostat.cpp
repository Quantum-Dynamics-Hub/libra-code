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

#include "Barostat.h"

namespace libdyn{
namespace libbarostat{


void Barostat::set(object at){
 set_value(is_Wg,      Wg,       at,"Wg");
 set_value(is_Nf_t,    Nf_t,     at,"Nf_t");
 set_value(is_Nf_r,    Nf_r,     at,"Nf_r");
 set_value(is_Nf_b,    Nf_b,     at,"Nf_b");
 set_value(is_nu_baro,nu_baro, at,"nu_baro");
 set_value(is_Pressure,    Pressure,    at,"Pressure");
 set_value(is_barostat_type, barostat_type, at, "barostat_type");
}

void Barostat::show_info(){

  std::cout<<"Barostat properties:"<<std::endl;
  if(is_Wg)   {std::cout<<"Wg  = "<<Wg<<" some units"<<std::endl;   }
  if(is_Nf_t) {std::cout<<"Nf_t = "<<Nf_t<<std::endl; }
  if(is_Nf_r) {std::cout<<"Nf_r = "<<Nf_r<<std::endl; }
  if(is_Nf_b) {std::cout<<"Nf_b = "<<Nf_b<<std::endl; }
  if(is_nu_baro)  {std::cout<<"nu_baro = "<<nu_baro<<" tau"<<std::endl; }
  if(is_Pressure) {std::cout<<"Pressure = "<<(Pressure/atm_to_int)<<" atmospheres "<<Pressure<<" internal units"<<std::endl;   }
  if(is_barostat_type){ std::cout<<"barostat_type = "<<barostat_type<<std::endl; }
  std::cout<<std::endl;
}

void Barostat::init_variables(){
  ksi_eps = 0.0;                is_ksi_eps = 1;
  G_eps = 0.0;                  is_G_eps = 1;      
  eps_iso = 0.0;                is_eps_iso = 1;  
  ksi_eps_iso = 0.0;            is_ksi_eps_iso = 1;  
  G_eps_iso = 0.0;              is_G_eps_iso = 1; 

  Wg       = 100.0;             is_Wg = 1;
  Nf_t = 1;                     is_Nf_t = 1;
  Nf_r = 0;                     is_Nf_r = 1;
  Nf_b = 1;                     is_Nf_b = 1;      // This defines isotropic case (uniform dilation)
  nu_baro = 1.0;                is_nu_baro = 1;
  Pressure = 1.0*atm_to_int;    is_Pressure = 1;  // Default: 1 atm
  barostat_type = "None";       is_barostat_type = 1;
}

void Barostat::copy_content(const Barostat& bar){

  if(bar.is_ksi_eps){ ksi_eps = bar.ksi_eps; is_ksi_eps = 1; }
  if(bar.is_G_eps){ G_eps = bar.G_eps; is_G_eps = 1; }
  if(bar.is_eps_iso){ eps_iso = bar.eps_iso; is_eps_iso = 1; }
  if(bar.is_ksi_eps_iso){ ksi_eps_iso = bar.ksi_eps_iso; is_ksi_eps_iso = 1; }
  if(bar.is_G_eps_iso){ G_eps_iso = bar.G_eps_iso; is_G_eps_iso = 1; }
  if(bar.is_Wg)  { Wg  = bar.Wg;   is_Wg  = 1;}
  if(bar.is_Nf_t){ Nf_t = bar.Nf_t; is_Nf_t = 1; }
  if(bar.is_Nf_r){ Nf_r = bar.Nf_r; is_Nf_r = 1; }
  if(bar.is_Nf_b){ Nf_b = bar.Nf_b; is_Nf_b = 1; }
  if(bar.is_nu_baro) { nu_baro = bar.nu_baro;  is_nu_baro = 1; }
  if(bar.is_Pressure) { Pressure = bar.Pressure;  is_Pressure = 1;}
  if(bar.is_barostat_type){ barostat_type = bar.barostat_type; is_barostat_type = 1;}

}

void Barostat::extract_dictionary(boost::python::dict d){
  std::string key;
  for(int i=0;i<len(d.values());i++){
    key = extract<std::string>(d.keys()[i]);
    if(key=="Wg") { Wg = extract<double>(d.values()[i]);  is_Wg = 1; }
    else if(key=="nu_baro") { nu_baro = extract<double>(d.values()[i]);  is_nu_baro = 1; }
    else if(key=="Pressure") { Pressure = extract<double>(d.values()[i]); 
                               Pressure *= atm_to_int;  // Convert atmospheres to internal units
                               is_Pressure = 1; }
    else if(key=="barostat_type") { barostat_type = extract<double>(d.values()[i]);  is_barostat_type = 1; }
  }
}

Barostat::Barostat(){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
}

Barostat::Barostat(boost::python::dict d){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
  extract_dictionary(d);
}

Barostat::Barostat(const Barostat& bar){
  /********************
    Copy constructor
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of bar object which is defined
  copy_content(bar);
}

Barostat& Barostat::operator=(const Barostat& bar){
  /********************
    Assignment operator
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of bar object which is defined
  copy_content(bar);
  return *this;
}


Barostat::~Barostat(){ ;; }


void Barostat::save(boost::property_tree::ptree& pt,std::string path){

  if(is_Nf_t){  ::save(pt,path+".Nf_t",Nf_t);    }
  if(is_Nf_r){  ::save(pt,path+".Nf_r",Nf_r);    }
  if(is_Nf_b){  ::save(pt,path+".Nf_b",Nf_b);    }
  if(is_ksi_eps){  ::save(pt,path+".ksi_eps",ksi_eps);    }
  if(is_G_eps){  ::save(pt,path+".G_eps",G_eps);    }
  if(is_eps_iso){  ::save(pt,path+".eps_iso",eps_iso);    }
  if(is_ksi_eps_iso){  ::save(pt,path+".ksi_eps_iso",ksi_eps_iso);    }
  if(is_G_eps_iso){  ::save(pt,path+".G_eps_iso",G_eps_iso);    }
  if(is_Wg){  ::save(pt,path+".Wg",Wg);    }
  if(is_nu_baro){  ::save(pt,path+".nu_baro",nu_baro);    }
  if(is_Pressure){  ::save(pt,path+".Pressure",Pressure);    }
  if(is_barostat_type){  ::save(pt,path+".barostat_type",barostat_type);    }


}

void save(boost::property_tree::ptree& pt,std::string path,vector<Barostat>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".Barostat"+rt);
  }
}

 
void Barostat::load(boost::property_tree::ptree& pt,std::string path,int& status){

  int st;
  status = 0;

  ::load(pt,path+".Nf_t",Nf_t,is_Nf_t); if(is_Nf_t==1) { status=1;}
  ::load(pt,path+".Nf_r",Nf_r,is_Nf_r); if(is_Nf_r==1) { status=1;}
  ::load(pt,path+".Nf_b",Nf_b,is_Nf_b); if(is_Nf_b==1) { status=1;}
  ::load(pt,path+".ksi_eps",ksi_eps,is_ksi_eps); if(is_ksi_eps==1) { status=1;}
  ::load(pt,path+".G_eps",G_eps,is_G_eps); if(is_G_eps==1) { status=1;}
  ::load(pt,path+".eps_iso",eps_iso,is_eps_iso); if(is_eps_iso==1) { status=1;}
  ::load(pt,path+".ksi_eps_iso",ksi_eps_iso,is_ksi_eps_iso); if(is_ksi_eps_iso==1) { status=1;}
  ::load(pt,path+".G_eps_iso",G_eps_iso,is_G_eps_iso); if(is_G_eps_iso==1) { status=1;}
  ::load(pt,path+".Wg",Wg,is_Wg); if(is_Wg==1) { status=1;}
  ::load(pt,path+".nu_baro",nu_baro,is_nu_baro); if(is_nu_baro==1) { status=1;}
  ::load(pt,path+".Pressure",Pressure,is_Pressure); if(is_Pressure==1) { status=1;}
  ::load(pt,path+".barostat_type",barostat_type,is_barostat_type); if(is_barostat_type==1) { status=1;}


}

void load(boost::property_tree::ptree& pt,std::string path,vector<Barostat>& vt,int& status){
  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      Barostat x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}



}// namespace libbarostat
}// namespace libdyn




