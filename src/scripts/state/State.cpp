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

#include "State.h"

/// liblibra namespace
namespace liblibra{

namespace libscripts{
namespace libstate{


void MD::init_variables(){

  ensemble = "NVT";              is_ensemble = 1;
  integrator = "DLML";           is_integrator = 1;
  dt = 0.02;                     is_dt = 1;
  n_medium = 1;                  is_n_medium = 1;
  n_fast = 1;                    is_n_fast = 1;
  n_outer = 1;                   is_n_outer = 1;
  max_time = 0.0;                is_max_time = 1;
  max_step = 0;                  is_max_step = 1;
  curr_time = 0.0;               is_curr_time = 1;
  curr_step = 0;                 is_curr_step = 1;

  use_vlist = 0;                 is_use_vlist = 1;
  vlist_upd_freq = 1;            is_vlist_upd_freq = 1;
  vlist_time = 0;                is_vlist_time = 1;

  terec_exp_size = 10;           is_terec_exp_size = 1;

}

void MD::copy_content(const MD& md){

  if(md.is_ensemble) { ensemble = md.ensemble;  is_ensemble = 1;}
  if(md.is_integrator) { integrator = md.integrator;  is_integrator = 1;}
  if(md.is_dt) { dt = md.dt;  is_dt = 1;}
  if(md.is_n_medium) { n_medium = md.n_medium; is_n_medium = 1; }
  if(md.is_n_fast) { n_fast = md.n_fast; is_n_fast = 1; }
  if(md.is_n_outer){ n_outer = md.n_outer; is_n_outer = 1; }
  if(md.is_max_time) { max_time = md.max_time;  is_max_time = 1;}
  if(md.is_max_step) { max_step = md.max_step;  is_max_step = 1;}
  if(md.is_curr_time) { curr_time = md.curr_time;  is_curr_time = 1;}
  if(md.is_curr_step) { curr_step = md.curr_step;  is_curr_step = 1;}

  if(md.is_use_vlist){ use_vlist = md.use_vlist; is_use_vlist = 1; }
  if(md.is_vlist_upd_freq) { vlist_upd_freq = md.vlist_upd_freq; is_vlist_upd_freq = 1; }
  if(md.is_vlist_time) { vlist_time = md.vlist_time; is_vlist_time = 1; }
  
  if(md.is_terec_exp_size){ terec_exp_size = md.terec_exp_size; is_terec_exp_size = 1; }
}

void MD::extract_dictionary(boost::python::dict d){
  std::string key;
  for(int i=0;i<len(d.values());i++){
    key = extract<std::string>(d.keys()[i]);

    if(key=="ensemble") { ensemble = extract<std::string>(d.values()[i]);  is_ensemble = 1; }
    else if(key=="integrator") { integrator = extract<std::string>(d.values()[i]);  is_integrator = 1; }
    else if(key=="dt") { dt = extract<double>(d.values()[i]);  is_dt = 1; }
    else if(key=="n_medium") { n_medium = extract<int>(d.values()[i]);  is_n_medium = 1; }
    else if(key=="n_fast") { n_fast = extract<int>(d.values()[i]);  is_n_fast = 1; }
    else if(key=="n_outer") { n_outer = extract<int>(d.values()[i]); is_n_outer = 1; }
    else if(key=="max_time") { max_time = extract<double>(d.values()[i]);  is_max_time = 1; }
    else if(key=="max_step") { max_step = extract<int>(d.values()[i]);  is_max_step = 1; }
    else if(key=="curr_time") { curr_time = extract<double>(d.values()[i]);  is_curr_time = 1; }
    else if(key=="curr_step") { curr_step = extract<int>(d.values()[i]);  is_curr_step = 1; }
    else if(key=="use_vlist") { use_vlist = extract<int>(d.values()[i]);  is_use_vlist = 1; }
    else if(key=="terec_exp_size"){ terec_exp_size = extract<int>(d.values()[i]); is_terec_exp_size = 1; }
  }

  if(!is_max_time){ if(is_max_step && is_dt){ is_max_time = dt*max_step; is_max_time = 1; } }

}


MD::MD(){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
}

MD::MD(boost::python::dict d){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
  extract_dictionary(d);
}

MD::MD(const MD& md){
  /********************
    Copy constructor
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(md);
}

MD& MD::operator=(const MD& md){
  /********************
    Assignment operator
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(md);
  return *this;
}


MD::~MD(){ }


/*
void MD::set(object at){
  set_value(is_ensemble,  ensemble,    at,"ensemble");
  set_value(is_integrator, integrator,   at,"integrator");
  set_value(is_dt,    dt,    at,"dt");
  set_value(is_max_time,    max_time,    at,"max_time");
  set_value(is_max_step,    max_step,    at,"max_step");
  set_value(is_curr_time,   curr_time,   at,"curr_time");
  set_value(is_curr_step,   curr_step,   at,"curr_step");


  set_value(is_use_vlist,use_vlist,at,"use_vlist");
  set_value(is_vlist_upd_freq, vlist_upd_freq, at,"vlist_upd_freq");

}
*/

void MD::show_info(){
  std::cout<<"MD properties:"<<std::endl;
  if(is_ensemble) {std::cout<<"ensemble = "<<ensemble<<" unitless"<<std::endl;   }
  if(is_integrator)   {std::cout<<"integrator = "<<integrator<<" unitless"<<std::endl;   }
  if(is_dt)   {std::cout<<"dt = "<<dt<<" tau"<<std::endl;   }
  if(is_n_medium)   {std::cout<<"n_medium = "<<n_medium<<" steps"<<std::endl;   }
  if(is_n_fast)   {std::cout<<"n_fast = "<<n_fast<<" steps"<<std::endl;   }
  if(is_n_outer) { std::cout<<"n_outer = "<<n_outer<<" steps"<<std::endl; }
  if(is_max_time) {std::cout<<"max_time = "<<max_time<<" fs"<<std::endl;   }
  if(is_max_step) {std::cout<<"max_step = "<<max_step<<" unitless"<<std::endl;   }
  if(is_curr_time) {std::cout<<"curr_time = "<<curr_time<<" fs"<<std::endl;   }
  if(is_curr_step) {std::cout<<"curr_step = "<<curr_step<<" unitless"<<std::endl;   }

  if(is_use_vlist){ std::cout<<"use_vlist = "<<use_vlist<<std::endl; }
  if(is_vlist_upd_freq){ std::cout<<"vlist_upd_freq = "<<vlist_upd_freq<<" steps"<<std::endl; }
  if(is_terec_exp_size){ std::cout<<"terec_exp_size = "<<terec_exp_size<<std::endl; }
  std::cout<<std::endl;
}

void MD::save(boost::property_tree::ptree& pt,std::string path){

  if(is_ensemble){  libio::save(pt,path+".ensemble",ensemble);    }
  if(is_integrator){  libio::save(pt,path+".integrator",integrator);    }
  if(is_dt){  libio::save(pt,path+".dt",dt);    }
  if(is_n_medium){  libio::save(pt,path+".n_medium",n_medium);    }
  if(is_n_fast){  libio::save(pt,path+".n_fast",n_fast);    }
  if(is_n_outer){ libio::save(pt,path+".n_outer",n_outer); }
  if(is_max_time){  libio::save(pt,path+".max_time",max_time);    }
  if(is_max_step){  libio::save(pt,path+".max_step",max_step);    }
  if(is_curr_time){  libio::save(pt,path+".curr_time",curr_time);    }
  if(is_curr_step){  libio::save(pt,path+".curr_step",curr_step);    }
  if(is_terec_exp_size){  libio::save(pt,path+".terec_exp_size",terec_exp_size);    }
  if(is_use_vlist){  libio::save(pt,path+".use_vlist",use_vlist);    }
  if(is_vlist_upd_freq){  libio::save(pt,path+".vlist_upd_freq",vlist_upd_freq);    }
  if(is_vlist_time){  libio::save(pt,path+".vlist_time",vlist_time);    }


}
 
void save(boost::property_tree::ptree& pt,std::string path,vector<MD>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".MD"+rt);
  }
}

void MD::load(boost::property_tree::ptree& pt,std::string path,int& status){
  int st;
  status = 0;

  libio::load(pt,path+".ensemble",ensemble,is_ensemble); if(is_ensemble==1) { status=1;}
  libio::load(pt,path+".integrator",integrator,is_integrator); if(is_integrator==1) { status=1;}
  libio::load(pt,path+".dt",dt,is_dt); if(is_dt==1) { status=1;}
  libio::load(pt,path+".n_medium",n_medium,is_n_medium); if(is_n_medium==1) { status=1;}
  libio::load(pt,path+".n_fast",n_fast,is_n_fast); if(is_n_fast==1) { status=1;}
  libio::load(pt,path+".n_outer",n_outer,is_n_outer); if(is_n_outer==1){ status = 1; }
  libio::load(pt,path+".max_time",max_time,is_max_time); if(is_max_time==1) { status=1;}
  libio::load(pt,path+".max_step",max_step,is_max_step); if(is_max_step==1) { status=1;}
  libio::load(pt,path+".curr_time",curr_time,is_curr_time); if(is_curr_time==1) { status=1;}
  libio::load(pt,path+".curr_step",curr_step,is_curr_step); if(is_curr_step==1) { status=1;}
  libio::load(pt,path+".terec_exp_size",terec_exp_size,is_terec_exp_size); if(is_terec_exp_size==1) { status=1;}
  libio::load(pt,path+".use_vlist",use_vlist,is_use_vlist); if(is_use_vlist==1) { status=1;}
  libio::load(pt,path+".vlist_upd_freq",vlist_upd_freq,is_vlist_upd_freq); if(is_vlist_upd_freq==1) { status=1;}
  libio::load(pt,path+".vlist_time",vlist_time,is_vlist_time); if(is_vlist_time==1) { status=1;}

}

void load(boost::property_tree::ptree& pt,std::string path,vector<MD>& vt,int& status){
  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      MD x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}


//====================== STATE ================================


void State::init_variables(){

  // Objects
  syst = NULL;
  thermostat = NULL;
  barostat = NULL;

  // Simulation parameters
  md = NULL;

  // State variables
  E_kin = 0.0;    is_E_kin = 1;
  E_kin_tr = 0.0; is_E_kin_tr = 1;
  E_kin_rot = 0.0;is_E_kin_rot = 1;
  E_pot = 0.0;    is_E_pot = 1;
  E_tot = 0.0;    is_E_tot = 1;
  H = 0.0;        is_H = 1;
  H_NHC = 0.0;    is_H_NHC = 1;
  H_NP  = 0.0;    is_H_NP  = 1;
  is_H0    = 0;
  curr_T = 0.0;   is_curr_T = 1;
  curr_V = 0.0;   is_curr_V = 1;
  curr_P = 0.0;   is_curr_P = 1;
  curr_P_tens = 0.0; is_curr_P_tens = 1;
  Nf_t = 0;       is_Nf_t = 1;
  Nf_r = 0;       is_Nf_r = 1;
  L_tot = 0.0;    is_L_tot = 1;
  P_tot = 0.0;    is_P_tot = 1;

  is_md_initialized = 0;
  
}

void State::copy_content(const State& st){

  if(st.syst != NULL){ syst = new System; *syst = *st.syst; }
  if(st.thermostat != NULL){ thermostat = new Thermostat; *thermostat = *st.thermostat; }
  if(st.barostat != NULL){ barostat = new Barostat; *barostat = *st.barostat; }

  if(st.md != NULL){ md = new MD; *md = *st.md; }

  if(st.is_E_kin){ E_kin = st.E_kin; is_E_kin = 1; }
  if(st.is_E_kin_tr){ E_kin_tr = st.E_kin_tr; is_E_kin_tr = 1; }
  if(st.is_E_kin_rot){ E_kin_rot = st.E_kin_rot; is_E_kin_rot = 1; }
  if(st.is_E_pot){ E_pot = st.E_pot; is_E_pot = 1; }
  if(st.is_E_tot){ E_tot = st.E_tot; is_E_tot = 1; }
  if(st.is_H){ H = st.H; is_H = 1; }
  if(st.is_H_NHC){ H_NHC = st.H_NHC; is_H_NHC = 1; }
  if(st.is_H_NP){ H_NP = st.H_NP; is_H_NP = 1; }
  if(st.is_H0){ H0 = st.H0; is_H0 = 1; }
  if(st.is_curr_T){ curr_T = st.curr_T; is_curr_T = 1; }
  if(st.is_curr_V){ curr_V = st.curr_V; is_curr_V = 1; }
  if(st.is_curr_P){ curr_P = st.curr_P; is_curr_P = 1; }
  if(st.is_curr_P_tens){ curr_P_tens = st.curr_P_tens; is_curr_P_tens = 1; }
  if(st.is_Nf_t){ Nf_t = st.Nf_t; is_Nf_t = 1; }
  if(st.is_Nf_r){ Nf_r = st.Nf_r; is_Nf_r = 1; }
  if(st.is_L_tot){ L_tot = st.L_tot; is_L_tot = 1; }
  if(st.is_P_tot){ P_tot = st.P_tot; is_P_tot = 1; }

  is_md_initialized = st.is_md_initialized;
}

void State::show_info(){ }

void State::set(object){ }

State::State(){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
}


State::State(const State& st){
  /********************
    Copy constructor
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(st);
}

State& State::operator=(const State& st){
  /********************
    Assignment operator
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(st);
  return *this;
}

State::~State(){  }



}// namespace libstate
}// namespace libscripts
}// liblibra



