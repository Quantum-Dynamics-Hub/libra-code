#include "../State.h"

void State::init_variables(){

  // Objects
  system = NULL;
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

  if(st.system != NULL){ system = new System; *system = *st.system; }
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

State::State(System& sys){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
  system = &sys;
  update();
}

State::State(System& sys,Thermostat& th){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
  system = &sys;
  thermostat = &th;
  update();
}

State::State(System& sys,Barostat& bar){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
  system = &sys;
  barostat = &bar;
  update();
}

State::State(System& sys,Thermostat& th,Barostat& bar){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
  system = &sys;
  thermostat = &th;
  barostat = &bar;
  update();
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

