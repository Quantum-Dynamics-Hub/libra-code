#include "../State.h"

void State::update(){
/********************************************
 This function updated state valiables using 
 corresponding parameters and/or data from its
 sub-objects (system,thermostat,etc.)
*********************************************/

 if(system!=NULL){
   curr_V = system->volume(); is_curr_V = 1;
 }
}

void State::cool(){
  //------------- System -------------------
  if(system!=NULL){
  for(int i=0;i<system->Number_of_fragments;i++){
    system->Fragments[i].Group_RB.scale_angular_(0.0);
    system->Fragments[i].Group_RB.scale_linear_(0.0);
  }
  system->zero_forces_and_torques();
  }
  if(thermostat!=NULL){
    if(thermostat->thermostat_type=="Nose-Poincare"){
      thermostat->Ps = 0.0;     thermostat->is_Ps = 1;
      thermostat->s_var = 1.0;  thermostat->is_s_var = 1;   
    }
    else if(thermostat->thermostat_type=="Nose-Hoover"){
      thermostat->init_nhc();
    }
  }
  if(barostat!=NULL){
    barostat->ksi_eps = 0.0;  barostat->is_ksi_eps = 1;
    barostat->ksi_eps_iso = 0.0;  barostat->is_ksi_eps_iso = 1;
  }
  //------------ State variables -----------
  E_kin = 0.0;
  E_pot = 0.0;
  H0 = E_pot; is_H0 = 0;
}

void State::init_velocities(double Temp){
  VECTOR TOT_P,TOT_L; TOT_P = 0.0; TOT_L = 0.0;
  init_velocities(Temp,TOT_P,TOT_L);
}
void State::init_velocities(double Temp,VECTOR TOT_P,VECTOR TOT_L){
  // This approach sets linear and angular momenta independently
  // and makes sure that the total linear and angular momenta are
  // proportional to corresponding arguments (or equal if this is consistent)
  // with the Temperature parameters
  // In any case the temperature is enforced to be constrained, while the
  // amplitude of the total momenta may be scaled if needed

  L_tot = 0.0; is_L_tot = 1;
  P_tot = 0.0; is_P_tot = 1;

  // Initialize random number generator
  srand( (unsigned)time(NULL) );
  int i;

  VECTOR* temp_p;
  VECTOR* temp_l;
  temp_p = new VECTOR[system->Number_of_fragments];
  temp_l = new VECTOR[system->Number_of_fragments];

  VECTOR tot_p; tot_p = 0.0;
  VECTOR tot_l; tot_l = 0.0;

  for(i=0;i<system->Number_of_fragments;i++){
    temp_p[i].x  = RANDOM(-0.5,0.5);
    temp_p[i].y  = RANDOM(-0.5,0.5);
    temp_p[i].z  = RANDOM(-0.5,0.5);
    temp_l[i].x  = RANDOM(-0.5,0.5);
    temp_l[i].y  = RANDOM(-0.5,0.5);
    temp_l[i].z  = RANDOM(-0.5,0.5);
    tot_p += temp_p[i];
  }

  double size = system->Number_of_fragments;
  tot_p = (TOT_P - tot_p)/size;
  for(i=0;i<system->Number_of_fragments;i++){   temp_p[i] +=  tot_p; }

  // Required temperature value scaling
  double temp_tr  = 0.0;
  double temp_rot = 0.0;
  for(i=0;i<system->Number_of_fragments;i++){
    RigidBody& top = system->Fragments[i].Group_RB;
    top.set_momentum(temp_p[i]);
    temp_tr += top.ekin_tr();
  }// for i

  // Rescale linear momenta to satisfy the translational kinetic energy
  double target_ekin_tr = 0.5*((double)(system->Nf_t))*(boltzmann*Temp);
  double scaling_factor_tr = (temp_tr==0.0)?0.0:sqrt(target_ekin_tr/temp_tr);
//  cout<<"scaling factor tr = "<<scaling_factor_tr<<endl;
  for(i=0;i<system->Number_of_fragments;i++){  temp_p[i] *= scaling_factor_tr; }


  // Angular momenta
  for(i=0;i<system->Number_of_fragments;i++){
    RigidBody& top = system->Fragments[i].Group_RB;
    VECTOR tmp; tmp.cross(top.rb_cm,temp_p[i]);
    tot_l += tmp;
  }
  tot_l = (TOT_L - tot_l)/size;
  for(i=0;i<system->Number_of_fragments;i++){  temp_l[i] =  system->Fragments[i].Group_RB.rb_A_I_to_e * tot_l;   }


  temp_tr = 0.0;
  for(i=0;i<system->Number_of_fragments;i++){
    RigidBody& top = system->Fragments[i].Group_RB;
    top.set_angular_momentum(temp_l[i]);
    top.set_momentum(temp_p[i]);
    temp_rot += top.ekin_rot();
    temp_tr += top.ekin_tr();
  }

//  double target_ekin_rot = 0.5*((double)(system->Nf_r))*(boltzmann*Temp);
//  double scaling_factor_rot = (temp_rot==0.0)?0.0:sqrt(target_ekin_rot/temp_rot);
//  cout<<"scaling factor rot = "<<scaling_factor_rot<<endl;
  double target_ekin = 0.5*((double)(system->Nf_r+system->Nf_t))*(boltzmann*Temp);
  double scaling_factor = (temp_rot==0.0)?0.0:sqrt(target_ekin/(temp_rot+temp_tr));


  for(i=0;i<system->Number_of_fragments;i++){
      temp_p[i] *= scaling_factor;
      temp_l[i] *= scaling_factor;
  }

  // Set corresponding variables
  E_kin = 0.0;
  L_tot = 0.0;
  P_tot = 0.0;
  for(i=0;i<system->Number_of_fragments;i++){
    RigidBody& top = system->Fragments[i].Group_RB; 
    top.set_momentum(temp_p[i]);
    top.set_angular_momentum(temp_l[i]);
    E_kin += (top.ekin_tr() + top.ekin_rot());

    VECTOR tmp; tmp.cross(top.rb_cm,top.rb_p);
    L_tot += top.rb_A_I_to_e_T * top.rb_l_e + tmp;
    P_tot += top.rb_p;

  }// for i - all fragments
  
//  cout<<"E_kin = "<<E_kin<<endl;
  curr_T =  2.0*E_kin/(((double)(system->Nf_t + system->Nf_r))*boltzmann);
  cout<<"in init_velocities...\n";
  cout<<"P_tot = "<<P_tot<<endl;
  cout<<"L_tot = "<<L_tot<<endl;
  cout<<"cutt_T = "<<curr_T<<endl;
  delete [] temp_p;
  delete [] temp_l;

}


