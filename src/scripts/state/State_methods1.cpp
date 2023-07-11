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


void State::init_md(Nuclear& mol, Electronic& el, Hamiltonian& ham, Random& rnd){

  //----- Checking type of ensemble and the existence of required variables ------
  int is_thermostat, is_barostat, is_system;
  is_system = (syst!=NULL);
  is_thermostat = (thermostat!=NULL);
  is_barostat   = (barostat!=NULL);
  cout<<"is_system = "<<is_system<<endl;
  cout<<"is_thermostat = "<<is_thermostat<<endl;
  cout<<"is_barostat = "<<is_barostat<<endl;

  int Nf_b = 0;
  if(!is_system){ std::cout<<"Fatal error: No system defined - nothing to simulate\n Exiting now...\n"; exit(1); }

  if(md->ensemble=="NVT"||md->ensemble=="NPT"){
    if(!is_thermostat){ std::cout<<"Error: No Thermostat object were set up. Can not run simulations in "<<md->ensemble<<" ensemble\n"; exit(1);}    
  }
  if(md->ensemble=="NPT"||md->ensemble=="NPH"){
    Nf_b = 1;
    if(!is_barostat){std::cout<<"Error: No Barostat object were set up. Can not run simulations in "<<md->ensemble<<" ensemble\n"; exit(1);}
  }
  if(md->ensemble=="NPT_FLEX"||md->ensemble=="NPH_FLEX"){
    Nf_b = 9;
    if(!is_barostat){std::cout<<"Error: No Barostat object were set up. Can not run simulations in "<<md->ensemble<<" ensemble\n"; exit(1);}
  }
//  if(md->ensemble=="NVE"){is_thermostat = 0; is_barostat = 0; }

  //---------------- Checking combinations -----------------------------------
  if(is_barostat&&is_thermostat){
    if(thermostat->thermostat_type!="Nose-Hoover" && md->ensemble=="NPT"){  
      std::cout<<"Error: Can not use thermostat other than Nose-Hoover for simulations in NPT ensemble\n Exiting now...\n"; exit(1);
    }
  }

  //------ Initialize thermostat and/or barostat if they exist ----------------
  if(is_thermostat){
    thermostat->set_Nf_t(syst->Nf_t); 
    thermostat->set_Nf_r(syst->Nf_r);
    if(is_barostat){  thermostat->set_Nf_b(Nf_b); }
    thermostat->init_nhc();
    syst->init_fragment_velocities(thermostat->Temperature, rnd);
  }
  if(is_barostat){   
    barostat->set_Nf_t(syst->Nf_t);
    barostat->set_Nf_r(syst->Nf_r);
    barostat->set_Nf_b(Nf_b);
    if(is_thermostat){    barostat->init(thermostat->Temperature); }
  }

  //---------------- Initialize system --------------------
  E_kin = 0.0;
  for(int i=0;i<syst->Number_of_fragments;i++){
    RigidBody& top = syst->Fragments[i].Group_RB;
    E_kin += (top.ekin_tr() + top.ekin_rot());
  }
  syst->zero_forces_and_torques();

//!!!!!!!!!!  E_pot = system->energy(); !!!!!!!!!!!!!!!!!
  syst->extract_atomic_q(mol.q); // syst -> mol
  E_pot = compute_potential_energy(mol, el, ham, 1); //  # 1 - FSSH forces
  compute_forces(mol, el, ham, 1); 
  syst->set_atomic_f(mol.f);    // mol -> syst

  syst->update_fragment_forces_and_torques();

  E_tot = E_pot + E_kin;
  H0 = E_tot; 
  H_NP = 0.0;
  if(is_thermostat){  H0 += thermostat->energy();  }
 
  if(md->integrator=="Terec"||md->integrator=="qTerec"){
    for(int i=0;i<syst->Number_of_fragments;i++){
      syst->Fragments[i].Group_RB.initialize_terec(md->terec_exp_size);
    }
  }
  is_md_initialized = 1;

  if(is_thermostat){ thermostat->show_info(); }
  if(is_barostat){  barostat->show_info(); }

}

void State::run_md(Nuclear& mol, Electronic& el, Hamiltonian& ham){
  int i;
  if(md==NULL) { std::cout<<"Error: MD parameters have not been defined\n"; exit(1);}
  if(!is_md_initialized){    std::cout<<"Error: Need to call init_md() first. MD is not initialized\n"; exit(2);   }

  int is_thermostat, is_barostat;
  is_thermostat = ((thermostat!=NULL) && ((md->ensemble=="NVT")||(md->ensemble=="NPT")||(md->ensemble=="NPT_FLEX")));
  is_barostat   = ((barostat!=NULL) && ((md->ensemble=="NPT")||(md->ensemble=="NPT_FLEX")||(md->ensemble=="NPH")||(md->ensemble=="NPH_FLEX")));

  double dt = md->dt;
  double dt_half = 0.5*md->dt;
  double Nf = syst->Nf_t + syst->Nf_r;
  int Nf_b = 0;
  if(is_barostat) {Nf_b = barostat->get_Nf_b();}
  double scl,sc3,sc4,ksi_r;
  MATRIX3x3 S,I,sc1,sc2;


  while(md->curr_step<md->max_step){

    // Operator NHCB(dt/2)
    if(is_thermostat){
      double ekin_baro = 0.0;
      if(is_barostat){  ekin_baro = barostat->ekin_baro(); }
      thermostat->update_thermostat_forces(syst->ekin_tr(),syst->ekin_rot(),ekin_baro);
      thermostat->propagate_nhc(dt_half,syst->ekin_tr(),syst->ekin_rot(),ekin_baro);
    }


    if(is_thermostat){  thermostat->propagate_sPs(dt_half);    }

    //bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
    // Operator B(dt/2)
    if(is_barostat){
      if(md->ensemble=="NPT"||md->ensemble=="NPH"){ barostat->update_barostat_forces(syst->ekin_tr(),syst->ekin_rot(),curr_V,curr_P);   }
      else if(md->ensemble=="NPT_FLEX"||md->ensemble=="NPH_FLEX"){ barostat->update_barostat_forces(syst->ekin_tr(),syst->ekin_rot(),curr_V,curr_P_tens);   }
      scl = 0.0; if(is_thermostat){ scl = thermostat->get_ksi_b();  }
      barostat->propagate_velocity(dt_half,scl);
    }
    //bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb


    double s_var,Ps,dt_half_s,dt_over_s,dt_over_s2;
    s_var = 1.0; Ps = 0.0;
    dt_half_s = dt_half;
    dt_over_s = dt;
    dt_over_s2 = dt;

    if(md->ensemble=="NVT"||md->ensemble=="NPT"||md->ensemble=="NPT_FLEX"||md->ensemble=="NPH"||md->ensemble=="NPH_FLEX"){
      if(is_thermostat){
        s_var = thermostat->get_s_var(); 
        dt_half_s = dt_half*s_var;
        dt_over_s = (dt/s_var);
        dt_over_s2 = (dt_over_s/s_var);
      }
    }

    //aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa    
    // Operator A(dt/2)
    //-------------------- Linear momentum propagation --------------------
    S = 0.0; I.identity();
    if(is_barostat){
      if(Nf_b==9){ S = barostat->ksi_eps + (barostat->ksi_eps.tr()/(barostat->get_Nf_t()/*+barostat->get_Nf_r()*/))*I; }
      else if(Nf_b==1){S = barostat->ksi_eps_iso * I + (3.0*barostat->ksi_eps_iso/(barostat->get_Nf_t()/*+barostat->get_Nf_r()*/))*I; }
    }
    if(is_thermostat){   S = S + thermostat->get_ksi_t() * I;      }
    sc1 = (exp_(S,-dt_half));//.symmetrized();
    sc2 = dt_half*(exp1_(S,-dt_half*0.5));//.symmetrized()*dt_half;

    //------------------- Angular momentum propagation -----------------------    
    if(is_thermostat){ ksi_r = thermostat->get_ksi_r();}else{ ksi_r = 0.0;}
    sc3 = exp(-dt_half*ksi_r);
    sc4 = dt_half*exp(-0.5*dt_half*ksi_r)*sinh_(0.5*dt_half*ksi_r);



    for(i=0;i<syst->Number_of_fragments;i++){
      RigidBody& top = syst->Fragments[i].Group_RB;
      //-------------------- Linear momentum propagation --------------------
      top.scale_linear_(sc1);
      top.apply_force(sc2);
      //------------------- Angular momentum propagation -----------------------
      top.scale_angular_(sc3);
      top.apply_torque(sc4);       
    }// for i
    

//    if(is_thermostat){
//      thermostat->propagate_Ps(-dt_half*E_pot);
//    }
    //aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

    if(is_thermostat){  thermostat->propagate_Ps(-dt_half*E_pot);    }

    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    //ccccccccccccccccccccccccccccccc Core part ccccccccccccccccccccccccccccccccccccc
    sc1.identity();
    sc2.identity();
    sc2 = sc2 * dt;
    if(is_barostat){
      sc1 = (barostat->pos_scale(dt));
      sc2 = dt*barostat->vpos_scale(dt);
    }

    for(i=0;i<syst->Number_of_fragments;i++){
      RigidBody& top = syst->Fragments[i].Group_RB;
      if(is_thermostat){  thermostat->propagate_Ps( 0.5*dt_over_s2*(top.ekin_rot()+top.ekin_tr()) ); }
      double Ps,s_var;
      s_var = 1.0;  if(is_thermostat){ s_var = thermostat->s_var; }
      Ps = 0.0;
      if(md->integrator=="Jacobi")    { top.propagate_exact_rb(dt_over_s); }
      else if(md->integrator=="DLML")  { top.propagate_dlml(dt_over_s,Ps); }
      else if(md->integrator=="Terec") { top.propagate_terec(dt_over_s);}
      else if(md->integrator=="qTerec") { top.propagate_qterec(dt_over_s);}
      else if(md->integrator=="NO_SQUISH"){ top.propagate_no_squish(dt_over_s);}
      else if(md->integrator=="KLN")   { top.propagate_kln(dt_over_s);}
      else if(md->integrator=="Omelyan"){ top.propagate_omelyan(dt_over_s);}

      if(is_thermostat){  thermostat->propagate_Ps( 0.5*dt_over_s2*(top.ekin_rot()+top.ekin_tr()) ); } 
      if(is_barostat) {  
        top.scale_position(sc1); 
        top.shift_position(sc2*top.rb_p*top.rb_iM);
      }
      else{
        top.shift_position(dt_over_s*top.rb_p*top.rb_iM);
      }        

    }// for i - all fragments

    if(is_thermostat){  thermostat->propagate_Ps(dt*( H0 - Nf*boltzmann*thermostat->Temperature*(log(thermostat->s_var)+1.0) ) ); }


    // Update cell shape
    if(is_barostat){ 
      if(syst->is_Box) {
        VECTOR t1_n,t2_n,t3_n; // new
        VECTOR t1_o,t2_o,t3_o; // old

        syst->Box  =  barostat->pos_scale(dt) * syst->Box;
/*
        system->Boxold.get_vectors(t1_o,t2_o,t3_o);
        system->Box.get_vectors(t1_n,t2_n,t3_n);

        // Find minimal heights of old and new boxes
        double h1_n,h2_n,h3_n;
        double h1_o,h2_o,h3_o;
        double V_n, V_o;  // Volume
        double max_n,min_n,max_o,min_o; // maximal and minimal heights new and old

        V_n = fabs(system->Box.Determinant());
        V_o = fabs(system->Boxold.Determinant());
        VECTOR S;
        S.cross(t2_n,t3_n); h1_n = V_n/fabs(S.length());
        S.cross(t3_n,t1_n); h2_n = V_n/fabs(S.length());
        S.cross(t1_n,t2_n); h3_n = V_n/fabs(S.length());
        max_n = (h1_n>=h2_n)?h1_n:h2_n;      min_n = (h1_n<=h2_n)?h1_n:h2_n;
        max_n = (max_n>=h3_n)?max_n:h3_n;    min_n = (min_n<=h3_n)?min_n:h3_n;

        S.cross(t2_o,t3_o); h1_o = V_o/fabs(S.length());
        S.cross(t3_o,t1_o); h2_o = V_o/fabs(S.length());
        S.cross(t1_o,t2_o); h3_o = V_o/fabs(S.length());
        max_o = (h1_o>=h2_o)?h1_o:h2_o;      min_o = (h1_o<=h2_o)?h1_o:h2_o;
        max_o = (max_o>=h3_o)?max_o:h3_o;    min_o = (min_o<=h3_o)?min_o:h3_o;

        int n_max = ceil((12.0+2.0)/min_n)+1;
        int o_max = ceil((12.0+2.0)/min_o)+1;
        n_max = (n_max>=o_max)?n_max:o_max;

        double dt1,dt2;
        dt1 = (max_n - min_o); dt1 = dt1*dt1;
        dt2 = (max_o - min_n); dt2 = dt2*dt2;
        system->dT_2 = 4.0*3.0*n_max*n_max*((dt1>=dt2)?dt1:dt2);
*/
      }
    }

    // Update atomic positions and calculate interactions
    for(i=0;i<syst->Number_of_fragments;i++){ syst->update_atoms_for_fragment(i);  }
    syst->zero_forces_and_torques();

//!!!!!!!!!!!!!!!!1    E_pot = system->energy(); !!!!!!!!!!!!!!!!1
    syst->extract_atomic_q(mol.q); // syst -> mol
//    E_pot = compute_potential_energy(mol, el, ham, 1); //  # 1 - FSSH forces
    E_pot = compute_forces(mol, el, ham, 1); 
    syst->set_atomic_f(mol.f);    // mol -> syst

    // Update rigid-body variables
    syst->update_fragment_forces_and_torques();


    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    
    //aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
    // Operator A(dt/2)
    E_kin = 0.0;
    //-------------------- Linear momentum propagation --------------------
    S = 0.0; I.identity();
    if(is_barostat){
      if(Nf_b==9){ S = barostat->ksi_eps + (barostat->ksi_eps.tr()/(barostat->get_Nf_t()/*+barostat->get_Nf_r()*/))*I; }
      else if(Nf_b==1){S = barostat->ksi_eps_iso * I + (3.0*barostat->ksi_eps_iso/(barostat->get_Nf_t()/*+barostat->get_Nf_r()*/))*I; }
    }
    if(is_thermostat){   S = S + thermostat->get_ksi_t() * I;      }
    sc1 = (exp_(S,-dt_half));//.symmetrized();
    sc2 = dt_half*(exp1_(S,-dt_half*0.5));//.symmetrized()*dt_half;

    //------------------- Angular momentum propagation -----------------------
    if(is_thermostat){ ksi_r = thermostat->get_ksi_r();}else{ ksi_r = 0.0;}
    sc3 = exp(-dt_half*ksi_r);
    sc4 = dt_half*exp(-0.5*dt_half*ksi_r)*sinh_(0.5*dt_half*ksi_r);


    for(i=0;i<syst->Number_of_fragments;i++){
      RigidBody& top = syst->Fragments[i].Group_RB;
      //-------------------- Linear momentum propagation --------------------
      top.scale_linear_(sc1);
      top.apply_force(sc2);
      //------------------- Angular momentum propagation -----------------------
      top.scale_angular_(sc3);
      top.apply_torque(sc4);
      E_kin += (top.ekin_rot() + top.ekin_tr());
    }// for i

    //aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

    if(is_thermostat){ thermostat->propagate_Ps( -dt_half*E_pot); }
   //------- Update state variables ------------
      curr_P_tens = syst->pressure_tensor();
      curr_P = (curr_P_tens.tr()/3.0);
      curr_V = syst->volume();
//   curr_P_tens = 0.0;
//   curr_P = 0.0;
//   curr_V = 1e+10;
   //-------------------------------------------

    //bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
    // Operator B(dt/2)
    if(is_barostat){
      if(md->ensemble=="NPT"||md->ensemble=="NPH"){ barostat->update_barostat_forces(syst->ekin_tr(),syst->ekin_rot(),curr_V,curr_P);   }
      else if(md->ensemble=="NPT_FLEX"||md->ensemble=="NPH_FLEX"){ barostat->update_barostat_forces(syst->ekin_tr(),syst->ekin_rot(),curr_V,curr_P_tens);   }
      scl = 0.0; if(is_thermostat){ scl = thermostat->get_ksi_b();  }
      barostat->propagate_velocity(dt_half,scl);
    }

    //bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb

    if(is_thermostat){thermostat->propagate_sPs(dt_half); }

    // Operator NHCB(dt/2)
    if(is_thermostat){
      double ekin_baro = 0.0;
      if(is_barostat){  ekin_baro = barostat->ekin_baro(); }
      thermostat->update_thermostat_forces(syst->ekin_tr(),syst->ekin_rot(),ekin_baro);
      thermostat->propagate_nhc(dt_half,syst->ekin_tr(),syst->ekin_rot(),ekin_baro);
    }


    if(md->ensemble=="NVE"){ s_var = 1.0;  Ps = 0.0;}
    if(is_thermostat){   E_kin/=(thermostat->s_var*thermostat->s_var); }

    E_kin_tr = syst->ekin_tr();
    E_kin_rot = syst->ekin_rot();
    E_tot = E_kin + E_pot;
//    if(!is_H0){ H0 = E_tot + thermostat->energy(); is_H0 = 1;}

    if(md->ensemble=="NVE"){  H_NP = E_tot; }
    else if(md->ensemble=="NVT"){
      if(is_thermostat){
        if(!is_H0){ H0 = E_tot + thermostat->energy(); is_H0 = 1;}
        if(thermostat->thermostat_type=="Nose-Poincare"){    H_NP = thermostat->s_var*(E_tot + thermostat->energy() - H0);   }
        else if(thermostat->thermostat_type=="Nose-Hoover"){ H_NP = E_tot + thermostat->energy();    }
      }
    }
    else if(md->ensemble=="NPH"||md->ensemble=="NPH_FLEX"){
      if(is_barostat){  H_NP = E_tot + barostat->ekin_baro() + curr_V * barostat->Pressure;   }
    }
    else if(md->ensemble=="NPT" || md->ensemble=="NPT_FLEX"){
      if(is_barostat){   H_NP = E_tot + barostat->ekin_baro() + curr_V * barostat->Pressure;  }
      if(is_thermostat){ H_NP += thermostat->energy();   }
    }

    curr_T = 2.0*E_kin/(Nf*(boltzmann/hartree));

    //------------- Angular velocity --------------
    L_tot = 0.0; 
    P_tot = 0.0;
    for(i=0;i<syst->Number_of_fragments;i++){
      RigidBody& top = syst->Fragments[i].Group_RB;
      VECTOR tmp; tmp.cross(top.rb_cm,top.rb_p);
      L_tot += top.rb_A_I_to_e_T * top.rb_l_e + tmp;
      P_tot += top.rb_p;
    }

    md->curr_step++;
    md->curr_time+=dt;
      
  }// for s
  md->curr_step = 0;
  md->curr_time = 0.0;
}


}// namespace libstate
}// namespace libscripts
}// liblibra

