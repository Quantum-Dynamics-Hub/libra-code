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
  \file Thermostat_methods.cpp
  \brief The file implements the methods to propagate Thermostat state
    
*/

#include "Thermostat.h"

/// liblibra namespace
namespace liblibra{


/// libdyn namespace
namespace libdyn{

/// libthermostat namespace
namespace libthermostat{


double Thermostat::energy(){
/** 
  \brief Return the energy of Thermostat

  So far only Nose-Hoover (including chain) and Nose-Poincare thermostats are possible
*/

  double comp = 0.0;
  int i;
  if(thermostat_type=="Nose-Hoover"){

  double kT = (boltzmann/hartree)*Temperature;
  if(Nf_t>0){
  for(i=0;i<NHC_size;i++){
      if(i==0){ comp += kT * Nf_t * s_t[i]; }
      else{     comp += kT * s_t[i]; }
      comp += ( 0.5*Q_t[i] * ksi_t[i] * ksi_t[i] );
  }  }

  if(Nf_r>0){
  for(i=0;i<NHC_size;i++){
      if(i==0){ comp += kT * Nf_r * s_r[i]; }
      else{     comp += kT * s_r[i]; }
      comp += ( 0.5*Q_r[i] * ksi_r[i] * ksi_r[i] );
  }  }

  if(Nf_b>0){
  for(i=0;i<NHC_size;i++){
      if(i==0){ comp += kT * Nf_b * s_b[i]; }
      else{     comp += kT * s_b[i]; }
      comp += ( 0.5*Q_b[i] * ksi_b[i] * ksi_b[i] );
  }  }

  }// if Nose-Hoover

  else if(thermostat_type=="Nose-Poincare"){
    comp = (0.5*Ps*Ps/Q)+ (Nf_t + Nf_r)*(boltzmann/hartree)*Temperature*log(s_var); 
  }

  return comp;
}


void Thermostat::propagate_sPs(double t){
/**
  \brief Propagate s and Ps variables - only for the case of Nose-Poincare thermostat
  \param[in] t The propagation time

  Action of the operator exp(t*D(H_3)), where 
  D(H_3) = (s*Ps/Q)(d/ds) - (Ps^2/2Q)(d/dPs), and
  H_3 = s*Ps^2/2Q
  on the state variables s_var and Ps
  according to:   Nose, S. "An Improved Symplectic
  Integrator for Nose-Poincare Thermostat" 
  2001, JPSJ, 70, 75-77

*/
  if(thermostat_type=="Nose-Poincare"){
    double tmp = (1.0+(0.5*t*Ps/Q));
    s_var *= (tmp*tmp);
    Ps = Ps/tmp;
  }
}

void Thermostat::propagate_Ps(double amnt){
/**
  \brief Propagate Ps variable - only for Nose-Poincare thermostat
  \param[in] amnt Amount of shift: Ps -> Ps + amnt
*/

  if(thermostat_type=="Nose-Poincare"){
    Ps += amnt;
  }
}

double Thermostat::vel_scale(double dt){
/**
  \brief Return velocity (momentum) scaling coeffieicent.
  \param[in] dt The propagation time
  
  Used to constructe NVT-MD algorithm. Only for Nose-Hoover thermostat
*/
  double res = 1.0;
  if(thermostat_type=="Nose-Hoover"){  res = exp(-dt*ksi_t[0]);  }
  return res;
}

double Thermostat::ang_vel_scale(double dt){
/**
  \brief Return angular velocity (angular momentum) scaling coeffieicent.
  \param[in] dt The propagation time
  
  Used to constructe NVT-MD algorithm. Only for Nose-Hoover thermostat
*/

  double res = 1.0;
  if(thermostat_type=="Nose-Hoover"){  res = exp(-dt*ksi_r[0]);  }
  return res;
}


void Thermostat::update_thermostat_forces(double ekin_tr, double ekin_rot,double ekin_baro){
/**
  \brief This function calculates the thermostat forces for Nose-Hoover chain thermostat
  \param[in] ekin_tr Translational kinetic energy of the system to which thermostat is coupled
  \param[in] ekin_rot Rotational kinetic energy of the system to which thermostat is coupled
  \param[in] ekin_baro Barostat kinetic energy of the system to which thermostat is coupled

  This procedure is applied only to Nose-Hoover thermostat

  As described in: Kamberaj, H.; Low, R. J; Neal, M. P. "Time reversible and symplectic
  integrators for molecular dynamics simulations of rigid molecules" J. Chem. Phys. 2005, 122, 224114-1 - 224114-30
  Note: The indexing for barostat chains is slightly different than in the original article
*/
  if(thermostat_type=="Nose-Hoover"){
  double kT = (boltzmann/hartree) * Temperature;

  for(int i=0;i<NHC_size;i++){
    if(i==0){
      if(Nf_t>0){G_t[0] = (2.0*ekin_tr  - Nf_t * kT)/Q_t[i];}//else{ G_t[0] = 0.0; }
      if(Nf_r>0){G_r[0] = (2.0*ekin_rot - Nf_r * kT)/Q_r[i];}//else{ G_r[0] = 0.0; }
      if(Nf_b>0){G_b[0] = (2.0*ekin_baro- Nf_b * kT)/Q_b[0];}//else{ G_b[0] = 0.0; }
    }else{
      if(Nf_t>0){G_t[i] = (Q_t[i-1] * ksi_t[i-1] * ksi_t[i-1] - kT)/Q_t[i];}//else{ G_t[i] = 0.0; } 
      if(Nf_r>0){G_r[i] = (Q_r[i-1] * ksi_r[i-1] * ksi_r[i-1] - kT)/Q_r[i];}//else{ G_r[i] = 0.0; }
      if(Nf_b>0){G_b[i] = (Q_b[i-1] * ksi_b[i-1] * ksi_b[i-1] - kT)/Q_b[i];}//else{ G_b[i] = 0.0; }
    }
  }// for i
  }// if Nose-Hoover
}

void Thermostat::update_thermostat_forces(double ekin_tr, double ekin_rot,double ekin_baro,int i){
/**
  \brief This function calculates the thermostat forces for Nose-Hoover chain thermostat
  \param[in] ekin_tr Translational kinetic energy of the system to which thermostat is coupled
  \param[in] ekin_rot Rotational kinetic energy of the system to which thermostat is coupled
  \param[in] ekin_baro Barostat kinetic energy of the system to which thermostat is coupled
  \param[in] i The index of the thermostat in the chain for which (thermostat) we want to compute the force

  This procedure is applied only to Nose-Hoover thermostat

  As described in: Kamberaj, H.; Low, R. J; Neal, M. P. "Time reversible and symplectic
  integrators for molecular dynamics simulations of rigid molecules" J. Chem. Phys. 2005, 122, 224114-1 - 224114-30
  Note: The indexing for barostat chains is slightly different than in the original article

  This function calculates the thermostat forces only for the i-th thermostat in the chain

*/

  if(thermostat_type=="Nose-Hoover"){
  double kT = (boltzmann/hartree) * Temperature;
    if(i==0){
      if(Nf_t>0){G_t[i] = (2.0*ekin_tr  - Nf_t * kT)/Q_t[i];}//else{ G_t[0] = 0.0; }
      if(Nf_r>0){G_r[i] = (2.0*ekin_rot - Nf_r * kT)/Q_r[i];}//else{ G_r[0] = 0.0; }
      if(Nf_b>0){G_b[i] = (2.0*ekin_baro- Nf_b * kT)/Q_b[i];}//else{ G_b[0] = 0.0; }
    }else{
      if(Nf_t>0){G_t[i] = (Q_t[i-1] * ksi_t[i-1] * ksi_t[i-1] - kT)/Q_t[i];}//else{ G_t[i] = 0.0; } 
      if(Nf_r>0){G_r[i] = (Q_r[i-1] * ksi_r[i-1] * ksi_r[i-1] - kT)/Q_r[i];}//else{ G_r[i] = 0.0; }
      if(Nf_b>0){G_b[i] = (Q_b[i-1] * ksi_b[i-1] * ksi_b[i-1] - kT)/Q_b[i];}//else{ G_b[i] = 0.0; }
    }
  }// if Nose-Hoover
}


void Thermostat::init_nhc(){
/**
  \brief This function calculates the masses for Nose-Hoover chain thermostat and allocates the memory for corresponding dynamic variables 
  As described in: Kamberaj, H.; Low, R. J; Neal, M. P. "Time reversible and symplectic
  integrators for molecular dynamics simulations of rigid molecules" J. Chem. Phys. 2005, 122, 224114-1 - 224114-30

  This function initializes variables, masses, and computes thermostat forces
  Note: For distinguishing between different barostats we use the number of degrees of freedom corresponding to barostat.
  
*/

  int i;
  double Qt,Qr,Qb;

  if(thermostat_type=="Nose-Hoover"){
  double kTt = ((boltzmann/hartree) * Temperature / (nu_therm * nu_therm));

  // Clear all variables
  if(s_t_size>0)  { s_t.clear();   s_t_size = 0; }
  if(s_r_size>0)  { s_r.clear();   s_r_size = 0; }
  if(s_b_size>0)  { s_b.clear();   s_b_size = 0; }
  if(G_t_size>0)  { G_t.clear();   G_t_size = 0; }
  if(G_r_size>0)  { G_r.clear();   G_r_size = 0; }
  if(G_b_size>0)  { G_b.clear();   G_b_size = 0; }
  if(Q_t_size>0)  { Q_t.clear();   Q_t_size = 0; }
  if(Q_r_size>0)  { Q_r.clear();   Q_r_size = 0; }
  if(Q_b_size>0)  { Q_b.clear();   Q_b_size = 0; }
  if(ksi_t_size>0){ ksi_t.clear(); ksi_t_size = 0;}
  if(ksi_r_size>0){ ksi_r.clear(); ksi_r_size = 0;}
  if(ksi_b_size>0){ ksi_b.clear(); ksi_b_size = 0;}

  // Now calculate masses, initialize s and ksi
  // and allocate memory for G
  if(Nf_t>0){
  for(i=0;i<NHC_size;i++){
    if(i==0){   Qt = Nf_t * kTt;}
    else{ Qt = kTt; }
    Q_t.push_back(Qt);
    s_t.push_back(0.0);
    ksi_t.push_back(0.0);
    G_t.push_back(0.0);
  }// for i
    Q_t_size = NHC_size;
    s_t_size = NHC_size;
    ksi_t_size = NHC_size;
    G_t_size = NHC_size;
  }// if Nf_t>o

  if(Nf_r>0){
  for(i=0;i<NHC_size;i++){
    if(i==0){   Qr = Nf_r * kTt;}
    else{ Qr = kTt; }
    Q_r.push_back(Qr);
    s_r.push_back(0.0);
    ksi_r.push_back(0.0);
    G_r.push_back(0.0);
  }// for i
    Q_r_size = NHC_size;
    s_r_size = NHC_size;
    ksi_r_size = NHC_size;
    G_r_size = NHC_size;
  }// if Nf_r>o

  if(Nf_b>0){
  for(i=0;i<NHC_size;i++){
    if(i==0){   Qb = Nf_b * kTt;}
    else{ Qb = kTt; }
    Q_b.push_back(Qb);
    s_b.push_back(0.0);
    ksi_b.push_back(0.0);
    G_b.push_back(0.0);
  }// for i
    Q_b_size = NHC_size;
    s_b_size = NHC_size;
    ksi_b_size = NHC_size;
    G_b_size = NHC_size;
  }// if Nf_b>o

  update_thermostat_forces(0.0,0.0,0.0);

  }// if Nose-Hoover


}

void Thermostat::propagate_nhc(double dt,double ekin_tr, double ekin_rot,double ekin_baro){
/**
  \brief This function provides an approximate solution for the Nose-Hoover chain thermostat via symplectic splitting 
  \param[in] dt The integration time
  \param[in] ekin_tr Translational kinetic energy of the system to which thermostat is coupled
  \param[in] ekin_rot Rotational kinetic energy of the system to which thermostat is coupled
  \param[in] ekin_baro Barostat kinetic energy of the system to which thermostat is coupled

  As described in: Kamberaj, H.; Low, R. J; Neal, M. P. "Time reversible and symplectic
  integrators for molecular dynamics simulations of rigid molecules" J. Chem. Phys. 2005, 122, 224114-1 - 224114-30

  Note: The indexing for barostat chains is slightly different than in original article

*******************************************************************/


  double argt,argr,argb;
  double et,er,eb;

  if(thermostat_type=="Nose-Hoover"){
  int M = NHC_size - 1;
  int k;

  update_thermostat_forces(ekin_tr,ekin_rot,ekin_baro,M);
  if(Nf_t>0){ ksi_t[M] = ksi_t[M] + 0.5*dt*G_t[M]; }
  if(Nf_r>0){ ksi_r[M] = ksi_r[M] + 0.5*dt*G_r[M]; }
  if(Nf_b>0){ ksi_b[M] = ksi_b[M] + 0.5*dt*G_b[M]; }

  for(k=1;k<=M;k++){
    update_thermostat_forces(ekin_tr,ekin_rot,ekin_baro,M-k);
    if(Nf_t>0){
      argt = 0.25*dt*ksi_t[M-k+1];  et = exp(-argt);  
      ksi_t[M-k] = et*(et*ksi_t[M-k] + 0.5*dt*G_t[M-k]*sinh_(argt));
    }
    if(Nf_r>0){
      argr = 0.25*dt*ksi_r[M-k+1];  er = exp(-argr);
      ksi_r[M-k] = er*(er*ksi_r[M-k] + 0.5*dt*G_r[M-k]*sinh_(argr));
    }
    if(Nf_b>0){
      argb = 0.25*dt*ksi_b[M-k+1];  eb = exp(-argb);
      ksi_b[M-k] = eb*(eb*ksi_b[M-k] + 0.5*dt*G_b[M-k]*sinh_(argb));
    }
  }// for k

  for(k=0;k<=M;k++){  
    if(Nf_t>0){ s_t[k] = s_t[k] + dt*ksi_t[k]; }
    if(Nf_r>0){ s_r[k] = s_r[k] + dt*ksi_r[k]; }
    if(Nf_b>0){ s_b[k] = s_b[k] + dt*ksi_b[k]; }
  }
  
  for(k=0;k<=(M-1);k++){
    update_thermostat_forces(ekin_tr,ekin_rot,ekin_baro,k);
    if(Nf_t>0){
      argt = 0.25*dt*ksi_t[k+1];  et = exp(-argt);
      ksi_t[k] = et*(et*ksi_t[k] + 0.5*dt*G_t[k]*sinh_(argt));
    }
    if(Nf_r>0){
      argr = 0.25*dt*ksi_r[k+1];  er = exp(-argr);
      ksi_r[k] = er*(er*ksi_r[k] + 0.5*dt*G_r[k]*sinh_(argr));
    }
    if(Nf_b>0){
      argb = 0.25*dt*ksi_b[k+1];  eb = exp(-argb);
      ksi_b[k] = eb*(eb*ksi_b[k] + 0.5*dt*G_b[k]*sinh_(argb));
    }
  }// for k

  update_thermostat_forces(ekin_tr,ekin_rot,ekin_baro,M);
  if(Nf_t>0){ ksi_t[M] = ksi_t[M] + 0.5*dt*G_t[M]; }
  if(Nf_r>0){ ksi_r[M] = ksi_r[M] + 0.5*dt*G_r[M]; }
  if(Nf_b>0){ ksi_b[M] = ksi_b[M] + 0.5*dt*G_b[M]; }

  }//if Nose-Hoover

}

void Thermostat::cool(){
/**
  \brief Cool down the thermostat

  This is done by setting Ps = 0 and s_var = 1 for Nose-Poincare thermostat,
  and by re-initializing NHC variables for the Nose-Hoover chain thermostats
*/

  if(thermostat_type=="Nose-Poincare"){
    Ps = 0.0;     is_Ps = 1;
    s_var = 1.0;  is_s_var = 1;   
  }
  else if(thermostat_type=="Nose-Hoover"){
    init_nhc();
  }

}



}// namespace libthermostat
}// namespace libdyn

}// liblibra

