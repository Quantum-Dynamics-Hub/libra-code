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
  \file Barostat_methods.cpp
  \brief The file implements basic methods of the Barostat class
    
*/

#include "Barostat.h"

/// liblibra namespace
namespace liblibra{


/// libdyn namespace
namespace libdyn{

/// libbarostat namespace
namespace libbarostat{


double Barostat::ekin_baro() const{
/**
  \brief Kinetic energy of the Barostat veriables

  Defferent schemes are used for isotropic (Nf_b = 1 ) and fully-flexible (Nf_b = 9) cells
*/

  double res = 0.0;
  if(Nf_b==9){ res = 0.5 * Wg * (ksi_eps.T() * ksi_eps).tr(); }
  else if(Nf_b==1){ res = 0.5 * Wg * ksi_eps_iso * ksi_eps_iso;}
  return res;
}


void Barostat::apply_barostat_force(double t){
/** 
  \brief Apply Barostat forces to change the Barostat velocities

  Different schemes are applied for isotropic (Nf_b = 1 ) and fully-flexible (Nf_b = 9) cells

  \param[in] t The force scaling parameter (usually related to the integration time step)
*/

  if(Nf_b==9){  ksi_eps = ksi_eps + t*G_eps;}
  else if(Nf_b==1){ ksi_eps_iso += t*G_eps_iso;}
}


void Barostat::scale_velocity(double sc){
/** 
  \brief Scale the barostat velocity by a factor sc
  
  Different schemes are applied for isotropic (Nf_b = 1 ) and fully-flexible (Nf_b = 9) cells

  \param[in] sc The scaling parameter
*/

  if(Nf_b==9) {ksi_eps = ksi_eps * sc; }
  else if(Nf_b==1){ksi_eps_iso *= sc; }
}

void Barostat::propagate_velocity(double dt,double ksi_b){
/**
  \brief Propagate barostat velocities 

  Different schemes are applied for isotropic (Nf_b = 1 ) and fully-flexible (Nf_b = 9) cells
  Based on: Kamberaj, H.; Low, R. J.; Neal, M. P. Time Reversible and Symplectic Integrators for Molecular Dynamics Simulations of Rigid Molecules. J. Chem. Phys. 2005, 122, 224114.
 
  \param[in] dt The integration time step
  \param[in] ksi_b Is the thermostat velocity (if coupled to Thermostat as well), 0 is no Thermostat

*/
  if(Nf_b==1){      ksi_eps_iso = exp(-dt*ksi_b)*ksi_eps_iso + dt*exp(-0.5*dt*ksi_b)*sinh_(0.5*dt*ksi_b)*G_eps_iso; }
  else if(Nf_b==9){ ksi_eps     = exp(-dt*ksi_b)*ksi_eps     + dt*exp(-0.5*dt*ksi_b)*sinh_(0.5*dt*ksi_b)*G_eps; }
}


void Barostat::update_barostat_forces(double ekin_tr,double ekin_rot,double curr_volume, double curr_pressure){
/**
  \brief Compute isotropic Barostat forces

  Based on: Kamberaj, H.; Low, R. J.; Neal, M. P. Time Reversible and Symplectic Integrators for Molecular Dynamics Simulations of Rigid Molecules. J. Chem. Phys. 2005, 122, 224114.
 
  \param[in] ekin_tr Translational kinetic energy
  \param[in] ekin_rot Rotational kinetic energy
  \param[in] curr_volume Current volume of the simulation cell
  \param[in] curr_pressure Current pressure of the system

*/

  G_eps_iso = (3.0/Wg)*(curr_volume * (curr_pressure - Pressure) + (ekin_tr/Nf_t) /*((ekin_tr + ekin_rot)/(Nf_t + Nf_r))*/ );
}

void Barostat::update_barostat_forces(double ekin_tr,double ekin_rot,double curr_volume, MATRIX3x3& curr_pressure_tensor){
/**
  \brief Compute isotropic Barostat forces, but from the anisotropic simulation cell variables

  Based on: Kamberaj, H.; Low, R. J.; Neal, M. P. Time Reversible and Symplectic Integrators for Molecular Dynamics Simulations of Rigid Molecules. J. Chem. Phys. 2005, 122, 224114.
  For stability, I symmetrize the pressure tensor before further calculations
 
  \param[in] ekin_tr Translational kinetic energy
  \param[in] ekin_rot Rotational kinetic energy
  \param[in] curr_volume Current volume of the simulation cell
  \param[in] curr_pressure_tensor Current pressure tensor of the system

*/

  // Symmetrize pressure tensor
  MATRIX3x3 current_pressure_tensor = 0.5*(curr_pressure_tensor + curr_pressure_tensor.T());
  MATRIX3x3 I; I.identity();
  double curr_pressure = (1.0/3.0)*current_pressure_tensor.tr();

  G_eps = (curr_volume * (curr_pressure_tensor - Pressure * I) + /*((ekin_tr + ekin_rot)/(Nf_t + Nf_r))*/ (ekin_tr/Nf_t)* I)/Wg;
  G_eps_iso = (3.0/Wg)*(curr_volume * (curr_pressure - Pressure) + (ekin_tr/Nf_t)/*((ekin_tr + ekin_rot)/(Nf_t + Nf_r))*/ );

}

void Barostat::init(double Temperature){
/**
  \brief Initialize the Barostat mass

  Based on: Kamberaj, H.; Low, R. J.; Neal, M. P. Time Reversible and Symplectic Integrators for Molecular Dynamics Simulations of Rigid Molecules. J. Chem. Phys. 2005, 122, 224114.
 
  \param[in] Temperature Is the target temperature (for thermostat)

*/

  double kTb = ((boltzmann/hartree) * Temperature / (nu_baro * nu_baro));
  Wg = (Nf_t + Nf_r + 3.0)*kTb/3.0;
}

MATRIX3x3 Barostat::pos_scale(double dt){
/**
  \brief Returns an operator to propagate Barostat position variable

  Different schemes are applied for isotropic (Nf_b = 1 ) and fully-flexible (Nf_b = 9) cells
  Based on: Kamberaj, H.; Low, R. J.; Neal, M. P. Time Reversible and Symplectic Integrators for Molecular Dynamics Simulations of Rigid Molecules. J. Chem. Phys. 2005, 122, 224114.
 
  \param[in] dt Is integration time step

*/

  MATRIX3x3 res; res.identity();
  if(Nf_b==9) {  res = (exp_(ksi_eps,dt));  }
  else if(Nf_b==1){ res = exp(ksi_eps_iso*dt)*res; }
  return res;
}

MATRIX3x3 Barostat::vpos_scale(double dt){
/**
  \brief Returns and operator to propagate Barostat velocity variable

  Different schemes are applied for isotropic (Nf_b = 1 ) and fully-flexible (Nf_b = 9) cells
  Based on: Kamberaj, H.; Low, R. J.; Neal, M. P. Time Reversible and Symplectic Integrators for Molecular Dynamics Simulations of Rigid Molecules. J. Chem. Phys. 2005, 122, 224114.
 
  \param[in] dt Is integration time step

*/

  MATRIX3x3 res; res.identity();
  if(Nf_b==9){ res = (exp1_(ksi_eps,dt));  }
  else if(Nf_b==1){   res = res * exp(0.5*dt*ksi_eps_iso)*sinh_(0.5*dt*ksi_eps_iso); /*exp1_(ksi_eps,dt)*/; }
  return res;
}

MATRIX3x3 Barostat::vel_scale(double dt,double ksi_t){
/**
  \brief Returns operator to scale Barostat velocity variable (based on translational velocity of Thermostat)

  Different schemes are applied for isotropic (Nf_b = 1 ) and fully-flexible (Nf_b = 9) cells
  Based on: Kamberaj, H.; Low, R. J.; Neal, M. P. Time Reversible and Symplectic Integrators for Molecular Dynamics Simulations of Rigid Molecules. J. Chem. Phys. 2005, 122, 224114.
 
  \param[in] dt Is integration time step
  \param[in] ksi_t Thermostat translational velocity (use 0 if there is no translational Thermostat applied)  

*/

  double Nf = Nf_t;// + Nf_r;
  MATRIX3x3 I,res; I.identity();
  if(Nf_b==9){ 
    MATRIX3x3 tmp = ( ksi_eps + ((ksi_eps.tr()/Nf) + ksi_t)*I);
    res = (exp_( tmp , -dt));
  }
  else if(Nf_b==1){
    res = exp(-dt*(ksi_t+ (1.0 + (3.0/Nf))*ksi_eps_iso)  )*I;
  }
  return res;
}

MATRIX3x3 Barostat::ang_vel_scale(double dt,double ksi_r){
/**
  \brief Returns operator to scale Barostat velocity variable (based on rotational velocity of Thermostat)

  Different schemes are applied for isotropic (Nf_b = 1 ) and fully-flexible (Nf_b = 9) cells
  Based on: Kamberaj, H.; Low, R. J.; Neal, M. P. Time Reversible and Symplectic Integrators for Molecular Dynamics Simulations of Rigid Molecules. J. Chem. Phys. 2005, 122, 224114.
 
  \param[in] dt Is integration time step
  \param[in] ksi_r Thermostat rotational velocity (use 0 if there is no rotational Thermostat applied)  

*/

  double Nf = Nf_t;// + Nf_r;
  MATRIX3x3 I,res; I.identity();
  if(Nf_b==9){ res =  exp(  -dt*(ksi_r/*+(ksi_eps.tr()/Nf)*/)  ) * I; }
  else if(Nf_b==1){ res = exp(-dt*(ksi_r/*+ (3.0/Nf)*ksi_eps_iso*/) ) * I; }
  return res;
}

void Barostat::cool(){
  ksi_eps = 0.0;     is_ksi_eps = 1;
  ksi_eps_iso = 0.0; is_ksi_eps_iso = 1;
}



}// namespace libbarostat
}// namespace libdyn
}// liblibra

