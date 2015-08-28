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


double Barostat::ekin_baro(){
  double res = 0.0;
  if(Nf_b==9){ res = 0.5 * Wg * (ksi_eps.T() * ksi_eps).tr(); }
  else if(Nf_b==1){ res = 0.5 * Wg * ksi_eps_iso * ksi_eps_iso;}
  return res;
}

void Barostat::apply_barostat_force(double t){
  if(Nf_b==9){  ksi_eps = ksi_eps + t*G_eps;}
  else if(Nf_b==1){ ksi_eps_iso += t*G_eps_iso;}
}

void Barostat::scale_velocity(double sc){
  if(Nf_b==9) {ksi_eps = ksi_eps * sc; }
  else if(Nf_b==1){ksi_eps_iso *= sc; }
}

void Barostat::propagate_velocity(double dt,double ksi_b){
  if(Nf_b==1){      ksi_eps_iso = exp(-dt*ksi_b)*ksi_eps_iso + dt*exp(-0.5*dt*ksi_b)*sinh_(0.5*dt*ksi_b)*G_eps_iso; }
  else if(Nf_b==9){ ksi_eps     = exp(-dt*ksi_b)*ksi_eps     + dt*exp(-0.5*dt*ksi_b)*sinh_(0.5*dt*ksi_b)*G_eps; }
}


void Barostat::update_barostat_forces(double ekin_tr,double ekin_rot,double curr_volume, double curr_pressure){
  G_eps_iso = (3.0/Wg)*(curr_volume * (curr_pressure - Pressure) + (ekin_tr/Nf_t) /*((ekin_tr + ekin_rot)/(Nf_t + Nf_r))*/ );
}

void Barostat::update_barostat_forces(double ekin_tr,double ekin_rot,double curr_volume, MATRIX3x3& curr_pressure_tensor){
  // Symmetrize pressure tensor
  MATRIX3x3 current_pressure_tensor = 0.5*(curr_pressure_tensor + curr_pressure_tensor.T());
  MATRIX3x3 I; I.identity();
  double curr_pressure = (1.0/3.0)*current_pressure_tensor.tr();

  G_eps = (curr_volume * (curr_pressure_tensor - Pressure * I) + /*((ekin_tr + ekin_rot)/(Nf_t + Nf_r))*/ (ekin_tr/Nf_t)* I)/Wg;
  G_eps_iso = (3.0/Wg)*(curr_volume * (curr_pressure - Pressure) + (ekin_tr/Nf_t)/*((ekin_tr + ekin_rot)/(Nf_t + Nf_r))*/ );

}

void Barostat::init(double Temperature){
  double kTb = ((boltzmann/hartree) * Temperature / (nu_baro * nu_baro));
  Wg = (Nf_t + Nf_r + 3.0)*kTb/3.0;
}

MATRIX3x3 Barostat::pos_scale(double dt){
  MATRIX3x3 res; res.identity();
  if(Nf_b==9) {  res = (exp_(ksi_eps,dt));  }
  else if(Nf_b==1){ res = exp(ksi_eps_iso*dt)*res; }
  return res;
}

MATRIX3x3 Barostat::vpos_scale(double dt){
  MATRIX3x3 res; res.identity();
  if(Nf_b==9){ res = (exp1_(ksi_eps,dt));  }
  else if(Nf_b==1){   res = res * exp(0.5*dt*ksi_eps_iso)*sinh_(0.5*dt*ksi_eps_iso); /*exp1_(ksi_eps,dt)*/; }
  return res;
}

MATRIX3x3 Barostat::vel_scale(double dt,double ksi_t){
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


