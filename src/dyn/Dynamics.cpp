/*********************************************************************************
* Copyright (C) 2019-2023 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dynamics.cpp
  \brief The file implements the general framework to run:
   * adiabatic dynamics (Verlet)
   * nonadiabatic Ehrenfest dynamics
   * nonadiabatic TSH dynamics 
   * include thermostat, if needed
   * include decoherence
   * include quantum nuclear effect (ETHD) 
   * include phase corrections/state tracking in NA-MD
   * the same framework for multiple trajectories
   * include coupled-trajectories methods (planned)
   * enable the NBRA-like calculations as well as non-NBRA
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"
#include "electronic/libelectronic.h"
#include "Dynamics.h"
#include "dyn_control_params.h"
#include "dyn_variables.h"
#include "dyn_ham.h"
#include "../calculators/NPI.h"


/// liblibra namespace
namespace liblibra{

using namespace libcalculators;

/// libdyn namespace
namespace libdyn{



namespace bp = boost::python;


void aux_get_transforms(CMATRIX** Uprev, nHamiltonian& ham){

  // For adiabatic representation only:
  // Save the previous orbitals info - in case we need to
  // do either phase correction of state tracking

  int ntraj = ham.children.size();

  for(int traj=0; traj<ntraj; traj++){
    *Uprev[traj] = ham.children[traj]->get_basis_transform();  
  }

}





vector<CMATRIX> compute_St(nHamiltonian& ham, nHamiltonian& ham_prev, int isNBRA){

  return compute_St(&ham, &ham_prev, isNBRA);
}

vector<CMATRIX> compute_St(nHamiltonian* ham, nHamiltonian* ham_prev, int isNBRA){
/**
  This function computes the time-overlap matrices for all trajectories

*/

  int nst = ham->nadi;
  int ntraj = ham->children.size();

  vector<CMATRIX> St(ntraj, CMATRIX(nst, nst));

  if(isNBRA==1){
    St[0] = ham_prev->children[0]->get_basis_transform().H() * ham->children[0]->get_basis_transform();    
    // AVA: temporarily comment on 11/25/2022
    ham->children[0]->set_time_overlap_adi_by_val(St[0]);
  }
  else{
    for(int traj=0; traj<ntraj; traj++){
      St[traj] = ham_prev->children[traj]->get_basis_transform().H() * ham->children[traj]->get_basis_transform();
      ham->children[traj]->set_time_overlap_adi_by_val(St[traj]);      
    }
  }
  return St;

}

vector<CMATRIX> compute_St(nHamiltonian& ham, nHamiltonian& ham_prev){

  int is_nbra = 0;
  return compute_St(&ham, &ham_prev, is_nbra);

}

vector<CMATRIX> compute_St(nHamiltonian* ham, int isNBRA){
/**
  This function computes the time-overlap matrices for all trajectories
*/

  int nst = ham->nadi;
  int ntraj = ham->children.size();

  vector<CMATRIX> St(ntraj, CMATRIX(nst, nst));
  if(isNBRA==1){
    St[0] = ham->children[0]->get_time_overlap_adi();
  }
  else{
    for(int traj=0; traj<ntraj; traj++){
      St[traj] = ham->children[traj]->get_time_overlap_adi();
    }
  }
  return St;

}


vector<CMATRIX> compute_St(nHamiltonian& ham, int isNBRA){

  return compute_St(&ham, isNBRA);
}

vector<CMATRIX> compute_St(nHamiltonian& ham){
  int is_nbra = 0;

  return compute_St(&ham, is_nbra);
}



void apply_afssh(dyn_variables& dyn_var, CMATRIX& C, vector<int>& act_states, MATRIX& invM,
                nHamiltonian& ham, bp::dict& dyn_params, Random& rnd  ){

  //cout<<"In apply_afssh\n";

  dyn_control_params prms;
  prms.set_parameters(dyn_params);

  int i;
  int ndof = invM.n_rows;
  int nst = C.n_rows;    
  int ntraj = C.n_cols;
  int traj, idof;
  int num_el = prms.num_electronic_substeps;
  double dt_el = prms.dt / num_el;

  // A-FSSH

    CMATRIX hvib_curr(nst, nst);
    CMATRIX force_full(nst, nst);
    CMATRIX force_diag(nst, nst);
    CMATRIX c_traj(nst, 1);
    CMATRIX dR_afssh(nst, nst);
    CMATRIX dP_afssh(nst, nst);



    //cout<<"Propagating moments...\n";
    //=========================== Propagate moments ===============
    for(traj=0; traj<ntraj; traj++){

      hvib_curr = ham.children[traj]->get_hvib_adi();
      c_traj = C.col(traj);

      double gamma_reset = 0.0;

      for(idof=0; idof<ndof; idof++){  

         force_full = -1.0 * ham.children[traj]->get_d1ham_adi(idof);

         for(i=0;i<nst;i++){  force_diag.set(i, i,  force_full.get(i,i) ); } 

         dR_afssh = *dyn_var.dR[traj][idof];
         dP_afssh = *dyn_var.dP[traj][idof];

         integrate_afssh_moments(dR_afssh, dP_afssh, hvib_curr, force_diag, 
                                 c_traj, 1.0/invM.get(idof,0), act_states[traj], dt_el, num_el);

         *dyn_var.dR[traj][idof] = dR_afssh; 
         *dyn_var.dP[traj][idof] = dP_afssh;         

      }// for idof
    }// for traj


    //cout<<"Computing reset and collapse probabilities...\n";

    //======================== Compute reset and collapse probabilities =========

    MATRIX gamma_reset(nst, ntraj);
    MATRIX gamma_collapse(nst, ntraj);


    for(traj=0; traj<ntraj; traj++){
      for(i=0;i<nst;i++){

        double gamma_reset_i = 0.0;
        double gamma_collapse_i = 0.0;

        for(idof=0; idof<ndof; idof++){  
        
          double dx_ii = dR_afssh.get(i, i).real();
          int as = act_states[traj];
          double f_i   = -ham.children[traj]->get_d1ham_adi(idof).get(i, i).real();
          double f_as = -ham.children[traj]->get_d1ham_adi(idof).get(as,as).real();

          gamma_reset_i -= 0.5*(f_i - f_as) * dx_ii;

          double f_ji   = -ham.children[traj]->get_d1ham_adi(idof).get(as, i).real();
          gamma_collapse_i += f_ji * dx_ii;

        }// for idof

        gamma_reset.set(i, traj, gamma_reset_i * prms.dt);
        gamma_collapse.set(i, traj, (gamma_reset_i -  2.0*fabs(gamma_collapse_i)) * prms.dt );

      }// for nst
    }// for traj


    //cout<<"Doing the collapses and resets...\n";
    //======================== Do the collapse and resets =======================

    complex<double> zero(0.0, 0.0);

    for(traj=0; traj<ntraj; traj++){

      for(i=0;i<nst;i++){

        if(i!=act_states[traj]){


          // Reset
          double ksi = rnd.uniform(0.0, 1.0);

          if(ksi < gamma_reset.get(i, traj)){
            for(idof=0;idof<ndof;idof++){
              dyn_var.dR[traj][idof]->scale(-1, i, zero);
              dyn_var.dR[traj][idof]->scale(i, -1, zero);    
              dyn_var.dP[traj][idof]->scale(-1, i, zero);
              dyn_var.dP[traj][idof]->scale(i, -1, zero);
            }// for j
          }// if ksi < gamma_reset


          // Collapse
          ksi = rnd.uniform(0.0, 1.0);
          if(ksi < gamma_collapse.get(i, traj)){
            collapse(C, traj, act_states[traj], prms.collapse_option);
          }// if ksi < gamma_collapse
          

        }// all non-active states

      }// for i
    }// for traj

    // cout<<"Done\n";
}


void apply_projectors(CMATRIX& C, vector<CMATRIX>& proj){

  int ntraj = proj.size();
  int nst = C.n_rows;

  int i;
  vector<int> t2(1,0);
  vector<int> t3(nst, 0); for(i=0;i<nst;i++){  t3[i] = i; }
  CMATRIX c_tmp(nst, 1);

  for(int traj=0; traj<ntraj; traj++){
    t2[0] = traj;
    pop_submatrix(C, c_tmp, t3, t2);
    c_tmp = proj[traj] * c_tmp;
    push_submatrix(C, c_tmp, t3, t2);
  }
}


CMATRIX Zhu_Liouvillian(double Etot, CMATRIX& Ham, CMATRIX& rho){
/**
  This is the Liouvillian based on Eqs. 3.34 and 3.35 of my Chapter, but withing the local diabatization
  approach, that is when NACs are dropped and Hvib is replaced by Hadi. 
*/

  int nst = Ham.n_cols;
  int sz = nst * nst;

  CMATRIX L(sz, sz);
  int ij;

  double Eeff = 0.0;
  for(int i=0; i<nst; i++){
    Eeff += Ham.get(i,i).real() * rho.get(i,i).real(); 
  }
  Eeff = Etot - Eeff;
  if(Eeff < 0.0){ Eeff = 0.0; } 
  else{ Eeff = sqrt(Eeff); }

  vector<double> sqE_Ei(nst, 0.0);
  for(int i=0; i<nst; i++){    
    sqE_Ei[i] = Etot - Ham.get(i,i).real(); 
    if(sqE_Ei[i]<0.0){ sqE_Ei[i] = 0.0; }
    else{ sqE_Ei[i] = sqrt(sqE_Ei[i]); }
    sqE_Ei[i] *= Eeff;
  }


  ij = 0;
  for(int i=0; i<nst; i++){
    for(int j=0; j<nst; j++){

      L.set(ij, ij, complex<double>(2.0*(sqE_Ei[j] - sqE_Ei[i]), 0.0) );

      ij++;
    }// for j
  }// for i

  return L;
  
}


MATRIX momenta_on_excited_states(dyn_variables& dyn_var, nHamiltonian* ham, int itraj){
/**
  Compute the momenta on all excited states given the momentum on the current active 
  state such that the total energy is conserved. All done for a given trajectory `itraj`.
  This is done in the adiabatic basis.
  
  The resulting momenta are taken to be parallel to the current momentum

  Returns: MATRIX(ndof, nadi)
*/

  int nst = dyn_var.nadi;
  int ndof = dyn_var.ndof;
  CMATRIX Ham(nst, nst);
  Ham = ham->get_ham_adi();

  int st_indx = dyn_var.act_states[itraj];
  double Ekin = dyn_var.compute_kinetic_energy(itraj);
  double Epot = Ham.get(st_indx, st_indx).real();

  MATRIX res(ndof, nst);
  MATRIX p_i(ndof, 1);
  MATRIX p(ndof, 1); p = dyn_var.p->col(itraj);

  for(int i=0; i<nst; i++){

    if(Ekin > 0.0){
      double Epot_i = Ham.get(i, i).real();
      double Ekin_i = Ekin + Epot - Epot_i; 
      if(Ekin_i<0.0){ Ekin_i = 0.0; }

      p_i = p * sqrt(Ekin_i / Ekin);
    }
    else{ p_i = 0.0; }
 
    for(int j=0; j<ndof; j++){  res.set(j, i, p_i.get(j, 0) );  }

  }// for i - all states

  return res;
}

MATRIX momenta_on_excited_states(dyn_variables& dyn_var, nHamiltonian& ham, int itraj){
  return momenta_on_excited_states(dyn_var, &ham, itraj);

}

void SSY_correction(CMATRIX& Ham, dyn_variables& dyn_var, nHamiltonian* ham, int itraj){
/**
  This is a correction of a Hamiltonian according to the phase-correction of Shenvi-Subotnik-Yang, 2011

  See my Chapter, Eq. 3.27
*/

  int ndof = dyn_var.ndof;
  int nst = dyn_var.nadi;
  
  MATRIX p(ndof, nst); 
  p = momenta_on_excited_states(dyn_var, ham, itraj);
  int st_indx = dyn_var.act_states[itraj]; // active state index
  MATRIX p_act(ndof, 1); p_act = p.col(st_indx);
  MATRIX p_tmp(ndof, 1);

   
  for(int i=0; i<nst; i++){
    p_tmp = p.col(i);
    p_tmp.dot_product( p_tmp,  *dyn_var.iM);
    p_tmp.dot_product( p_tmp, p_act);   
    double sm = p_tmp.sum(); 

    Ham.set(i, i, complex<double>(-sm, 0.0) );

  }// for i

}

void SSY_correction(CMATRIX& Ham, dyn_variables& dyn_var, nHamiltonian& ham, int itraj){
  SSY_correction(Ham, dyn_var, &ham, itraj);
}


void propagate_electronic(dyn_variables& dyn_var, nHamiltonian& ham, nHamiltonian& ham_prev, dyn_control_params& prms){

  propagate_electronic(dyn_var, &ham, &ham_prev, prms);

}

void apply_thermal_correction(dyn_variables& dyn_var, nHamiltonian& ham, nHamiltonian& ham_old, 
                              vector<int> old_states, dyn_control_params& prms, Random& rnd){
/**
  Computes the thermal corrections the adiabatic NACs according to the experimental method
  and scales NACs according to it
*/
  
  //======= Parameters of the dyn variables ==========
  int ndof = dyn_var.ndof;
  int ntraj = dyn_var.ntraj;
  int nadi = dyn_var.nadi;
  int ndia = dyn_var.ndia;

  int nst, a, b;
  if(prms.rep_tdse==0 || prms.rep_tdse==2 ){ nst = ndia; }
  else if(prms.rep_tdse==1 || prms.rep_tdse==3 || prms.rep_tdse == 4 ){ nst = nadi; }
  else{ nst = -1; } // to throw an error


  for(int itraj=0; itraj<ntraj; itraj++){
  
    int i = dyn_var.act_states[itraj];
    double E_tot = prms.total_energy; 
    double T0 = ham.children[itraj]->gs_kinetic_energy;
    double Ei = ham.children[itraj]->get_ham_adi().get(i,i).real();

    // the first call of the thermal correction function after initialization
    if( fabs(dyn_var.tcnbra_ekin[itraj] - (-1000.0)) < 1e-4){  
      dyn_var.tcnbra_ekin[itraj] = (E_tot - Ei );  
      dyn_var.tcnbra_thermostats[itraj].nu_therm = prms.tcnbra_nu_therm;
      dyn_var.tcnbra_thermostats[itraj].NHC_size = prms.tcnbra_nhc_size;

      dyn_var.tcnbra_thermostats[itraj].Nf_t = dyn_var.ndof - prms.constrained_dofs.size(); //For now, we assume  prms.quantum_dofs.size(); 

      dyn_var.tcnbra_thermostats[itraj].init_nhc();
    }
   

      int prev_i = old_states[itraj];
      double E_prev_i = ham_old.children[itraj]->get_ham_adi().get(prev_i, prev_i).real();
      double dE = Ei - E_prev_i; // adiabatic change of potential energy
      dyn_var.tcnbra_ekin[itraj] -= dE; // reflect it in the kinetic energy

      
    if(dyn_var.tcnbra_ekin[itraj] > 0.0){

      double scl = dyn_var.tcnbra_thermostats[itraj].vel_scale(0.5*prms.dt);
      dyn_var.tcnbra_ekin[itraj] *= (scl * scl);

      dyn_var.tcnbra_thermostats[itraj].propagate_nhc(prms.dt, dyn_var.tcnbra_ekin[itraj], 0.0, 0.0);  

      scl = dyn_var.tcnbra_thermostats[itraj].vel_scale(0.5*prms.dt);
      dyn_var.tcnbra_ekin[itraj] *= (scl * scl);
    }


    if(prms.tcnbra_do_nac_scaling==1){
      double alp = dyn_var.tcnbra_ekin[itraj] / T0;

      if(alp>0.0){  alp = sqrt(alp); }
      else{  alp = 1.0; }

      dyn_var.thermal_correction_factors[itraj] = alp;

      // Scale hvib_adi, and nac_adi by the alp parameter
      for(a=0; a<nst; a++){
        for(b=0; b<nst; b++){

          if(a!=b){
            ham.children[itraj]->nac_adi->scale(a, b, alp);
            ham.children[itraj]->hvib_adi->scale(a, b, alp);
            ham.children[itraj]->time_overlap_adi->scale(a, b, alp);
          }       
        }// for b
      }// for a
    }// if do scaling
    else{  dyn_var.thermal_correction_factors[itraj] = 1.0; }

  }// for itraj

}

void remove_thermal_correction(dyn_variables& dyn_var, nHamiltonian& ham, dyn_control_params& prms){

  int ntraj = dyn_var.ntraj;
  int nadi = dyn_var.nadi;
  int ndia = dyn_var.ndia;

  int nst, a,b;
  if(prms.rep_tdse==0 || prms.rep_tdse==2 ){ nst = ndia; }
  else if(prms.rep_tdse==1 || prms.rep_tdse==3 ){ nst = nadi; }

  for(int itraj=0; itraj<ntraj; itraj++){

    double alp = dyn_var.thermal_correction_factors[itraj];

    // Scale hvib_adi, and nac_adi by the alp parameter
    for(a=0; a<nst; a++){
      for(b=0; b<nst; b++){

        if(a!=b){
          ham.children[itraj]->nac_adi->scale(a, b, 1.0/alp);
          ham.children[itraj]->hvib_adi->scale(a, b, 1.0/alp);
          ham.children[itraj]->time_overlap_adi->scale(a, b, 1.0/alp);
//          ham->children[itraj]->ovlp_adi->scale(a, b, alp)
        }
      }// for b
    }// for a

  }// for itraj

}

void update_wp_width(dyn_variables& dyn_var, dyn_control_params& prms){
/**
  Updates the wave packet width in XF based on the Gaussian approximation
*/
  
  //======= Parameters of the dyn variables ==========
  int ndof = dyn_var.ndof;
  int ntraj = dyn_var.ntraj;

  if (prms.use_td_width == 0){
    for(int itraj=0; itraj<ntraj; itraj++){
      for(int idof=0; idof<ndof; idof++){
        dyn_var.wp_width->set(idof, itraj, prms.wp_width->get(idof, 0));
      }
    }
  }
  else if (prms.use_td_width == 1){
    double elapsed_time, s2;
    
    elapsed_time = prms.dt*dyn_var.timestep;

    for(int itraj=0; itraj<ntraj; itraj++){
      for(int idof=0; idof<ndof; idof++){
        s2 = pow(prms.wp_width->get(idof, 0), 2.0) + pow(prms.wp_v->get(idof, 0)*elapsed_time, 2.0);
        dyn_var.wp_width->set(idof, itraj, sqrt(s2));
      }
    }
  }
  else if (prms.use_td_width == 2){
    double elapsed_time, s2;

    //elapsed_time = prms.dt*dyn_var.timestep;
    
    for(int itraj=0; itraj<ntraj; itraj++){
      for(int idof=0; idof<ndof; idof++){
        s2 = pow(prms.wp_width->get(idof,0), 2.0);
        dyn_var.wp_width->set(idof, itraj, 4*M_PI/(s2*dyn_var.p->get(idof, itraj)) );
      }
    }
  }
  else if (prms.use_td_width == 3){
    int nadi = dyn_var.nadi;

    // wp_width is computed by pairwise widths based on auxiliary trajectories
    dyn_var.wp_width->set(-1, -1, 0.0);
    MATRIX sum_inv_w2(ndof, 1);
    MATRIX w2_temp(ndof, 1);

    for(int traj=0; traj<ntraj; traj++){
      vector<int>& is_mixed = dyn_var.is_mixed[traj];
      vector<int>& is_first = dyn_var.is_first[traj];
      
      sum_inv_w2.set(-1, 0, 0.0);
      w2_temp.set(-1, 0, 0.0);
    
      int check_mixing = 0;
      for(int i=0; i<nadi; i++){
        for(int j=0; j<nadi; j++){
          if(i>=j){continue;}

          if(is_mixed[i] == 0 or is_mixed[j] == 0){continue; }
          check_mixing = 1;

          if(is_first[i] == 1 or is_first[j] == 1){
            w2_temp.set(-1, 0, 1.0e+10); // At initial, an auxiliary pair does not contribute to wp_width
          }
          else{
            for(int idof=0; idof<ndof; idof++){
              w2_temp.set(idof, 0, fabs(dyn_var.q_aux[traj]->get(i, idof) - dyn_var.q_aux[traj]->get(j, idof))/
                fabs(dyn_var.p_aux[traj]->get(i, idof) - dyn_var.p_aux[traj]->get(j, idof)) );
            }
          }

          for(int idof=0; idof<ndof; idof++){
            sum_inv_w2.add(idof, 0, 1.0/w2_temp.get(idof));
          }

        } //j
      } //i
    
      if(check_mixing == 1){
        for(int idof=0; idof<ndof; idof++){
          dyn_var.wp_width->set(idof, traj, sqrt((nadi - 1)* 1.0/sum_inv_w2.get(idof, 0)) );
        }
      }

    } // traj
  }
}

int need_active_states_diff_rep(dyn_control_params& prms){
  int res = 0;

  int list_size = bp::len(prms.properties_to_save);

  if (prms.rep_sh == 0){
    for (int i = 0; i < list_size; ++i) {
        // Extract each element and convert it to a std::string
        std::string property = bp::extract<std::string>(prms.properties_to_save[i]);
        if (property == "sh_pop_adi") { res = 1;}
    }
  }
  else if(prms.rep_sh == 1){
    for (int i = 0; i < list_size; ++i) {
        // Extract each element and convert it to a std::string
        std::string property = bp::extract<std::string>(prms.properties_to_save[i]);
        if (property == "sh_pop_dia") { res = 1;}
    }
  }

  return res;
}

void propagate_electronic(dyn_variables& dyn_var, nHamiltonian* Ham, nHamiltonian* Ham_prev, dyn_control_params& prms){

  int itraj;

  int num_el = prms.num_electronic_substeps;
  double dt = prms.dt / num_el;
  int rep = prms.rep_tdse;
  int method = prms.electronic_integrator;
  int is_ssy = prms.do_ssy;

  //======= Parameters of the dyn variables ==========
  //int ndof = dyn_var.ndof;
  int ntraj = dyn_var.ntraj;
  int nadi = dyn_var.nadi;
  int ndia = dyn_var.ndia;

  int nst;
  if(prms.rep_tdse==0 || prms.rep_tdse==2 ){ nst = ndia; }
  else if(prms.rep_tdse==1 || prms.rep_tdse==3 || prms.rep_tdse == 4 ){ nst = nadi; }
  else{ nst = -1; } // to throw an error 

  //int ampl_transformation_method = prms.ampl_transformation_method;


  CMATRIX C(nst, 1);
  CMATRIX vRHO(nst*nst, 1); // vectorized DM
  CMATRIX RHO(nst, nst);    // DM
  CMATRIX L(nst*nst, nst*nst); // Liouvillian
  CMATRIX Coeff(nst, ntraj);
  if(prms.rep_tdse==0){ Coeff = *dyn_var.ampl_dia; }
  else if(prms.rep_tdse==1){ Coeff = *dyn_var.ampl_adi; }

//  dyn_var.compute_kinetic_energies()

  ///======================== Now do the integration of the TD-SE ===================
  for(itraj=0; itraj<ntraj; itraj++){

    C = Coeff.col(itraj);

    int traj1 = itraj;  if(method >=100 && method <200){ traj1 = 0; }

    nHamiltonian* ham = Ham->children[traj1];
    nHamiltonian* ham_prev = Ham_prev->children[traj1];
/*

    //================= Basis re-expansion ===================
    CMATRIX P(ham->nadi, ham->nadi);
    //proj[itraj]->load_identity();
    CMATRIX T(*dyn_var.proj_adi[itraj]);  T.load_identity();
    CMATRIX T_new(ham->nadi, ham->nadi);

    P = ham->get_time_overlap_adi();  // U_old.H() * U;

//    if(prms.thermally_corrected_nbra){
//      P = thermally_corrected_time_overlap(dyn_var, ham, prms, itraj);
//    }

    // More consistent with the new derivations:
    FullPivLU_inverse(P, T_new);
    T_new = orthogonalized_T( T_new );
    
    *dyn_var.proj_adi[itraj] = T_new;

*/
    CMATRIX T_new(*dyn_var.proj_adi[itraj]);
    if(prms.assume_always_consistent){    T_new.identity();   }


//    if ampl_transformation_method
     

  if(rep==0){  // diabatic
    CMATRIX Hvib(ham->ndia, ham->ndia);
    CMATRIX Sdia(ham->ndia, ham->ndia);

    if(method==-1){ ;;  } // No evolution
    else if(method==0 || method==100){
      // Based on Lowdin transformations, using mid-point Hvib
      Hvib = 0.5 * (ham->get_hvib_dia() + ham_prev->get_hvib_dia());
      Sdia = ham->get_ovlp_dia();
      propagate_electronic_eig(dt, C, Hvib, Sdia); // in this case C - diabatic coeffs
    }
    else if(method==1 || method==101){
      Hvib = 0.5 * (ham->get_hvib_dia() + ham_prev->get_hvib_dia());
      Sdia = ham->get_ovlp_dia();
      propagate_electronic_qtag(dt, C, Hvib, Sdia); // in this case C - diabatic coeffs
    }
    else if(method==2 || method==102){
      Hvib = ham->get_ham_dia();
      Sdia = ham->get_ovlp_dia();
      CMATRIX Hvib_old(ham->ndia, ham->ndia);   Hvib_old = ham_prev->get_ham_dia();
      CMATRIX Sdia_old(ham->ndia, ham->ndia);   Sdia_old = ham_prev->get_ovlp_dia();

      propagate_electronic_qtag2(dt, C, Hvib, Hvib_old, Sdia, Sdia_old);
    }
    else if(method==3 || method==103){
      // Using exp(S^-1 * Hvib_dia * dt)
      Hvib = 0.5 * (ham->get_hvib_dia() + ham_prev->get_hvib_dia());
      Sdia = ham->get_ovlp_dia();
      CMATRIX invS(ham->ndia, ham->ndia);
      FullPivLU_inverse(Sdia, invS);
      Hvib = invS * Hvib;
      propagate_electronic_nonHermitian(dt, C, Hvib);
    }

  }// rep == 0 // diabatic

  else if(rep==1){  // adiabatic
    /**
        -1              -  No

         0              -  ld, with crude splitting,  with exp_
         1              -  ld, with symmetric splitting, with exp_
         2              -  ld, original, with exp_       
         3              -  1-point, Hvib integration, with exp_
         4              -  2-points, Hvib integration, with exp_
         5              -  2-points, Hvib, integration with the second-point correction of Hvib, with exp_
         6              -  2-points, Hvib, no reordering

        10              -  same as 0, but with rotations
        11              -  same as 1, but with rotations
        12              -  same as 2, but with rotations
        13              -  same as 3, but with rotations
        14              -  same as 4, but with rotations
        15              -  same as 5, but with rotations

    */
    CMATRIX Hvib(ham->nadi, ham->nadi);
    CMATRIX Hvib_old(ham->nadi, ham->nadi);
    CMATRIX A(ham->nadi, ham->nadi); /// this is A(t)
    CMATRIX B(ham->nadi, ham->nadi); /// this is actually A(t+dt)

    if(method==-1){ ;;  } // No evolution
    else if(method==0 || method==100){
      // A crude factorization of the Hamiltonian operator
      Hvib = ham_prev->get_ham_adi();  // T is the identity matrix
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham_prev, itraj); }
      A = exp_(Hvib, complex<double>(0.0, -0.5*dt) );

      Hvib = ham->get_ham_adi();
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham, itraj); }
      B = exp_(Hvib, complex<double>(0.0, -0.5*dt) );

      C = B * T_new * A * C;
    }// method == 0 && 100

    else if(method==1 || method==101){
      // Trotter-based symmetric splitting
      Hvib = ham_prev->get_ham_adi();  // T is the identity matrix
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham_prev, itraj); }
      B = exp_(Hvib, complex<double>(0.0, -0.25*dt) );

      Hvib = ham->get_ham_adi();
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham, itraj); }
      A = exp_(Hvib, complex<double>(0.0, -0.5*dt) );

      C = T_new * B * T_new.H() * A * T_new * B * C;
    }// method == 1 && 101

    else if(method==2 || method==102){
      // The local diabatization approach
      Hvib_old = ham_prev->get_ham_adi();  // T is the identity matrix      
      if(is_ssy){ SSY_correction(Hvib_old, dyn_var, ham_prev, itraj); }

      Hvib = ham->get_ham_adi();
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham, itraj); }
      Hvib = T_new.H() * Hvib * T_new;      
      Hvib += Hvib_old;

      A = exp_(Hvib, complex<double>(0.0, -0.5*dt) );
      C = T_new * A * C;
//      C = A * C;

    }// method == 2 && 102

    else if(method==3 || method==103){
      // 1-point with vibronic Hamiltonian
      Hvib = ham_prev->get_hvib_adi();

      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham_prev, itraj); }

      A = exp_(Hvib, complex<double>(0.0, -dt) );
      C = T_new * A * C;
    }// method == 3 && 103

    else if(method==4 || method==104){
      // 2-point with vibronic Hamiltonian
      Hvib_old = ham_prev->get_hvib_adi();
      if(is_ssy){ SSY_correction(Hvib_old, dyn_var, ham_prev, itraj); }

      Hvib = ham->get_hvib_adi();
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham_prev, itraj); }

      Hvib += Hvib_old;

      A = exp_(Hvib, complex<double>(0.0, -0.5*dt) );
      C = T_new * A * C;
    }// method == 4 && 104

    else if(method==5 || method==105){
      // 2-point Hvib with vibronic Hamiltonian and correction
      Hvib_old = ham_prev->get_hvib_adi();  // T is the identity matrix
      if(is_ssy){ SSY_correction(Hvib_old, dyn_var, ham_prev, itraj); }

      Hvib = ham->get_hvib_adi();
      Hvib = T_new.H() * Hvib * T_new;
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham, itraj); }

      Hvib += Hvib_old;
      A = exp_(Hvib, complex<double>(0.0, -0.5*dt) );
      C = T_new * A * C;
    }// method == 5 && 105

    else if(method==6 || method==106){
      // 2-point with vibronic Hamiltonian, no reordering 
      Hvib_old = ham_prev->get_hvib_adi();
      if(is_ssy){ SSY_correction(Hvib_old, dyn_var, ham_prev, itraj); }

      Hvib = ham->get_hvib_adi();
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham_prev, itraj); }

      Hvib += Hvib_old;

      A = exp_(Hvib_old, complex<double>(0.0, -0.5*dt) );
      C = A * C;

      //dyn_var.proj_adi[itraj]->load_identity();

      //T_new.identity();
    }// method == 6 && 106

    else if(method==7 || method==107){
      // 2-point Hvib with vibronic Hamiltonian and correction - same as 5, but different order
      // Propagate old coefficients using old Hamiltonian, only half-way
      Hvib_old = ham_prev->get_hvib_adi();  // T is the identity matrix
      if(is_ssy){ SSY_correction(Hvib_old, dyn_var, ham_prev, itraj); }
      A = exp_(Hvib_old, complex<double>(0.0, -0.5*dt) );
      C = A * C;

      // Reproject the coefficients
      C = T_new.H() * C;

      // Propagate the new coefficients using the transformed new Hamiltonian, half-way
      Hvib = ham->get_hvib_adi();
      Hvib = T_new.H() * Hvib * T_new;
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham, itraj); }

      A = exp_(Hvib, complex<double>(0.0, -0.5*dt) );
      C = A * C;

    }// method == 7 && 107

    else if(method==8 || method==108){
      // new LD: 2-point Ham - same as 2 but different order

      // Propagate with the old Ham
      Hvib_old = ham_prev->get_hvib_adi();  // T is the identity matrix
      if(is_ssy){ SSY_correction(Hvib_old, dyn_var, ham_prev, itraj); }
      A = exp_(Hvib, complex<double>(0.0, -0.5*dt) );
      C = A * C;

      // Reorder:
      //C = T_new.H() * C;

      // Propagate with the new Ham
      Hvib = ham->get_hvib_adi();
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham, itraj); }
      Hvib = T_new.H() * Hvib * T_new;
      A = exp_(Hvib, complex<double>(0.0, -0.5*dt) );
      C = A * C;

      // Reorder back:
      //C = T_new * C;

    }// method == 8 && 108




    //========================== Rotation-based ============================
    else if(method==10 || method==110){
      // Same as 0 or 100, but with rotations

      Hvib = ham_prev->get_ham_adi();  // T is the identity matrix
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham_prev, itraj); }
      propagate_electronic_rot(0.5*dt, C, Hvib);
      C = T_new * C;

      Hvib = ham->get_ham_adi(); 
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham, itraj); }
      propagate_electronic_rot(0.5*dt, C, Hvib);

    }// method == 10 && 110

    else if(method==11 || method==111){
      // Trotter-based symmetric splitting, but with rotations
      Hvib_old = ham_prev->get_ham_adi();  // T is the identity matrix
      if(is_ssy){ SSY_correction(Hvib_old, dyn_var, ham_prev, itraj); }
      Hvib = ham->get_ham_adi();
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham, itraj); }

      propagate_electronic_rot(0.25*dt, C, Hvib_old);
      C = T_new * C;
      propagate_electronic_rot(0.5*dt, C, Hvib);
      C = T_new.H() * C;
      propagate_electronic_rot(0.25*dt, C, Hvib_old);
      C = T_new * C;

    }// method == 11 && 111

    else if(method==12 || method==12){
      // The local diabatization approach, but with rotation
      Hvib_old = ham_prev->get_ham_adi();  // T is the identity matrix
      if(is_ssy){ SSY_correction(Hvib_old, dyn_var, ham_prev, itraj); }

      Hvib = ham->get_ham_adi();
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham, itraj); }
      Hvib = T_new.H() * Hvib * T_new;
      Hvib += Hvib_old;

      propagate_electronic_rot(0.5*dt, C, Hvib);
      C = T_new * C;

    }// method == 12 && 112

    else if(method==13 || method==113){
      // 1-point with vibronic Hamiltonian
      Hvib = ham_prev->get_hvib_adi();
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham_prev, itraj); }

      propagate_electronic_rot(dt, C, Hvib);
      C = T_new * C;
    }// method == 13 && 113

    else if(method==14 || method==114){
      // 2-point with vibronic Hamiltonian
      Hvib_old = ham_prev->get_hvib_adi();
      if(is_ssy){ SSY_correction(Hvib_old, dyn_var, ham_prev, itraj); }

      Hvib = ham->get_hvib_adi();
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham_prev, itraj); }

      Hvib += Hvib_old;

      propagate_electronic_rot(0.5*dt, C, Hvib);
      C = T_new * C;
    }// method == 14 && 114

    else if(method==15 || method==115){
      // 2-point Hvib with vibronic Hamiltonian and correction
      Hvib_old = ham_prev->get_hvib_adi();  // T is the identity matrix
      if(is_ssy){ SSY_correction(Hvib_old, dyn_var, ham_prev, itraj); }

      Hvib = ham->get_hvib_adi();
      Hvib = T_new.H() * Hvib * T_new;
      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham, itraj); }

      Hvib += Hvib_old;

      propagate_electronic_rot(0.5*dt, C, Hvib);
      C = T_new * C;

    }// method == 15 && 115

  }// rep == 1 - adiabatic

  else if(rep==2){  // diabatic, density matrix formalism

    if(method==0 || method==100){
      // Based on Lowdin transformations, using mid-point Hvib

      CMATRIX Hvib(ham->ndia, ham->ndia);
      Hvib = 0.5 * (ham->get_hvib_dia() + ham_prev->get_hvib_dia());

      L = make_Liouvillian(Hvib);
      vRHO = vectorize_density_matrix( dyn_var.dm_dia[itraj] ); 

      vRHO = exp_(L, complex<double>(0.0, -dt) ) * vRHO;
//      propagate_electronic_rot(dt, vRHO, L);

      *dyn_var.dm_dia[itraj] = unvectorize_density_matrix( vRHO ); 

    }// method == 0 

  }// rep == 2 - diabatic, density matrix

  else if(rep==3){  // adiabatic, density matrix formalism
    /**
         0              -  mid-point Hvib with the second-point correction of Hvib
         1              -  Zhu Liouvillian

        10              -  same as 0, but with rotations
    */

    if(method==0 || method==100){
      // Based on Lowdin transformations, using mid-point Hvib
      CMATRIX Hvib(ham->nadi, ham->nadi);
      Hvib = 0.5 * (T_new.H() * ham->get_ham_adi() * T_new + ham_prev->get_ham_adi()); // "raw" to dyn-const

      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham_prev, itraj); }

      L = make_Liouvillian(Hvib);

      //RHO = T_new.H() * (*dyn_var.dm_adi[itraj]) * T; // "raw" to dyn-const
      RHO = *dyn_var.dm_adi[itraj];
      vRHO = vectorize_density_matrix( RHO );
      vRHO = exp_(L, complex<double>(0.0, -dt) ) * vRHO;

      RHO = unvectorize_density_matrix( vRHO );
      RHO = T_new * RHO * T_new.H(); // dyn-const to "raw"
      *dyn_var.dm_adi[itraj] = RHO;
    }// method == 0 or 100

    else if(method==1 || method==101){
      // Here goes the Zhu's method, Eqs. 3.34-3.35 of my Chapter
      // run this with a TSH versions, not Ehrenfest, cause the total energies
      // are different in these groups of methods

      CMATRIX Hvib(ham->nadi, ham->nadi);
      Hvib = ham_prev->get_ham_adi(); // + ham_prev->get_ham_adi());

      int st_indx = dyn_var.act_states[itraj];
      double Epot = Hvib.get(st_indx, st_indx).real(); 
      double Ekin = dyn_var.compute_kinetic_energy(itraj);

      L = Zhu_Liouvillian(Epot + Ekin, Hvib, *dyn_var.dm_adi[itraj] );
//      L = make_Liouvillian(Hvib);
      vRHO = vectorize_density_matrix( dyn_var.dm_adi[itraj] );
//      propagate_electronic_rot(dt, vRHO, L);
      vRHO = exp_(L, complex<double>(0.0, -dt) ) * vRHO;

      RHO = unvectorize_density_matrix( vRHO );
      RHO = T_new.H() * RHO * T_new;
      *dyn_var.dm_adi[itraj] = RHO; 

    }// method == 1 or 101

    else if(method==10 || method==110){
      // Same as 0, but with rotations
      CMATRIX Hvib(ham->nadi, ham->nadi);
      Hvib = 0.5 * (T_new.H() * ham->get_ham_adi() * T_new + ham_prev->get_ham_adi());  // "raw" to dyn-const

      if(is_ssy){ SSY_correction(Hvib, dyn_var, ham_prev, itraj); }

      L = make_Liouvillian(Hvib);
      //RHO = T_new.H() * (*dyn_var.dm_adi[itraj]) * T; // "raw" to dyn-const
      RHO = *dyn_var.dm_adi[itraj];
      vRHO = vectorize_density_matrix( RHO );

      propagate_electronic_rot(dt, vRHO, L);

      RHO = unvectorize_density_matrix( vRHO );
      RHO = T_new * RHO * T_new.H(); // dyn-const to "raw"
      *dyn_var.dm_adi[itraj] = RHO;
    }// method == 10 or 110


  }// rep == 3 - adiabatic, density matrix  



//  *dyn_var.proj_adi[itraj] = T_new;

  // Insert the propagated result back
  for(int st=0; st<nst; st++){  Coeff.set(st, itraj, C.get(st, 0));  }


  }// for itraj - all trajectories


  if(prms.rep_tdse==0){ *dyn_var.ampl_dia = Coeff; }
  else if(prms.rep_tdse==1){ *dyn_var.ampl_adi = Coeff; }


}

void renormalize_hopping_probabilities(vector< vector<double> >& g, vector< vector<vector<double> > >& coherence_factors, vector<int>& act_states){
/**
  Scale the hopping probabilities by the factors and ensure the sum of all 
  hopping probabilities is still 1

  g[itraj][i] - hopping probability from the current active state for the trajectory itraj to the state i

  act_states[itraj] - the active state for the trajectory itraj

  coherence_factors[itraj][i][j] - the scaling of the hopping probability from state i to state j for the trajectory itraj


  Result: the g vector is modified
*/
  int ntraj = g.size();
  int nstates = coherence_factors[0].size();

  for(int itraj=0; itraj<ntraj; itraj++){
    int a = act_states[itraj];
    double sum = 0.0;
    for(int j=0;j<nstates; j++){
      if(j!=a){ 
        double val = g[itraj][j] * coherence_factors[itraj][a][j];
        sum += val;
        g[itraj][j] =  val;
      } 
    }// for j
    // The self-element: 
    g[itraj][a] = 1.0 - sum;
    if(g[itraj][a]<0){  g[itraj][a] = 0.0; }

  }// for itraj  

}

void reset_coherence_factors(vector< vector<vector<double> > >& coherence_factors, vector<int>& act_states, vector<int>& old_states){

  int ntraj = coherence_factors.size();
  int nstates = coherence_factors[0].size();

  for(int itraj=0; itraj<ntraj; itraj++){
    if( act_states[itraj] != old_states[itraj]){
      // If the hop has happened for this trajectory, we start 
      // new evolution of the wavepackets on all the surfaces
      // so we have to reset all the coherence factors to 1.0
      for(int i=0; i<nstates; i++){
        for(int j=0; j<nstates; j++){
          coherence_factors[itraj][i][j] = 1.0;
        }// for j
      }// for i

    }// if
  }// for itraj

}


void compute_dynamics(dyn_variables& dyn_var, bp::dict dyn_params,
              nHamiltonian& ham, nHamiltonian& ham_aux, bp::object py_funct, bp::dict params,  Random& rnd,
              vector<Thermostat>& therm){
/**
  \brief One step of the TSH algorithm for electron-nuclear DOFs for one trajectory

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] C [nadi x ntraj]  or [ndia x ntraj] matrix containing the electronic coordinates. The amplitudes
   are assumed to be dynamically-consistent
  \param[in,out] projectors [ntraj CMATRIX(nadi, nadi)] - the projector matrices that account for the state tracking and 
  phase correction. These matrices should be considered as the dynamical varibles, similar to quantum amplitudes. Except
  their evolution does not necessarily follow from some equations of motion, but rather from various ad hoc schemes.
  \param[in,out] act_states - vector of ntraj indices of the physical states in which each of the trajectories
  initially is (active states). 
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) 
  - its internal variables (well, actually the variables it points to) are changed during the compuations
  \param[in] py_funct Python function object that is called when this algorithm is executed. The called Python function does the necessary 
  computations to update the diabatic Hamiltonian matrix (and derivatives), stored externally.
  \param[in] params The Python object containing any necessary parameters passed to the "py_funct" function when it is executed.
  \param[in] params1 The Python dictionary containing the control parameters passed to this function
  \param[in] rnd The Random number generator object

  Return: propagates C, q, p and updates state variables

*/

//  cout<<"In compute_dynamics\n";
  //======== General variables =======================
  int i, j, cdof, traj, dof, idof, ntraj1;

  //========= Control parameters variables ===========
  dyn_control_params prms;
  prms.set_parameters(dyn_params);

  int num_el = prms.num_electronic_substeps;

  //======= Parameters of the dyn variables ==========
  int ndof = dyn_var.ndof; 
  int ntraj = dyn_var.ntraj;   
  int nadi = dyn_var.nadi;
  int ndia = dyn_var.ndia;

  int nst; 
  if(prms.rep_tdse==0 || prms.rep_tdse==2 ){ nst = ndia; }
  else if(prms.rep_tdse==1 || prms.rep_tdse==3 || prms.rep_tdse == 4 ){ nst = nadi; }
  else{ nst = -1; } // to throw an error


  vector<int> act_states(dyn_var.act_states); // = dyn_var.act_states;
  vector<int> act_states_dia(dyn_var.act_states_dia);
  MATRIX p(*dyn_var.p);
  MATRIX& invM = *dyn_var.iM;
  

  //======== General variables ======================= 
  // Defining ntraj1 as a reference for making these matrices
  // ntraj is defined as q.n_cols as since it would be large in NBRA
  // we can define another variable like ntraj1 and build the matrices based on that.
  // We can make some changes where q is generated but this seems to be a bit easier
  if(prms.isNBRA==1){   ntraj1 = 1;  }
  else{  ntraj1 = ntraj;  }


  // Defining matrices based on ntraj1
  vector<CMATRIX> Eadi(ntraj1, CMATRIX(nst, nst));  
  vector<MATRIX> decoherence_rates(ntraj1, MATRIX(nst, nst)); 
  vector<double> Ekin(ntraj1, 0.0);  
  MATRIX gamma(ndof, ntraj);
  MATRIX p_traj(ndof, 1);
  vector<int> t1(ndof, 0); for(dof=0;dof<ndof;dof++){  t1[dof] = dof; }
  vector<int> t2(1,0);
  vector<int> t3(nst, 0); for(i=0;i<nst;i++){  t3[i] = i; }

  //============ Sanity checks ==================
  int n_therm_dofs = 1; // default 
  if(prms.ensemble==1){  
    n_therm_dofs = therm[0].Nf_t + therm[0].Nf_r;
    if(n_therm_dofs != prms.thermostat_dofs.size()){
      cout<<"Error in compute_dynamics: The number of thermostat DOFs ( currently "<<n_therm_dofs<<") must be \
      equal to the number of thermostat dofs set up by the `thermostat_dofs` parameter ( currently "
      <<prms.thermostat_dofs.size()<<")\nExiting...\n";
      exit(0);
    }
  }


  //***************************** Coherent dynamics *******************************
  //============== Nuclear propagation ===================
  // NVT dynamics
  if(prms.ensemble==1){  
    for(idof=0; idof<n_therm_dofs; idof++){
      dof = prms.thermostat_dofs[idof];
      for(traj=0; traj<ntraj; traj++){
        dyn_var.p->scale(dof, traj, therm[traj].vel_scale(0.5*prms.dt));
      }// traj
    }// idof 
  }

  *dyn_var.p = *dyn_var.p + 0.5 * prms.dt * (*dyn_var.f);
  //if(prms.use_qtsh==1){ *dyn_var.p = *dyn_var.p + 0.5 * compute_dkinemat(dyn_var, ham); }
  if(prms.use_qtsh==1){ *dyn_var.p = *dyn_var.p + 0.5 * prms.dt * (*dyn_var.qtsh_f_nc); }

  // Kinetic constraint
  for(cdof = 0; cdof < prms.constrained_dofs.size(); cdof++){   
    dyn_var.p->scale(prms.constrained_dofs[cdof], -1, 0.0); 
  }

  if(prms.entanglement_opt==22){
    gamma = ETHD3_friction(*dyn_var.q, *dyn_var.p, invM, prms.ETHD3_alpha, prms.ETHD3_beta);
  }
  // Update coordinates of nuclei for all trajectories
  for(traj=0; traj<ntraj; traj++){
    for(dof=0; dof<ndof; dof++){  
      dyn_var.q->add(dof, traj,  invM.get(dof,0) * dyn_var.p->get(dof,traj) * prms.dt ); 
      
      if(prms.entanglement_opt==22){
        dyn_var.q->add(dof, traj,  invM.get(dof,0) * gamma.get(dof,traj) * prms.dt ); 
      }
    }
  }



  // Use the current as the old
  ham_aux.copy_content(ham);

  // Recompute diabatic/adiabatic states, time-overlaps, NAC, Hvib, etc. in response to change of q
  update_Hamiltonian_variables(prms, dyn_var, ham, ham_aux, py_funct, params, 0);

  // Copy diabatic-to-adiabatic basis transformation to the dynamical variable
  dyn_var.update_basis_transform(ham);

  // Recompute the orthogonalized reprojection matrices, stored in dyn_var.proj_adi
  // this calculaitons used ham.children[i].time_overlap matrix, updated in the previous step
  update_proj_adi(prms, dyn_var, ham, ham_aux); 

  // Recompute NAC, Hvib, etc. in response to change of p
  update_Hamiltonian_variables(prms, dyn_var, ham, ham_aux, py_funct, params, 1);

  // Update wave packet width in XF algorithm
  if(prms.decoherence_algo == 5 or prms.decoherence_algo == 6){
    dyn_var.get_current_timestep(params);
    update_wp_width(dyn_var, prms);  
  }

  // Propagate electronic coefficients in the [t, t + dt] interval, this also updates the 
  // basis re-projection matrices 
  for(i=0; i<num_el; i++){
    if(prms.decoherence_algo == 5 or prms.decoherence_algo == 6){
      update_ham_xf(dyn_var);
      propagate_half_xf(dyn_var, ham, prms);
    }

    propagate_electronic(dyn_var, ham, ham_aux, prms);

    if(prms.decoherence_algo == 5 or prms.decoherence_algo == 6){
      rotate_nab_phase(dyn_var, ham, prms);
      update_ham_xf(dyn_var);
      propagate_half_xf(dyn_var, ham, prms);
    }
  }

  // Recompute density matrices in response to the updated amplitudes  
  dyn_var.update_amplitudes(prms);
  dyn_var.update_density_matrix(prms);
 
  vector<int> old_states( dyn_var.act_states);

  // In the interval [t, t + dt], we may have experienced the basis reordering, so we need to 
  // change the active adiabatic state
  if(prms.tsh_method == 3 or prms.tsh_method == 4 or prms.rep_sh == 0 ){  
  // Don't update states based on amplitudes, in the LZ method
    ;;
  }
  else{  dyn_var.update_active_states(1, 0); // 1 - forward; 0 - only active state
  }

  // For now, this function also accounts for the kinetic energy adjustments to reflect the adiabatic evolution
  if(prms.thermally_corrected_nbra==1){    apply_thermal_correction(dyn_var, ham, ham_aux, old_states, prms, rnd); }

  update_forces(prms, dyn_var, ham);
 
  if(prms.decoherence_algo == 6 and prms.use_xf_force == 1){
    update_forces_xf(dyn_var, ham, ham_aux);
  }

  // NVT dynamics
  if(prms.ensemble==1){  
    for(traj=0; traj<ntraj; traj++){
      t2[0] = traj; 
      pop_submatrix(p, p_traj, t1, t2);
      double ekin = compute_kinetic_energy(p_traj, invM, prms.thermostat_dofs);
      therm[traj].propagate_nhc(prms.dt, ekin, 0.0, 0.0);
    }

  }

  *dyn_var.p = *dyn_var.p + 0.5*prms.dt* (*dyn_var.f);
  //if(prms.use_qtsh==1){ *dyn_var.p = *dyn_var.p + 0.5* compute_dkinemat(dyn_var, ham); }
  if(prms.use_qtsh==1){ *dyn_var.p = *dyn_var.p + 0.5 * prms.dt * (*dyn_var.qtsh_f_nc); }

  // Kinetic constraint
  for(cdof=0; cdof<prms.constrained_dofs.size(); cdof++){   
    dyn_var.p->scale(prms.constrained_dofs[cdof], -1, 0.0); 
  }

  // NVT dynamics
  if(prms.ensemble==1){  
    for(idof=0; idof<n_therm_dofs; idof++){
      dof = prms.thermostat_dofs[idof];
      for(traj=0; traj<ntraj; traj++){
        dyn_var.p->scale(dof, traj, therm[traj].vel_scale(0.5*prms.dt));
      }// traj
    }// idof 
  }

  //ham_aux.copy_content(ham);
  update_Hamiltonian_variables(prms, dyn_var, ham, ham_aux, py_funct, params, 1);
  

  // Update amplitudes and density matrices in response to regular evolution
  dyn_var.update_amplitudes(prms);
  dyn_var.update_density_matrix(prms);
 
  //============== Begin the TSH part ===================

  //================= Update decoherence rates & times ================
  /// Effectively turn off decoherence effects
  if(prms.decoherence_times_type==-1){
    for(traj=0; traj<ntraj1; traj++){   decoherence_rates[traj] = 0.0;   }
  }
  /// Just use the plain times given from the input, usually the mSDM formalism
  else if(prms.decoherence_times_type==0){
    for(traj=0; traj<ntraj1; traj++){   decoherence_rates[traj] = *prms.decoherence_rates;   }
  }
  /// Compute the dephasing rates according the original energy-based formalism
  else if(prms.decoherence_times_type==1){
    Eadi = get_Eadi(ham);
    Ekin = dyn_var.compute_kinetic_energies();  
    decoherence_rates = edc_rates(Eadi, Ekin, prms.decoherence_C_param, prms.decoherence_eps_param, prms.isNBRA);       
  }

  else if(prms.decoherence_times_type==2){
    decoherence_rates = schwartz_1(prms, *dyn_var.ampl_adi, ham, *prms.schwartz_decoherence_inv_alpha); 
  }

  else if(prms.decoherence_times_type==3){
    decoherence_rates = schwartz_2(prms, ham, *prms.schwartz_decoherence_inv_alpha); 
  }

  else if(prms.decoherence_times_type==4){
    decoherence_rates = schwartz_1(prms, *dyn_var.ampl_adi, *dyn_var.p, ham, *prms.schwartz_interaction_width); 
  }

  else if(prms.decoherence_times_type==5){
    decoherence_rates = Gu_Franco(prms, *dyn_var.ampl_adi); 
  }

  ///== Optionally, apply the dephasing-informed correction ==
  if(prms.dephasing_informed==1){
    Eadi = get_Eadi(ham); 
    MATRIX ave_gaps(*prms.ave_gaps);
    dephasing_informed_correction(decoherence_rates, Eadi, ave_gaps, prms.isNBRA);
  }

  //============ Apply pre-TSH decoherence corrections ==================

  // SDM and alike methods - only in the adiabatic rep
  if(prms.decoherence_algo==0){
    if(prms.rep_tdse==1){
      *dyn_var.ampl_adi = sdm(*dyn_var.ampl_adi, prms.dt, dyn_var.act_states, 
                              decoherence_rates, prms.sdm_norm_tolerance, prms.isNBRA);
    }
    else{ cout<<"ERROR: SDM/mSDM requires rep_tdse = 1\nExiting now...\n"; exit(0); }
  }
  else if(prms.decoherence_algo==1){ /* ID-A. Nothing to do here, this is done in the surface hopping section */ }
  else if(prms.decoherence_algo==2){ /* A-FSSH. Nothing to do here, this is done in the surface hopping section */ }
  // BCSH
  else if(prms.decoherence_algo==3){ 
    if(prms.rep_tdse==1){
      wp_reversal_events(dyn_var, ham, prms.dt);
      *dyn_var.ampl_adi = bcsh(*dyn_var.ampl_adi, prms.dt, dyn_var.act_states, *dyn_var.reversal_events);
    }
    else{ cout<<"ERROR: BCSH requires rep_tdse = 1\nExiting now...\n"; exit(0); }
  }
  // MFSD
  else if(prms.decoherence_algo==4){
    if(prms.rep_tdse==1){
      p = *dyn_var.p;
      //cout<<"p before mfsd\n"; dyn_var.p->show_matrix();
      *dyn_var.ampl_adi = mfsd(p, *dyn_var.ampl_adi, *dyn_var.iM, prms.dt, dyn_var.act_states, decoherence_rates, ham, rnd, prms.isNBRA);
      *dyn_var.p = p;
      //cout<<"p after mfsd\n"; dyn_var.p->show_matrix();
       

      // Recompute NAC, Hvib, etc. in response to change of p
      update_Hamiltonian_variables(prms, dyn_var, ham, ham_aux, py_funct, params, 1);
    }
    else{ cout<<"ERROR: MFSD requires rep_tdse = 1\nExiting now...\n"; exit(0); }
  }
  // SHXF
  else if(prms.decoherence_algo==5){
    if(prms.rep_tdse==1){shxf(dyn_var, ham, ham_aux, prms);}
    else{ cout<<"ERROR: SHXF requires rep_tdse = 1\nExiting now...\n"; exit(0); }
  }
  // MQCXF
  else if(prms.decoherence_algo==6){
    if(prms.rep_tdse==1){mqcxf(dyn_var, ham, ham_aux, prms);}
    else{ cout<<"ERROR: MQCXF requires rep_tdse = 1\nExiting now...\n"; exit(0); }
  }

  // DISH, rev2023
  else if(prms.decoherence_algo==7){  dish_rev2023(dyn_var, ham,  decoherence_rates, prms, rnd);   }

  // Simple decoherence
  else if(prms.decoherence_algo==9){
    for(traj=0; traj<ntraj; traj++){
      for(i=0; i<nadi; i++){
        for(j=0; j<nadi; j++){
          double argg = decoherence_rates[traj].get(i,j) * prms.dt;
          dyn_var.coherence_factors[traj][i][j] *= exp( - argg * argg );
        }// for j
      }// for i
    }// for traj
  }// simple decoherence


  // Update amplitudes and density matrices in response to decoherence corrections
  dyn_var.update_amplitudes(prms);
  dyn_var.update_density_matrix(prms);


  //************************************ TSH options ****************************************
  // Adiabatic dynamics
  if(prms.tsh_method==-1){ ;; } 

  // FSSH, GFSH, MSSH, LZ, ZN, DISH, MASH, FSSH2, FSSH3
  else if(prms.tsh_method == 0 || prms.tsh_method == 1 || prms.tsh_method == 2 || prms.tsh_method == 3 
       || prms.tsh_method == 4 || prms.tsh_method == 5 || prms.tsh_method == 6 || prms.tsh_method == 7
       || prms.tsh_method == 8 || prms.tsh_method == 9 ){


    vector<int> old_states(dyn_var.act_states); 
    vector<int> old_states_dia(dyn_var.act_states_dia); 
    //========================== Hop proposal and acceptance ================================

    // FSSH (0), GFSH (1), MSSH (2), LZ(3), ZN (4), MASH(6), FSSH2(7), FSSH3(8), GFSH original (9)
    if(prms.tsh_method == 0 || prms.tsh_method == 1 || prms.tsh_method == 2 || prms.tsh_method == 3  
    || prms.tsh_method == 4 || prms.tsh_method == 6 || prms.tsh_method == 7 || prms.tsh_method == 8
    || prms.tsh_method == 9){

      /// Compute hop proposal probabilities from the active state of each trajectory to all other states 
      /// of that trajectory
      vector< vector<double> > g;
      g = hop_proposal_probabilities(prms, dyn_var, ham, ham_aux);

      // simple decoherence correction
      if(prms.decoherence_algo==9){   // intended only for adiabatic states for now 
        renormalize_hopping_probabilities(g, dyn_var.coherence_factors, dyn_var.act_states); 
      }

      // Propose new discrete states for all trajectories
      vector<int> prop_states(ntraj, 0);
      if(prms.rep_sh==1){
        prop_states = propose_hops(g, dyn_var.act_states, rnd);
    
        // Decide if to accept the transitions (and then which)        
        act_states = accept_hops(dyn_var, ham, prop_states, dyn_var.act_states, prms, rnd);
      }
      else if(prms.rep_sh==0){
        prop_states = propose_hops(g, dyn_var.act_states_dia, rnd);
        act_states_dia = accept_hops(dyn_var, ham, prop_states, dyn_var.act_states_dia, prms, rnd);
      }

      //=== Post-hop decoherence options ===

      // Instantaneous decoherence
      if(prms.decoherence_algo==1){
        if(prms.rep_tdse==1){
          instantaneous_decoherence(*dyn_var.ampl_adi, act_states, prop_states, old_states,
                                    prms.instantaneous_decoherence_variant, prms.collapse_option);
        }
        else{ cout<<"ERROR: Instantaneous Decoherence requires rep_tdse = 1\nExiting now...\n"; exit(0); }
      }
      else if(prms.decoherence_algo==2){
        /// Temporarily commented AVA 11/7/2022
        ///apply_afssh(dyn_var, Coeff, act_states, invM, ham, dyn_params, rnd);
      }// AFSSH
      else if(prms.decoherence_algo==3){ ;; } // BCSH of Linjun Wang, nothing to do here
      else if(prms.decoherence_algo==4){ ;; } // MF-SD of Bedard-Hearn, Larsen, Schwartz, nothing to do here 
      else if(prms.decoherence_algo==5 or prms.decoherence_algo==6){  
        if(prms.rep_sh==1){
          xf_hop_reset(dyn_var, act_states, old_states);   
        }
        else{ cout<<"ERROR: Independent-trajectory XF methods require rep_sh = 1\nExiting now...\n"; exit(0); }
      } // SHXF or MQCXF
      else if(prms.decoherence_algo==7){ ;; } // DISH rev 2023, nothing to do here

      // Experimental: instantaneous decoherence in diabatic basis | TODO: check the compatibility with rep_sh==0
      else if(prms.decoherence_algo==8){
        if(prms.rep_tdse==1){
          instantaneous_decoherence_dia(*dyn_var.ampl_adi, ham, act_states, prop_states, old_states,
                                    prms.instantaneous_decoherence_variant, prms.collapse_option);
        }
        else{ cout<<"ERROR: Instantaneous Decoherence requires rep_tdse = 1\nExiting now...\n"; exit(0); }
      }// algo == 8

      else if(prms.decoherence_algo==9){ // simple decoherence method
        reset_coherence_factors(dyn_var.coherence_factors, act_states, old_states);   
      }// algo == 9 

    }
    // DISH - the old one
    else if(prms.tsh_method == 5){
      if(prms.decoherence_algo==-1){ ;; }
      else{ cout<<"ERROR: DISH method should be used only with decoherence_algo = -1\nExiting now...\n"; exit(0); }
      act_states = dish(dyn_var, ham, decoherence_rates, prms, rnd);
    }// DISH



    //====================== Momenta adjustment after successful/frustrated hops ===================
    // Velocity rescaling: however here we may be changing velocities
    if(prms.rep_sh==1){
      handle_hops_nuclear(dyn_var, ham, act_states, old_states, prms);
      dyn_var.act_states = act_states;
    }
    else{
      handle_hops_nuclear(dyn_var, ham, act_states_dia, old_states_dia, prms);
      dyn_var.act_states_dia = act_states_dia;
    }

    // Set the active states in the other representation through the tranformation matrix
    if(need_active_states_diff_rep(prms) == 1){ dyn_var.set_active_states_diff_rep(prms.rep_sh, rnd); }

    // Re-scale (back) couplings and time-overlaps, if the TC-NBRA was used
    if(prms.thermally_corrected_nbra==1 && prms.tcnbra_do_nac_scaling==1){  remove_thermal_correction(dyn_var, ham, prms);  }
    
    // Update vib Hamiltonian to reflect the change of the momentum
    update_Hamiltonian_variables(prms, dyn_var, ham, ham_aux, py_funct, params, 1); 

  }// tsh_method == 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

  else{   cout<<"tsh_method == "<<prms.tsh_method<<" is undefined.\nExiting...\n"; exit(0);  }


  // Update the amplitudes and DM, in response to state hopping and other changes in the TSH part
  // so that we have them consistent in the output
  dyn_var.update_amplitudes(prms);
  dyn_var.update_density_matrix(prms);


  // Saves the current density matrix into the previous - needed for FSSH2 and GFSH (original)
  dyn_var.save_curr_dm_into_prev();

  

}




}// namespace libdyn
}// liblibra

