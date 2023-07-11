/*********************************************************************************
* Copyright (C) 2018-2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_decoherence_time.cpp
  \brief The file implements the methods to compute dephasing rates and decoherence intervals
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{


MATRIX edc_rates(CMATRIX& Hvib, double Ekin, double C_param, double eps_param, int isNBRA){
/**
    This function computes the decoherence rates matrix used in the 
    energy-based decoherence scheme of Granucci-Persico and Truhlar
    Reference: Granucci, G.; Persico, M. J. Chem. Phys. 2007, 126, 134114
    
    \param[in]       Hvib  [ CMATRIX ] Vibronic Hamiltonian matrix. 
    \param[in]        Ekin [ float ] The classical kinetic energy of nuclei. Units = Ha
    \param[in]     C_param [ float ] The method parameter, typically set to 1.0 Ha
    \param[in]   eps_param [ float ] The method parameter, typically set to 0.1 Ha
    \param[in]      isNBRA [ int   ] The method for considering NBRA-type calculations
*/

  int i,j;
  int nst = Hvib.n_cols;
  MATRIX decoh_rates(nst, nst);

  for(i=0; i<nst; i++){
    for(j=0; j<nst; j++){

      double itau = fabs( Hvib.get(i,i).real() - Hvib.get(j,j).real()) / ( C_param + (eps_param/Ekin) );
      decoh_rates.set(i,j, itau);

    }
  }

  return decoh_rates;

}


MATRIX edc_rates(CMATRIX& Hvib, double Ekin, double C_param, double eps_param){
  int is_nbra = 0;
  return edc_rates(Hvib, Ekin, C_param, eps_param, is_nbra);
}



vector<MATRIX> edc_rates(vector<CMATRIX>& Hvib, vector<double>& Ekin, double C_param, double eps_param, int isNBRA){

  int ntraj = Hvib.size();
  int nst = Hvib[0].n_cols;
  // A new variable instead of ntraj
  int ntraj1;

  if(isNBRA==1){
  ntraj1 = 1;
  }
  else{
  ntraj1 = Hvib.size();
  }

  vector<MATRIX> res(ntraj1, MATRIX(nst, nst));
  if(isNBRA==1){
  if(Ekin.size()!=ntraj){
    cout<<"ERROR in edc_rates: the sizes of the input variables Hvib and Ekin are inconsistent\n";
    cout<<"Hvib.size() = "<<Hvib.size()<<"\n";
    cout<<"Ekin.size() = "<<Ekin.size()<<"\n";
    cout<<"exiting...\n";
    exit(0);
  }
  }
  for(int traj=0; traj<ntraj1; traj++){
    res[traj] = edc_rates(Hvib[traj], Ekin[traj], C_param, eps_param, isNBRA);
  }
  return res;

}

vector<MATRIX> edc_rates(vector<CMATRIX>& Hvib, vector<double>& Ekin, double C_param, double eps_param){
  int is_nbra = 0;
  return edc_rates(Hvib, Ekin, C_param, eps_param, is_nbra); 
}




void dephasing_informed_correction(MATRIX& decoh_rates, CMATRIX& Hvib, MATRIX& ave_gaps, int isNBRA){
/**
    This function computes the corrected dephasing rates  
    The correction is based on the instnataneous energy levels and on the average
    of absolute values of the energy gaps

    Reference: Sifain, A. E.; Wang, L.; Teritiak, S.; Prezhdo, O. V. J. Chem. Phys. 2019, 150, 194104
    
    \param[in, out]  decoh_rates [ MATRIX ] uncorrected decoherence rates [units: a.u.t.^-1]
    \param[in]             Hvib  [ CMATRIX ] Instantaneous vibronic Hamiltonian [units: Ha]
    \param[in]          ave_gaps [ MATRIX ] time-averaged module of the energy level gaps: <|E_i - E_j|>  [units: Ha]
    \param[in]            isNBRA [ int    ] The method for considering NBRA-type calculations

*/

  int i,j;
  int nst = Hvib.n_cols;

  for(i=0; i<nst; i++){
    for(j=0; j<nst; j++){

      double dE_ij = fabs( Hvib.get(i,i).real() - Hvib.get(j,j).real());

      if(ave_gaps.get(i,j) > 0.0){

        decoh_rates.scale(i,j, dE_ij / ave_gaps.get(i,j) );

      }
      else{

       decoh_rates.set(i,j, 1e+25);

      }

    }// for j
  }// for i

}


void dephasing_informed_correction(MATRIX& decoh_rates, CMATRIX& Hvib, MATRIX& ave_gaps){

  int is_nbra = 0;
  dephasing_informed_correction(decoh_rates, Hvib, ave_gaps, is_nbra);
}



void dephasing_informed_correction(vector<MATRIX>& decoh_rates, vector<CMATRIX>& Hvib, MATRIX& ave_gaps, int isNBRA){

  int ntraj = Hvib.size();

  if(isNBRA==1){
    dephasing_informed_correction(decoh_rates[0], Hvib[0], ave_gaps, isNBRA);
  }
  else{
  if(decoh_rates.size()!=ntraj){
    cout<<"ERROR in dephasing_informed_correction: the sizes of the input variables \
    decoh_rates and Hvib are inconsistent\n";
    cout<<"decoh_rates.size() = "<<decoh_rates.size()<<"\n";
    cout<<"Hvib.size() = "<<Hvib.size()<<"\n";
    cout<<"exiting...\n";
    exit(0);
  }

  for(int traj=0; traj<ntraj; traj++){

    dephasing_informed_correction(decoh_rates[traj], Hvib[traj], ave_gaps, isNBRA);

  }
  }
}


void dephasing_informed_correction(vector<MATRIX>& decoh_rates, vector<CMATRIX>& Hvib, MATRIX& ave_gaps){

  int is_nbra = 0;
  dephasing_informed_correction(decoh_rates, Hvib, ave_gaps, is_nbra);

}



MATRIX coherence_intervals(CMATRIX& Coeff, MATRIX& rates){
/**
  This function computes the time-dependent (and population-dependent) coherence intervals
  (the time after which different states should experience a decoherence event)
  as described by Eq. 11 in:
  Jaeger, H. M.; Fischer, S.; Prezhdo, O. V. Decoherence-Induced Surface Hopping. J. Chem. Phys. 2012, 137, 22A545.

  1/tau_i  (t) =  sum_(j!=i)^nstates {  rho_ii(t) * rate_ij }


  \param[in] Coeff Amplitudes of the electronic states
  \param[in] rates A matrix containing the decoherence rates (inverse of the
  decoherence time for each given pair of states)

  Returns: A matrix of the coherence intervals for each state

*/
  int nstates = Coeff.n_rows; 

  CMATRIX denmat(nstates, nstates);   
  //denmat = (Coeff * Coeff.H() ).conj();
  denmat = Coeff * Coeff.H();

  MATRIX tau_m(nstates, 1);   

  for(int i=0;i<nstates;i++){

    double summ = 0.0;
    for(int j=0;j<nstates;j++){

      if(j!=i){
        summ += denmat.get(j,j).real() * rates.get(i,j); 
      }// if

    }// for j

    if(summ>0.0){   tau_m.set(i, 0, 1.0/summ); }
    else        {   tau_m.set(i, 0, 1.0e+25);  } // infinite coherence interval
    
     
  }// for i

//  delete denmat;

  return tau_m;
}


MATRIX coherence_intervals(CMATRIX& Coeff, vector<MATRIX>& rates){
/**
  This function computes the time-dependent (and population-dependent) coherence intervals
  (the time after which different states should experience a decoherence event)
  as described by Eq. 11 in:
  Jaeger, H. M.; Fischer, S.; Prezhdo, O. V. Decoherence-Induced Surface Hopping. J. Chem. Phys. 2012, 137, 22A545.

  1/tau_i  (t) =  sum_(j!=i)^nstates {  rho_ii(t) * rate_ij }


  \param[in] Coeff - CMATRIX(nstates, ntraj) Amplitudes of the electronic states
  \param[in] rates - ntraj matrices of nstates x nstates - Matrices containing the 
  decoherence rates (inverse of the decoherence time for each given pair of states) for each trajectory

  Returns: A matrix of the coherence intervals for each state for each trajectory

*/
  int i, traj;
  int nstates = Coeff.n_rows; 
  int ntraj = Coeff.n_cols;

  if(ntraj!=rates.size()){
    cout<<"ERROR in coherence_intervals: the ntraj dimensions do not agree for Coeff and rates\n";
    cout<<"ntraj = "<<ntraj<<endl;
    cout<<"rates.size() = "<<rates.size()<<endl;
    cout<<"Exiting...\n";
    exit(0);
  }

  MATRIX res(nstates, ntraj);
  CMATRIX coeff(nstates, 1);
  MATRIX tau_m(nstates, 1);

  vector<int> stenc_x(nstates, 0); for(i=0;i<nstates;i++){  stenc_x[i] = i; }
  vector<int> stenc_y(1, 0); 

  for(traj=0; traj<ntraj; traj++){

    stenc_y[0] = traj;
    pop_submatrix(Coeff, coeff, stenc_x, stenc_y);

    tau_m = coherence_intervals(coeff, rates[traj]);

    push_submatrix(res, tau_m, stenc_x, stenc_y);

  }

  return res;
}


vector<MATRIX> schwartz_1(dyn_control_params& prms, CMATRIX& amplitudes, nHamiltonian& ham, MATRIX& inv_alp){
/**
  Compute decoherence rates 1/tau_i for all states and all trajectories according to Schwartz prescription 
  
  amplitudes  - CMATRIX(nstates, ntraj)
  inv_alp - MATRIX(ndof, 1)

  Return:

  MATRIX(nstates, ntraj) - 1/tau - decoherence rates for all states and trajectories
*/

  int ndof = ham.nnucl;
  int nstates = ham.nadi; 
  int ntraj = ham.children.size();


  MATRIX F_mf(ndof, ntraj);
  MATRIX F_st(ndof, ntraj);
  MATRIX tmp(ndof, ntraj);

  vector<MATRIX> res(ntraj, MATRIX(nstates, nstates));

  vector<int> act_states(ntraj, 0);
  /// AVA - commented for now, 12/7/2022
  ///F_mf = ham.Ehrenfest_forces_adi(amplitudes, 1).real();  //aux_get_forces(prms_mf, amplitudes, projectors, act_states, ham);

  int option = 0; // default value for NAC-based integrators
  if(prms.electronic_integrator==0 ||  prms.electronic_integrator==1 ||
     prms.electronic_integrator==2 ||  prms.electronic_integrator==10 ||
     prms.electronic_integrator==11 || prms.electronic_integrator==12 
    ){ option = 1; }

  F_mf = ham.Ehrenfest_forces_adi(amplitudes, 1, option).real();

  for(int i=0;i<nstates; i++){
    vector<int> act_states(ntraj, i);
    
    F_st = ham.forces_adi(act_states).real();  // aux_get_forces(prms_mf, amplitudes, projectors, act_states, ham);

    for(int itraj=0; itraj<ntraj; itraj++){

      double tau_inv2 = 0.0;
      for(int idof=0; idof<ndof; idof++){
        double dF = F_mf.get(idof, itraj) - F_st.get(idof, itraj);
        
        tau_inv2 += 0.25 * inv_alp.get(idof, 0) * dF * dF; 
      }

      double tau_inv = sqrt(tau_inv2);

      res[itraj].set(i, i, tau_inv);

    }// for itraj
    
  }// for i

  return res;
}




vector<MATRIX> schwartz_2(dyn_control_params& prms, nHamiltonian& ham, MATRIX& inv_alp){
/**
  Compute decoherence rates 1/tau_ij for all pairs of states and all trajectories according to Schwartz state-pair prescription 
  
  inv_alp - MATRIX(ndof, 1)

  Return:

  MATRIX(nstates, nstates) x ntraj  - 1/tau_ij - decoherence rates for all pairs of states and trajectories
*/


  int ndof = ham.nnucl;
  int nstates = ham.nadi; 
  int ntraj = ham.children.size();

//  dyn_control_params prms_st(prms);  prms_st.force_method = 1;  prms_st.rep_force = 1; /// adiabatic, state-resolved force


  // Precompute state-resolved forces
  //CMATRIX amplitudes(nstates, ntraj);
  vector<MATRIX> F(nstates, MATRIX(ndof, ntraj));

  for(int i=0; i<nstates; i++){
    vector<int> act_states_i(ntraj, i);
    F[i] = ham.forces_adi(act_states_i).real(); //  aux_get_forces(prms_st, amplitudes, projectors, act_states_i, ham);

  }// for i


  vector<MATRIX> res(ntraj, MATRIX(nstates, nstates));

  for(int i=0; i<nstates; i++){
    for(int j=i+1; j<nstates; j++){
    
      for(int itraj=0; itraj<ntraj; itraj++){

        double tau_inv2 = 0.0;
        for(int idof=0; idof<ndof; idof++){
          double dF = F[i].get(idof, itraj) - F[j].get(idof, itraj);
        
          tau_inv2 += 0.25 * inv_alp.get(idof, 0) * dF * dF; 
        }
        double tau_inv = sqrt(tau_inv2);
        res[itraj].set(i, j, tau_inv);
        res[itraj].set(j, i, tau_inv);

      }// for itraj    
    }// for j
  }// for i

  return res;

}



}// namespace libdyn
}// liblibra

