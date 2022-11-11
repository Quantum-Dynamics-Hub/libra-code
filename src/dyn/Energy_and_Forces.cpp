/*********************************************************************************
* Copyright (C) 2015-2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Energy_and_Forces.cpp
  \brief The file implements the functions for dynamics-immediate energy calculations

  The "dynamics-immediate" means the energies and forces computed and organized to 
  be used in dynamical calculations of different type - classical, quantum, quantum-classical.
    
*/

#include "Energy_and_Forces.h"
#include "dyn_projectors.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace 
namespace libdyn{


double compute_kinetic_energy(MATRIX& p, MATRIX& invM, vector<int>& which_dofs){
/**
  \brief Compute kinetic energy of nuclear DOFs
  \param[in] p [ndof x ntraj] Momenta of ntraj replicas of the system with ndof 
  nuclear DOFs
  \param[in] invM [ndof x 1] Matrix of inverted masses of all DOFs

  In the case of more than 1 trajectories supplied into this function, the 
  total kinetic energy of all trajectories is returned.

  This is the classical nuclear kinetic energy
*/
  int ndof = which_dofs.size(); ///p.n_rows;
  int ntraj = p.n_cols;

  double Ekin = 0.0;
  
    
  for(int idof=0; idof < ndof; idof++){ 

    int dof = which_dofs[idof];

    double sum = 0.0;
    for(int traj=0; traj < ntraj; traj++){    
      sum = p.get(dof, traj) * p.get(dof, traj);
    }
    Ekin += sum * invM.get(dof, 0);
  }
  Ekin *= 0.5;

  return Ekin;

}

double compute_kinetic_energy(MATRIX& p, MATRIX& invM){

  int ndof = p.n_rows;
  vector<int> which_dofs(ndof);
  for(int i = 0; i < ndof; i++){  which_dofs[i] = i; }

  return compute_kinetic_energy(p, invM, which_dofs);

}


vector<double> compute_kinetic_energies(MATRIX& p, MATRIX& invM, vector<int>& which_dofs){
/**
  \brief Compute trajectory-resolved kinetic energies of nuclear DOFs

  \param[in] p [ndof x ntraj] Momenta of ntraj replicas of the system with ndof 
  nuclear DOFs
  \param[in] invM [ndof x 1] Matrix of inverted masses of all DOFs

  This is the classical nuclear kinetic energy
*/
  int ndof = which_dofs.size(); ///p.n_rows;
  int ntraj = p.n_cols;

  vector<double> Ekin(ntraj, 0.0);
  
  for(int traj=0; traj < ntraj; traj++){    

    double sum = 0.0;

    for(int idof=0; idof < ndof; idof++){ 

      int dof = which_dofs[idof];
      sum += p.get(dof, traj) * invM.get(dof, 0) * p.get(dof, traj);
    }
    sum *= 0.5;

    Ekin[traj] = sum;

  }// for traj

  return Ekin;

}



vector<double> compute_kinetic_energies(MATRIX& p, MATRIX& invM){

  int ndof = p.n_rows;
  vector<int> which_dofs(ndof);
  for(int i = 0; i < ndof; i++){  which_dofs[i] = i; }

  return compute_kinetic_energies(p, invM, which_dofs);

}





CMATRIX tsh_indx2ampl(vector<int>& res, int nstates){
/**
  Convert the active state index to the vector with occupation numbers
  This is done for 1 or many trajectories
*/

  int ntraj = res.size();
  CMATRIX ampl(nstates, ntraj);

  
  for(int traj = 0; traj < ntraj; traj++){

    ampl.set(res[traj], traj, complex<double>(1.0, 0.0));

  }

  return ampl;

}


double average_potential_energy(dyn_control_params& prms, dyn_variables& dyn_vars, nHamiltonian& ham){

  vector<double> res(dyn_vars.ntraj, 0.0);
  res = potential_energies(prms, dyn_vars, ham);

  double ave = 0.0;
  for(int itraj = 0; itraj < dyn_vars.ntraj; itraj++){  ave += res[itraj]; }
  ave /= double(dyn_vars.ntraj);
  
  return ave;
}

double average_potential_energy(bp::dict params, dyn_variables& dyn_vars, nHamiltonian& ham){

  dyn_control_params prms;
  prms.set_parameters(params);

  return average_potential_energy(prms, dyn_vars, ham);
}


vector<double> potential_energies(dyn_control_params& prms, dyn_variables& dyn_vars, nHamiltonian& ham){

  int itraj;
  int ntraj = dyn_vars.ntraj;
  vector<int> id(2,0);

  vector<double> res(ntraj, 0.0);

  if(prms.force_method==0){  // No forces

    // Don't compute forces at all - e.g. in NBRA
    for(itraj=0; itraj<ntraj; itraj++){ res[itraj] = 0.0; }

  }// NBRA

  else if(prms.force_method==1){   // State-specific forces

    // TSH or adiabatic (including excited states)
    // state-specific forces

    vector<int> effective_states(dyn_vars.act_states);

    if(prms.enforce_state_following==1){ 
      for(itraj=0; itraj<ntraj; itraj++){ effective_states[itraj] = prms.enforced_state_index; }
    }

    // Diabatic 
    if(prms.rep_force==0){  
      for(itraj=0; itraj<ntraj; itraj++){
        id[1] = itraj;
        int ist = effective_states[itraj];
        res[itraj] = ham.get_ham_dia(id).get(ist, ist).real();   
      }
    }
    // Adiabatic 
    else if(prms.rep_force==1){  
      for(itraj=0; itraj<ntraj; itraj++){
        id[1] = itraj;
        int ist = effective_states[itraj];
        res[itraj] = ham.get_ham_adi(id).get(ist, ist).real();   
      }
    }

  }// TSH && adiabatic

  else if(prms.force_method==2){  // Ehrenfest forces

    // Diabatic 
    if(prms.rep_force==0){  
      for(itraj=0; itraj<ntraj; itraj++){
        id[1] = itraj;
        res[itraj] = ham.Ehrenfest_energy_dia(*dyn_vars.ampl_dia, id).real();
      }
    }
    // Adiabatic 
    else if(prms.rep_force==1){  
      for(itraj=0; itraj<ntraj; itraj++){
        id[1] = itraj;
        res[itraj] = ham.Ehrenfest_energy_adi(*dyn_vars.ampl_adi, id).real();
      }
    }
  
  }// Ehrenfest


  return res;
}

vector<double> potential_energies(bp::dict params, dyn_variables& dyn_vars, nHamiltonian& ham){

  dyn_control_params prms;
  prms.set_parameters(params);

  return potential_energies(prms, dyn_vars, ham);

}

MATRIX aux_get_forces(dyn_control_params& prms, dyn_variables& dyn_vars, nHamiltonian& ham){
  /**
    Compute the force depending on the method used
  */

  int ndof = ham.nnucl;
  int nst = ham.nadi;
  int ntraj = dyn_vars.ntraj;

  MATRIX F(ndof, ntraj);
//  CMATRIX _amplitudes(nst, ntraj); // CMATRIX version of "act_states"

  if(prms.force_method==0){  // No forces

    // Don't compute forces at all - e.g. in NBRA

  }// NBRA

  else if(prms.force_method==1){   // State-specific forces

    // TSH or adiabatic (including excited states)
    // state-specific forces

    vector<int> effective_states(dyn_vars.act_states);

    if(prms.enforce_state_following==1){ // NBRA-like enforcement: adiabatic dynamics, in terms of forces 
      for(int itraj=0; itraj<ntraj; itraj++){ effective_states[itraj] = prms.enforced_state_index; }
    }

    // Diabatic 
    if(prms.rep_force==0){  F = ham.forces_dia(effective_states).real();   }
    // Adiabatic 
    else if(prms.rep_force==1){  F = ham.forces_adi(effective_states).real();  }

  }// TSH && adiabatic

  else if(prms.force_method==2){  // Ehrenfest forces

    // Diabatic 
    if(prms.rep_force==0){  
//      CMATRIX _amplitudes(nst, ntraj);
//      _amplitudes = *dynvars.ampl_dia;
      F = ham.Ehrenfest_forces_dia(*dyn_vars.ampl_dia, 1).real();
    }

    // Adiabatic 
    else if(prms.rep_force==1){  

      //CMATRIX _amplitudes(nst, ntraj); 
      //_amplitudes = dynconsyst_to_raw(amplitudes, projectors);
      //_amplitudes = *dynvars.ampl_adi;
      F = ham.Ehrenfest_forces_adi(*dyn_vars.ampl_adi, 1).real();
    }
  
  }// Ehrenfest


  return F;
}

MATRIX aux_get_forces(bp::dict params, dyn_variables& dyn_vars, nHamiltonian& ham){

  dyn_control_params prms;
  prms.set_parameters(params);

  return aux_get_forces(prms, dyn_vars, ham);
  
}



vector<CMATRIX> get_Eadi(nHamiltonian& ham){

  int nst = ham.nadi;
  int ntraj = ham.children.size();

  vector<CMATRIX> Eadi(ntraj, CMATRIX(nst, nst));

  for(int traj=0; traj<ntraj; traj++){
    Eadi[traj] = ham.children[traj]->get_ham_adi(); 
  }

  return Eadi;

}







double compute_kinetic_energy(Nuclear* mol){
/**
  \brief Compute kinetic energy of Nuclear object 
  \param[in] mol The pointer to the Nuclear DOF object

  This is the classical nuclear kinetic energy
*/

  double Ekin = 0.0;

  for(int i=0;i<mol->nnucl;i++){  Ekin += mol->p[i]*mol->p[i]/mol->mass[i];   }
  Ekin *= 0.5;

  return Ekin;

}

double compute_kinetic_energy(Nuclear& mol){  
/**
  \brief Compute kinetic energy of Nuclear object - Python-friendly
  \param[in] mol The reference to the Nuclear DOF object

  This is the classical nuclear kinetic energy
*/

  return compute_kinetic_energy(&mol); 
}

double compute_kinetic_energy(Ensemble& ens){
/**
  \brief Compute average kinetic energy of an ensemble of Nuclear objects - Python-friendly
  \param[in] ens The reference to the Ensemble object

  This is the classical nuclear kinetic energy
*/


  double res = 0.0;
  for(int traj=0;traj<ens.ntraj;traj++){
    res += compute_kinetic_energy(&ens.mol[traj]);
  }
  res /= (double)ens.ntraj;
  
  return res;
}



double compute_potential_energy(Nuclear* mol, Electronic* el, Hamiltonian* ham, int opt){
/**
  \brief Compute potential energy of a quantum-classical system
  \param[in,out] mol Describes the nuclear DOF
  \param[in] el Describes electronic DOF
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0
*/

  int i,j;
  double Heff = 0.0;
  double Epot = 0.0;


  // Calculate all surfaces,
  ham->set_q(mol->q);
  ham->compute();


  // Mean-field/Ehrenfest mixing
  if(opt==0){

    // Potential energy - electronic Hamiltonian
    // working in a.u., so hbar = 1
    for(i=0;i<el->nstates;i++){
      for(j=0;j<el->nstates;j++){
        Heff += 0.5*ham->Hvib(i,j).real() * (el->q[i] * el->q[j] + el->p[i] * el->p[j]);
        Heff += ham->Hvib(i,j).imag() * el->p[i] * el->q[j]; 
      }// for j
    }// for i

    Epot = 2.0*Heff;  // in a.u. - hbar is assumed = 1
  }// algo = "MF"

  // FSSH mixing 
  else if(opt==1){   
    i = el->istate; // current electornic state
    Epot = ham->H(i,i).real();  
  }

  return Epot; 

}

double compute_potential_energy(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt){
/**
  \brief Compute potential energy of a quantum-classical system - Python-friendly
  \param[in,out] mol Describes the nuclear DOF
  \param[in] el Describes electronic DOF
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0

*/

  return compute_potential_energy(&mol, &el, &ham, opt);

}

double compute_potential_energy(Ensemble& ens, int opt){
/**
  \brief Compute potential energy of a quantum-classical system - Python-friendly
  \param[in,out] ens The reference to the ensemble of trajectories
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0

  The returned potential energy is also averaged over all trajectories included in the ensemble object ens.

*/

  double res = 0.0;
  for(int traj=0;traj<ens.ntraj;traj++){
    res += compute_potential_energy(&ens.mol[traj], &ens.el[traj], ens.ham[traj], opt);
  }
  res /= (double)ens.ntraj;
  
  return res;
}



double compute_forces(Nuclear* mol, Electronic* el, Hamiltonian* ham, int opt){
/**
  \brief Compute potential energy and forces (an all classical DOF) of a quantum-classical system
  \param[in,out] mol Describes the nuclear DOF, is modified during calculations to update forces
  \param[in] el Describes electronic DOF
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0
  The computed forces are stored in the mol object
*/


  int i,j,k;
  double Heff = 0.0;
  double Epot = 0.0;


  // Calculate all surfaces, if needed
  ham->set_q(mol->q);
  ham->compute();

  // Zero all forces in mol
  for(i=0;i<mol->nnucl;i++){ mol->f[k] = 0.0; }

  
  // Mean-field/Ehrenfest mixing
  if(opt==0){

    // Potential energy - electronic Hamiltonian
    // working in a.u., so hbar = 1
    for(i=0;i<el->nstates;i++){
      for(j=0;j<el->nstates;j++){

        double cij_re = (el->q[i] * el->q[j] + el->p[i] * el->p[j]);
        double cij_im = el->p[i] * el->q[j];

        Heff += 0.5*ham->Hvib(i,j).real() * cij_re;
        Heff += ham->Hvib(i,j).imag() * cij_im; 

        for(k=0;k<mol->nnucl;k++){
          
          mol->f[k] -= 2.0 * 0.5*ham->dHdq(i,j,k).real() * cij_re;
          mol->f[k] -= 2.0 * ham->dHdq(i,j,k).imag() * cij_im; 

        }// for k

      }// for j
    }// for i

    Epot = 2.0*Heff;  // in a.u. - hbar is assumed = 1

  }// algo = "MF"

  // FSSH mixing 
  else if(opt==1){

    i = el->istate; // current electornic state
    Epot = ham->H(i,i).real();  

    for(k=0;k<mol->nnucl;k++){  
      mol->f[k] = -ham->dHdq(i,i,k).real(); 
    }// for k

  }// opt == 1

  return Epot;

}

double compute_forces(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt){
/**
  \brief Compute potential energy and forces (an all classical DOF) of a quantum-classical system - Python-friendly
  \param[in,out] mol Describes the nuclear DOF, is modified during calculations to update forces
  \param[in] el Describes electronic DOF
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0
  The computed forces are stored in the mol object
*/

  return compute_forces(&mol, &el, &ham, opt);

}

double compute_forces(Ensemble& ens, int opt){
/**
  \brief Compute potential energy and forces (an all classical DOF) of a quantum-classical system
  \param[in,out] ens The reference to the ensemble of trajectories
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0
  The computed forces are stored in the mol object

  The returned potential energy is also averaged over all trajectories included in the ensemble object ens.
  No averaging of forces is needed.

*/


  double epot = 0.0;

  for(int traj=0;traj<ens.ntraj;traj++){
    epot += compute_forces(&ens.mol[traj], &ens.el[traj], ens.ham[traj], opt);
  }
  epot = epot/float(ens.ntraj);

  return epot;
  
}



void compute_energies(Ensemble* ens, double& Epot, double& Ekin, double& Etot,int opt){
/**
  \brief Compute ensemble-averaged potential, kinetic and total energy
  \param[in,out] ens The pointer to Ensemble object for which we want to compute properties
  \param[out] Epot The computed averaged potential energy will be stored here
  \param[out] Ekin The computed averaged kinetic energy will be stored here
  \param[out] Etot The computed averaged total energy will be stored here
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  
*/
  Epot = 0.0;
  Ekin = 0.0;
  Etot = 0.0;
  
  for(int traj=0;traj<ens->ntraj;traj++){

    double ek = compute_kinetic_energy(ens->mol[traj]);
    double ep = compute_potential_energy(&ens->mol[traj], &ens->el[traj], ens->ham[traj], opt);

    Epot += ep;
    Ekin += ek;
  } 

  Epot /= ((double)ens->ntraj);
  Ekin /= ((double)ens->ntraj);
  Etot = Ekin + Epot;

}





}// namespace libdyn
}// liblibra
