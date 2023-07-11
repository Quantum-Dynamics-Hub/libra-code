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
  \file Ensemble.cpp
  \brief The file implements Python export function
    
*/

#include "Ensemble.h"

/// liblibra namespace
namespace liblibra{


/// libdyn namespace
namespace libdyn{

/// libensemble namespace
namespace libensemble{



Ensemble::~Ensemble(){
/** 
  \brief Destructor

  Doesn't do anything

*/


}



void Ensemble::_init(int _ntraj, int _nelec, int _nnucl){
/**
  \brief Ensemble object initializer 
  
  Allocate memory for an ensemble of ntraj trajectories, with all electronic components 
  represented in a basis of nstates electronic state and with nuclear component having the
  dimensionality of nnucl

  \param[in] _ntraj The number of trajectories in the ensemble
  \param[in] _nelec The number of electronic (quantum) states for each trajectory
  \param[in] _nnucl The number of nuclear DOF for each trajectory

*/ 

  int i;
  ntraj = _ntraj;
  nelec = _nelec;
  nnucl = _nnucl;

  // Allocate electronic part
  el.resize(ntraj);
  for(i=0;i<ntraj;i++){   el[i] = Electronic(nelec, 0);   }  

  // Allocate nuclear part
  mol.resize(ntraj);
  for(i=0;i<ntraj;i++){ mol[i] = Nuclear(nnucl); }

  // Allocate Hamiltonian handlers
  ham.resize(ntraj);

  // Activate all trajectories
  is_active = vector<int>(ntraj,1);


/*
  // Allocate memory for statistical data  
  ave_q = vector<double>(nnucl, 0.0);
  ave_p = vector<double>(nnucl, 0.0);
  sigma_q = vector<double>(nnucl, 0.0);
  sigma_p = vector<double>(nnucl, 0.0);
*/

}


void Ensemble::ham_set_ham(int i, Hamiltonian& _ham){ 
/**
  \brief Set Hamiltonian for the given trajectory

  This function binds the external Hamiltonian (or derived) object to the one stored and used internally.
  This is the pointer assignment, so eventually, both objects (their pointers) point to the same address.
  
  \param[in] i Is the index of the trajectory for which we setup the Hamiltonian 
  \param[in] _ham The external object, which can be modified outside, but will still affect the Ensemble computations
*/
  ham[i] = &_ham;
}

void Ensemble::ham_set_ham(int i, std::string opt, int mopt){
/**
  \brief Set Hamiltonian for the given trajectory

  This function creates a new Hamiltonian internally. So far, works only for model Hamiltonians.
  
  \param[in] i Is the index of the trajectory for which we setup the Hamiltonian 
  \param[in] opt is the option to select the type of Hamiltonian. So far, can use only opt = "model" to create 
  an array of model Hamiltonians
  \param[in] mopt The option to select the type of model Hamiltonian when creating the array of such Hamiltonians.
*/

  if(opt=="model"){
    ham[i] = new Hamiltonian_Model(mopt);    
  }// model
  /*
  else if(opt=="QM"){

    ham[i] = new Hamiltonian_Atomistic(nelec, nnucl);
    ham[i].set_Hamiltonian_type("QM");

  }// atomistic
  */
}

void Ensemble::ham_set_ham(std::string opt, int mopt){ 
/**
  \brief Set Hamiltonians for all trajectories

  This function creates a new Hamiltonians internally. So far, works only for model Hamiltonians.
  
  \param[in] opt is the option to select the type of Hamiltonian. So far, can use only opt = "model" to create 
  an array of model Hamiltonians
  \param[in] mopt The option to select the type of model Hamiltonian when creating the array of such Hamiltonians.
*/

  for(int i=0;i<ntraj;i++){  ham_set_ham(i,opt,mopt); }

}


void Ensemble::ham_set_rep(int i, int _rep){
/**
  \brief Set a representation for given trajectory

  \param[in] i Is the index of the trajectory for which we setup the representation
  \param[in] _rep The representation. Possible options: 0 (for diabatic) and 1 (for adiabatic)
*/


  ham[i]->set_rep(_rep);
}

void Ensemble::ham_set_rep(int _rep){
/**
  \brief Set a representation for all trajectories in the ensemble

  \param[in] _rep The representation. Possible options: 0 (for diabatic) and 1 (for adiabatic)
*/

  for(int i=0;i<ntraj;i++){  ham[i]->set_rep(_rep); } 
}


void Ensemble::ham_set_params(int i, vector<double>& params_){ 
/**
  \brief Set the parameters of the Hamiltonian associated with a given trajectory

  So far, the model Hamiltonians are implied
  \param[in] i Is the index of the trajectory for which we setup Hamiltonian's parameters
  \param[in] params_ The vector of double-valued parameters to feed into the Hamiltonian
*/

  ham[i]->set_params(params_); 
}

void Ensemble::ham_set_params(int i, boost::python::list params_){  
/**
  \brief Set the parameters of the Hamiltonian associated with a given trajectory (Python-friendly)

  So far, the model Hamiltonians are implied
  \param[in] i Is the index of the trajectory for which we setup Hamiltonian's parameters
  \param[in] params_ The list of parameters (could be of different types) to feed into the Hamiltonian
*/

  ham[i]->set_params(params_); 
}

void Ensemble::ham_set_params(vector<double>& params_){
/**
  \brief Set the parameters of the Hamiltonian associated with all trajectories in the ensemble

  So far, the model Hamiltonians are implied
  \param[in] params_ The vector of double-valued parameters to feed into the Hamiltonian
*/

  for(int i=0;i<ntraj;i++){ ham[i]->set_params(params_); } 
}

void Ensemble::ham_set_params(boost::python::list params_){
/**
  \brief Set the parameters of the Hamiltonian associated with all trajectories in the ensemble (Python-friendly)

  So far, the model Hamiltonians are implied
  \param[in] params_ The list of parameters (could be of different types) to feed into the Hamiltonian
*/

  for(int i=0;i<ntraj;i++){ ham[i]->set_params(params_);  } 
}


void Ensemble::ham_set_q(int i, vector<double>& q){ 
/**
  \brief Update Hamiltonian coordinates (all are real-valued scalars)

  \param[in] i Is the index of the trajectory for which we setup Hamiltonian's coordinates
  \param[in] q The vector of real-valued coordinates to be used for Hamiltonian calculations.
*/

  ham[i]->set_q(q); 
}

void Ensemble::ham_set_q(int i, boost::python::list q){ 
/**
  \brief Update Hamiltonian coordinates (all are real-valued scalars - the components of the Python list) - Python-friendly

  \param[in] i Is the index of the trajectory for which we setup Hamiltonian's coordinates
  \param[in] q The list of real-valued coordinates to be used for Hamiltonian calculations.
*/

  ham[i]->set_q(q);
}


void Ensemble::ham_set_v(int i, vector<double>& v){ 
/**
  \brief Update Hamiltonian velocities (all are real-valued scalars)

  \param[in] i Is the index of the trajectory for which we setup Hamiltonian's velocities
  \param[in] v The vector of real-valued velocities to be used for Hamiltonian calculations.

  The velocities are only needed for vibronic Hamiltonian (adiabatic representation) calculations. 
  Otherwise, they are not used.
*/

  ham[i]->set_v(v);
}

void Ensemble::ham_set_v(int i, boost::python::list v){
/**
  \brief Update Hamiltonian velocities (all are real-valued scalars -the components of Python list) - Python-friendly 

  \param[in] i Is the index of the trajectory for which we setup Hamiltonian's velocities
  \param[in] v The vector of real-valued velocities to be used for Hamiltonian calculations.

  The velocities are only needed for vibronic Hamiltonian (adiabatic representation) calculations. 
  Otherwise, they are not used.
*/

  ham[i]->set_v(v); 
}

void Ensemble::ham_set_v(){
/**
  \brief Update Hamiltonian velocities for all trajectories

  The velocities are computed using the nuclear momenta (extracted from internally-stored mol object) and the masses of
  all DOFs, also internally-stored in the mol object.

  The velocities are only needed for vibronic Hamiltonian (adiabatic representation) calculations. 
  Otherwise, they are not used.
*/


  vector<double> v(nnucl,0.0); 
  for(int i=0;i<ntraj;i++){ 
    // Velocities for all DOFs for given trajectory
    for(int n=0;n<nnucl;n++){ v[n] = mol[i].p[n]/mol[i].mass[n];  }
    ham[i]->set_v(v);
  }// for i

}

void Ensemble::ham_compute(int i){ 
/**
  \brief Perform actual Hamiltonian computations for the given trajectory

  The computations of either diabatiatic or adiabatic or both Hamiltonians are invoked, depending
  on the representation set up for this Hamiltonian and on the state of the computations of such 
  Hamiltonians (so, if no change of position/velocity has been made since the last computation of 
  given Hamiltonian, no actually computations will be carryied out, not to do usefull work). Also,
  computations of adiabatic Hamiltonians (if adiabatic representation is set up) may call computation
  of the diabatic Hamiltonians, since they may be required. On the contrary, if the diabatic Hamiltonian
  is selected, the adiabatic Hamiltonian is not updated.

  Note, that just updating momenta and positions will not lead to automatic recomputation of the 
  Hamiltonians and derivatives

*/

  ham[i]->compute();
}

void Ensemble::ham_compute_diabatic(int i){
/**
  \brief Perform actual Hamiltonian computations (only diabatic) for the given trajectory

  The computations of either diabatiatic Hamiltonian is invoked, depending
  on the representation set up for this Hamiltonian and on the state of the computations of such 
  Hamiltonian (so, if no change of position/velocity has been made since the last computation of 
  given Hamiltonian, no actually computations will be carryied out, not to do usefull work). 

  Note, that just updating momenta and positions will not lead to automatic recomputation of the 
  Hamiltonians and derivatives

*/

  ham[i]->compute_diabatic(); 
}

void Ensemble::ham_compute_adiabatic(int i){
/**
  \brief Perform actual Hamiltonian computations (only adiabatic) for the given trajectory

  The computations of  adiabatiatic Hamiltonian is invoked, depending
  on the representation set up for this Hamiltonian and on the state of the computations of such 
  Hamiltonian (so, if no change of position/velocity has been made since the last computation of 
  given Hamiltonian, no actually computations will be carryied out, not to do usefull work). 

  Note, that just updating momenta and positions will not lead to automatic recomputation of the 
  Hamiltonians and derivatives

*/

  ham[i]->compute_adiabatic();
}


std::complex<double> Ensemble::ham_H(int traj, int i,int j){ 
/**
  \brief Return electronic Hamiltonian for given trajectory

  The returned Hamiltonian depends on the selected representation - can be either diabatic or adiabatic.
  This function does not invoke actual computation - it only returns whatever exists in the internal variables.

  \param[in] traj Index of the trajectory
  \param[in] i index of electronic state
  \param[in] j index of electronic state

*/

  return ham[traj]->H(i,j); 
}

std::complex<double> Ensemble::ham_dHdq(int traj, int i,int j,int n){ 
/**
  \brief Return the derivative of electronic Hamiltonian w.r.t. nuclear DOF for given trajectory

  The returned Hamiltonian depends on the selected representation - can be either diabatic or adiabatic.
  This function does not invoke actual computation - it only returns whatever exists in the internal variables.

  \param[in] traj Index of the trajectory
  \param[in] i index of electronic state
  \param[in] j index of electronic state
  \param[in] n index of nuclear DOF w.r.t. which the differentiation is performed

*/

  return ham[traj]->dHdq(i,j,n); 
}

std::complex<double> Ensemble::ham_D(int traj, int i,int j,int n){ 
/**
  \brief Return the derivative coupling w.r.t. nuclear DOF for given trajectory

  The returned coupling depends on the selected representation - can be either diabatic or adiabatic.
  This function does not invoke actual computation - it only returns whatever exists in the internal variables.

  D = <i|d/dR_n|j> 

  \param[in] traj Index of the trajectory
  \param[in] i index of electronic state
  \param[in] j index of electronic state
  \param[in] n index of nuclear DOF w.r.t. which the coupling is computed

*/

  return ham[traj]->D(i,j,n); 
}

std::complex<double> Ensemble::ham_nac(int traj,int i,int j){ 
/**
  \brief Return the nonadiabatic coupling for given trajectory

  The returned coupling depends on the selected representation - can be either diabatic or adiabatic.
  This function does not invoke actual computation - it only returns whatever exists in the internal variables.

  nac = sum_n { dR_n/dt * <i|d/dR_n|j> }

  \param[in] traj Index of the trajectory
  \param[in] i index of electronic state
  \param[in] j index of electronic state

*/

  return ham[traj]->nac(i,j); 
}

std::complex<double> Ensemble::ham_Hvib(int traj, int i,int j){ 
/**
  \brief Return the vibronic Hamiltonian for given trajectory

  The returned Hamiltonian depends on the selected representation - can be either diabatic or adiabatic.
  This function does not invoke actual computation - it only returns whatever exists in the internal variables.

  \param[in] traj Index of the trajectory
  \param[in] i index of electronic state
  \param[in] j index of electronic state

*/

  return ham[traj]->Hvib(i,j); 
}





void Ensemble::el_propagate_electronic(int i,double dt){
/**
  \brief Propagate electronic DOF for a given trajectory

  \param[in] i Index of the trajectory
  \param[in] dt Integration time step (integration duration)

*/
             
  el[i].propagate_electronic(dt, ham[i]);

}

void Ensemble::el_propagate_electronic(double dt){
/**
  \brief Propagate electronic DOF for all trajectories

  \param[in] dt Integration time step (integration duration)

*/


  for(int i=0;i<ntraj;i++){     el[i].propagate_electronic(dt, ham[i]);  }

}

void Ensemble::mol_propagate_q(int i,double dt){
/**
  \brief Propagate nuclear coordinates for given trajectory

  \param[in] i Index of the trajectory
  \param[in] dt Integration time step (integration duration)

*/


  mol[i].propagate_q(dt);

}
void Ensemble::mol_propagate_q(double dt){
/**
  \brief Propagate nuclear coordinates for all trajectories

  \param[in] dt Integration time step (integration duration)

*/


  for(int i=0;i<ntraj;i++){    mol[i].propagate_q(dt);   }

}

void Ensemble::mol_propagate_p(int i,double dt){
/**
  \brief Propagate nuclear momenta for given trajectory

  \param[in] i Index of the trajectory
  \param[in] dt Integration time step (integration duration)

*/

  mol[i].propagate_p(dt);

}
void Ensemble::mol_propagate_p(double dt){
/**
  \brief Propagate nuclear momenta for all trajectories

  \param[in] dt Integration time step (integration duration)

*/

  for(int i=0;i<ntraj;i++){    mol[i].propagate_p(dt);   }

}




Ensemble::Ensemble(int _ntraj, int _nstates, int _nnucl){
/**
  \brief Constructor with parameters

  \param[in] _ntraj The number of the trajectories of the ensemble 
  \param[in] _nstates The number of electronic DOF (states) for each trajectoriy
  \param[in] _nnucl The number of nuclear DOF for each trajectory

*/

// Constructor
   _init(_ntraj,_nstates,_nnucl);
}


void Ensemble::se_pop(vector<double>& pops,double xmin, double xmax){
/**
  \brief Compute Schrodinger equation (coherent) populations

  The computed population take into consideration the spatial localization of trajectories -
  only those trajectories that are in the given window between xmin and xmax.
  The trajectory is considered to be in that window if ALL nuclear DOF are larger than xmin
  but are smaller than xmax.

  pops[i] - will contain the fraction of the wavefunction (as given by the average weight)
  of i-th electronic basis state that is enclosed by a box [xmin, xmax] in all dimensions (x,y,z)
  and for all nuclear degrees of freedom

  \param[out] pops The computed populations are saved into this vector
  \param[in] xmin The minimal boundary of spatial window
  \param[in] xmax The maximal boundary of spatial window

*/

  int i,traj,n;

  // In case array is of wrong size
  if(pops.size()!=nelec){
    pops.reserve(nelec);
    pops.resize(nelec,0.0);  
  }


  for(i=0;i<nelec;i++){  // all electronic states
    pops[i] = 0.0;

    for(traj=0;traj<ntraj;traj++){

      if(el[traj].istate == i){
        
        int res = 1;
        for(n=0;n<nnucl;n++){
          if(mol[traj].q[n]<xmin || mol[traj].q[n] > xmax ){  res = 0;  }
        }// for n

        if(res==1){ 
          double q = el[traj].q[i];
          double p = el[traj].p[i];
          pops[i] += (q*q + p*p);       
        }


      }// if right state
    }// for j - all trajectories

    pops[i] /= (double)ntraj;   // normalize

  }// for i - all electronic states


   

}

void Ensemble::se_pop(vector<double>& pops){
/**
  \brief Compute Schrodinger equation (coherent) populations

  Wavefunction population of all states, without regard to nuclear wavefunction localization
  This is done by taking very large box
  If you are out of this box - this is kinda strange (you need to reduce time of simulation, maybe)

  \param[out] pops The computed populations are saved into this vector

*/

  se_pop(pops,-1000000.0,1000000.0);

}

boost::python::list Ensemble::se_pop(){
/**
  \brief Compute Schrodinger equation (coherent) populations - Python-friendly version

  Wavefunction population of all states, without regard to nuclear wavefunction localization
  This is done by taking very large box
  If you are out of this box - this is kinda strange (you need to reduce time of simulation, maybe)

  The computed vector of populations will be returned as the list of the floating-point values
*/


  vector<double> pops;
  se_pop(pops);

  boost::python::list res;
  for(int i=0;i<pops.size();i++){
    res.append(pops[i]);
  }

  return res;
}

boost::python::list Ensemble::se_pop(double xmin, double xmax){
/**
  \brief Compute Schrodinger equation (coherent) populations - Python-friendly version

  The computed population take into consideration the spatial localization of trajectories -
  only those trajectories that are in the given window between xmin and xmax.
  The trajectory is considered to be in that window if ALL nuclear DOF are larger than xmin
  but are smaller than xmax.

  pops[i] - will contain the fraction of the wavefunction (as given by the average weight)
  of i-th electronic basis state that is enclosed by a box [xmin, xmax] in all dimensions (x,y,z)
  and for all nuclear degrees of freedom

  \param[in] xmin The minimal boundary of spatial window
  \param[in] xmax The maximal boundary of spatial window

  The computed vector of populations will be returned as the list of the floating-point values
*/

  vector<double> pops;
  se_pop(pops, xmin, xmax);

  boost::python::list res;
  for(int i=0;i<pops.size();i++){
    res.append(pops[i]);
  }

  return res;
}




void Ensemble::sh_pop(vector<double>& pops,double xmin, double xmax){
/**
  \brief Compute surface hopping (incoherent) populations

  The computed population take into consideration the spatial localization of trajectories -
  only those trajectories that are in the given window between xmin and xmax.
  The trajectory is considered to be in that window if ALL nuclear DOF are larger than xmin
  but are smaller than xmax.

  pops[i] - will contain the fraction of the wavefunction (as given by the average weight)
  of i-th electronic basis state that is enclosed by a box [xmin, xmax] in all dimensions (x,y,z)
  and for all nuclear degrees of freedom

  \param[out] pops The computed populations are saved into this vector
  \param[in] xmin The minimal boundary of spatial window
  \param[in] xmax The maximal boundary of spatial window

*/


  int i,traj,n;

  // In case array is of wrong size
  if(pops.size()!=nelec){
    pops.reserve(nelec);
    pops.resize(nelec,0.0);  
  }
  
  for(i=0;i<nelec;i++){  // all electronic states
    pops[i] = 0.0;

    for(traj=0;traj<ntraj;traj++){

      if(el[traj].istate == i){
        
        int res = 1;
        for(n=0;n<nnucl;n++){
          if(mol[traj].q[n]<xmin || mol[traj].q[n] > xmax ){  res = 0;  }
        }// for n

        if(res==1){    pops[i] += 1.0;        }


      }// if right state
    }// for j - all trajectories

    pops[i] /= (double)ntraj;   // normalize

  }// for i - all electronic states
  

}

void Ensemble::sh_pop(vector<double>& pops){
/**
  \brief Compute surface hopping (incoherent) populations

  Wavefunction population of all states, without regard to nuclear wavefunction localization
  This is done by taking very large box
  If you are out of this box - this is kinda strange (you need to reduce time of simulation, maybe)

  \param[out] pops The computed populations are saved into this vector

*/

  sh_pop(pops,-1000000.0,1000000.0);

}

boost::python::list Ensemble::sh_pop(){
/**
  \brief Compute surface hopping (incoherent) populations - Python-friendly version

  Wavefunction population of all states, without regard to nuclear wavefunction localization
  This is done by taking very large box
  If you are out of this box - this is kinda strange (you need to reduce time of simulation, maybe)

  The computed vector of populations will be returned as the list of the floating-point values
*/


  vector<double> pops;
  sh_pop(pops);

  boost::python::list res;
  for(int i=0;i<pops.size();i++){
    res.append(pops[i]);
  }

  return res;
}

boost::python::list Ensemble::sh_pop(double xmin, double xmax){
/**
  \brief Compute surface hopping (incoherent) populations - Python-friendly version

  The computed population take into consideration the spatial localization of trajectories -
  only those trajectories that are in the given window between xmin and xmax.
  The trajectory is considered to be in that window if ALL nuclear DOF are larger than xmin
  but are smaller than xmax.

  pops[i] - will contain the fraction of the wavefunction (as given by the average weight)
  of i-th electronic basis state that is enclosed by a box [xmin, xmax] in all dimensions (x,y,z)
  and for all nuclear degrees of freedom

  \param[in] xmin The minimal boundary of spatial window
  \param[in] xmax The maximal boundary of spatial window

  The computed vector of populations will be returned as the list of the floating-point values
*/

  vector<double> pops;
  sh_pop(pops, xmin, xmax);

  boost::python::list res;
  for(int i=0;i<pops.size();i++){
    res.append(pops[i]);
  }

  return res;
}



void Ensemble::sh_pop1(vector<double>& pops,double xmin, double xmax){
/**
  \brief Compute surface hopping (incoherent) populations

  This is specially for Marcus spin-boson problem

  The computed population take into consideration the spatial localization of trajectories -
  only those trajectories that are in the given window between xmin and xmax.
  The trajectory is considered to be in that window if ALL nuclear DOF are larger than xmin
  but are smaller than xmax.

  pops[i] - will contain the fraction of the wavefunction (as given by the average weight)
  of i-th electronic basis state that is enclosed by a box [xmin, xmax] in all dimensions (x,y,z)
  and for all nuclear degrees of freedom

  \param[out] pops The computed populations are saved into this vector
  \param[in] xmin The minimal boundary of spatial window
  \param[in] xmax The maximal boundary of spatial window

*/

  int i,j;
  vector<double> v_(ntraj);

  // In case array is of wrong size
  if(pops.size()!=nelec){
    pops.reserve(nelec);
    pops.resize(nelec,0.0);  
  }

  for(i=0;i<nelec;i++){  // all electronic states
    pops[i] = 0.0;
  }

 
  for(j=0;j<ntraj;j++){ // for all trajectories

//    ham[j]->set_status(0);
    ham[j]->set_q(mol[j].q); 
    for(i=0;i<mol[j].nnucl;i++){ v_[i] = mol[j].p[i]/mol[j].mass[i]; }  ham[j]->set_v(v_);
    ham[j]->compute();

    double E0, E1, H0, H1, V;
    ham[j]->set_rep(1);
    E0 = ham[j]->H(0,0).real();   // adiabatic
    E1 = ham[j]->H(1,1).real();

    ham[j]->set_rep(0);
    H0 = ham[j]->H(0,0).real(); // diabatic
    H1 = ham[j]->H(1,1).real();
    V  = ham[j]->H(0,1).real();

    ham[j]->set_rep(1); // return back to adiabatic

    for(i=0;i<nelec;i++){  // all electronic states

      if(mol[j].q[0] > xmin && mol[j].q[0] < xmax){ // only those trajectories that are in given range

        if(el[j].istate == 0){ // we are in 0-th adiabatic state
        
          if(i==0){             // Probability to be on 0-th (left) diabat - f0
           
            pops[0] += V*V/((H0 - E0)*(H0 - E0) + V*V);  // f0^2

          }
          else if(i==1){             // Probability to be on 1-th (right) diabat - g0
           
            pops[1] += (H0 - E0)*(H0 - E0)/((H0 - E0)*(H0 - E0) + V*V);  // g0^2

          }


        }// 0-th adiabatic state
        else if(el[j].istate == 1){

          if(i==0){             // Probability to be on 0-th (left) diabat - f1
         
            pops[0] += V*V/((H0 - E1)*(H0 - E1) + V*V);  // f1

          }
          else if(i==1){             // Probability to be on 1-th (right) diabat - g1

            pops[1] += (H0 - E1)*(H0 - E1)/((H0 - E1)*(H0 - E1) + V*V);  // g1

          }


        }// 1-th adiabatic state

      }// if min < x < max

    }// for i - all electronic states
  }// for j - all trajectories



  for(i=0;i<nelec;i++){  // all electronic states

    pops[i] /= ((double)ntraj);   // normalize

  }
  

}


void Ensemble::sh_pop1(vector<double>& pops){
/**
  \brief Compute surface hopping (incoherent) populations

  This is specially for Marcus spin-boson problem

  Wavefunction population of all states, without regard to nuclear wavefunction localization
  This is done by taking very large box
  If you are out of this box - this is kinda strange (you need to reduce time of simulation, maybe)

  \param[out] pops The computed populations are saved into this vector

*/

  sh_pop1(pops,-100000000.0,100000000.0);

}

boost::python::list Ensemble::sh_pop1(){
/**
  \brief Compute surface hopping (incoherent) populations - Python-friendly version

  This is specially for Marcus spin-boson problem

  Wavefunction population of all states, without regard to nuclear wavefunction localization
  This is done by taking very large box
  If you are out of this box - this is kinda strange (you need to reduce time of simulation, maybe)

  The computed vector of populations will be returned as the list of the floating-point values
*/

  vector<double> pops;
  sh_pop1(pops);

  boost::python::list res;
  for(int i=0;i<pops.size();i++){
    res.append(pops[i]);
  }

  return res;
}

boost::python::list Ensemble::sh_pop1(double xmin, double xmax){
/**
  \brief Compute surface hopping (incoherent) populations - Python-friendly version

  This is specially for Marcus spin-boson problem

  The computed population take into consideration the spatial localization of trajectories -
  only those trajectories that are in the given window between xmin and xmax.
  The trajectory is considered to be in that window if ALL nuclear DOF are larger than xmin
  but are smaller than xmax.

  pops[i] - will contain the fraction of the wavefunction (as given by the average weight)
  of i-th electronic basis state that is enclosed by a box [xmin, xmax] in all dimensions (x,y,z)
  and for all nuclear degrees of freedom

  \param[in] xmin The minimal boundary of spatial window
  \param[in] xmax The maximal boundary of spatial window

  The computed vector of populations will be returned as the list of the floating-point values
*/

  vector<double> pops;
  sh_pop1(pops, xmin, xmax);

  boost::python::list res;
  for(int i=0;i<pops.size();i++){
    res.append(pops[i]);
  }

  return res;
}



void Ensemble::print_map(std::string prefix, double Xmin, double Xmax, double dx, double Ymin, double Ymax, double dy, int snap){
/**
  \brief for 2D projections on XY plane. Inactive now.

*/

/*
 // for 2D projections on XY plane

  std::string snaps,filename;
  stringstream ss(stringstream::in | stringstream::out);
  stringstream ss1(stringstream::in | stringstream::out);
  stringstream ss3(stringstream::in | stringstream::out);

  ss << snap;  ss >> snaps;
  
  int Nx = floor((Xmax-Xmin)/dx);
  int Ny = floor((Ymax-Ymin)/dy);

  filename = prefix+".frame"+snaps;

  ofstream out(filename.c_str(),ios::out);

  for(int nx=0;nx<Nx;nx++){
    double X = Xmin + nx*dx;

    for(int ny=0;ny<Ny;ny++){
      double Y = Ymin + ny*dy;

      // Compute density of trajectories: number of trajectories in the box [X,X+dx] x [Y,Y+dy]
      double dens = 0.0;
      double cnt = 0.0;
      for(int traj=0;traj<size;traj++){
        for(int n=0;n<mol[traj]->Nnucl;n++){

          if( (X<mol[traj]->R[n].x) && (mol[traj]->R[n].x<=(X+dx))  &&
              (Y<mol[traj]->R[n].y) && (mol[traj]->R[n].y<=(Y+dy))
            ){
             dens += 1.0;
           }
          cnt += 1.0;
        }// for n
      }// for traj
      dens /= cnt;

      out<<X<<"  "<<Y<<"  "<<dens<<endl;
    }// for ny
    out<<"\n";
  }// for nx

  out.close();

*/
}// print map


void Ensemble::integral_flux(vector< vector<double> >& Int_flx, double Xmin, double Xmax, double dx, 
                                                                double Ymin, double Ymax, double dy){
/**
  \brief for 2D projections on XY plane. Inactive now.

*/


 // for 2D projections on XY plane
/*
  int Nx = Int_flx.size();
  int Ny = Int_flx[0].size();

  //
  double denom = 0.0;
  for(int traj=0;traj<size;traj++){
    for(int n=0;n<mol[traj]->Nnucl;n++){
      denom += 1.0;
    }
  }
  denom = 1.0/denom;



  for(int traj=0;traj<size;traj++){
    for(int n=0;n<mol[traj]->Nnucl;n++){

      int nx = floor((mol[traj]->R[n].x - Xmin)/dx);
      int ny = floor((mol[traj]->R[n].y - Ymin)/dy);

      if(nx>=0 && nx<Nx && ny>=0 && ny<Ny){
        Int_flx[nx][ny] += denom;
      }

    }
  }

*/

}// print map





void Ensemble::compute_averages(){
/**
  \brief Computing the ensemble-averaged value. Inactive now.

*/


/*
  ave_x = 0.0; 
  ave_x2 = 0.0; 
  ave_y = 0.0; 
  ave_y2 = 0.0; 

  ave_px = 0.0; 
  ave_px2 = 0.0; 
  ave_py = 0.0; 
  ave_py2 = 0.0; 

  ave_xpx = 0.0; 
  ave_ypy = 0.0; 


  for(int traj=0;traj<size;traj++){

    ave_x  += mol[traj]->R[0].x;
    ave_x2 += (mol[traj]->R[0].x * mol[traj]->R[0].x);
    ave_y  += mol[traj]->R[0].y;
    ave_y2 += (mol[traj]->R[0].y * mol[traj]->R[0].y);

    ave_px  += mol[traj]->P[0].x;
    ave_px2 += (mol[traj]->P[0].x * mol[traj]->P[0].x);
    ave_py  += mol[traj]->P[0].y;
    ave_py2 += (mol[traj]->P[0].y * mol[traj]->P[0].y);

    ave_xpx = mol[traj]->R[0].x * mol[traj]->P[0].x; 
    ave_ypy = mol[traj]->R[0].y * mol[traj]->P[0].y; 

  }

  ave_x  /= (float(size));
  ave_x2 /= (float(size));
  ave_y  /= (float(size));
  ave_y2 /= (float(size));

  ave_px  /= (float(size));
  ave_px2 /= (float(size));
  ave_py  /= (float(size));
  ave_py2 /= (float(size));

  ave_xpx /= (float(size));
  ave_ypy /= (float(size));


  sx = sqrt(ave_x2 - ave_x*ave_x); 
  sy = sqrt(ave_y2 - ave_y*ave_y); 

  psx = (ave_xpx - ave_x * ave_px)/sx;
  psy = (ave_ypy - ave_y * ave_py)/sy;
*/

}


/*
double Ensemble::Epot(){
  
  double res = 0.0; 
  for(int traj=0;traj<size;traj++){
    double ep = ham[traj]->H(mol[traj],0,0).real();
    res += ep;
  } 
  res /= ((double)size);

  return res;

}

double Ensemble::Ekin(){

  double res = 0.0; 
  for(int traj=0;traj<size;traj++){
    double ek = compute_kinetic_energy(mol[traj]);
    res += ek;
  }
  res /= ((double)size);

  return res;
}

double Ensemble::Etot(){

  return (Ekin() + Epot());

}
*/


}// namespace libensemble
}// namespace libdyn
}// liblibra
