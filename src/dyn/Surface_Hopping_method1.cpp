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

#include "Surface_Hopping.h"
#include "Surface_Hopping_method1.h"
#include "Energy_and_Forces.h"

namespace libdyn{



void compute_hopping_probabilities_esh(Ensemble& ens, MATRIX* g, double dt, int use_boltz_factor,double T){

// esh - entangled surface hopping

// Assume ham internal processings are done outside
// assume matrix g is allocated


  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k,n, traj;
  double sum,g_ij,argg;


  // Verify if memory is allocated, if not generate
  if(g->num_of_elems!=ens.nelec * ens.nelec){cout<<"Matrix g is not allocated\n"; exit(0);  }


  // Now calculate the hopping probabilities
  // denominator
  vector<double> a_ii(ens.nelec, 0.0);

  for(i=0;i<ens.nelec;i++){
    for(traj=0; traj<ens.ntraj; traj++){

      a_ii[i] += ens.el[traj].q[i] * ens.el[traj].q[i] + ens.el[traj].p[i] * ens.el[traj].p[i]; // c_i^* * c_i

    }// for traj
    a_ii[i] /= ((double)ens.ntraj);

  }// for i

  // Kinetic energy spread:
  // Average values
  vector<double> p_ave(ens.nnucl,0.0);
  for(n=0; n<ens.nnucl; n++){  
    for(traj=0; traj<ens.ntraj; traj++){  

      p_ave[n] += ens.mol[traj].p[n];
    
    }// for traj
    p_ave[n] /= ((double)ens.ntraj);
  }// for n

  // Sigmas
  vector<double> p_sigma2(ens.nnucl,0.0);
  for(n=0; n<ens.nnucl; n++){  
    for(traj=0; traj<ens.ntraj; traj++){  

      p_sigma2[n] += (ens.mol[traj].p[n] - p_ave[n]) * (ens.mol[traj].p[n] - p_ave[n]);
    
    }// for traj
    p_sigma2[n] =  p_sigma2[n]/ ((double)ens.ntraj);
    
  }// for n

//  total_sigma = sqrt(0.1 * 2.0*2000.0*T*kb); // mass = 2000.0


  // nominator
  for(i=0;i<ens.nelec;i++){    
    double sum = 0.0;

    for(j=0;j<ens.nelec;j++){      
      if(i!=j){ // according to formula the diagonal probability P(i->i) should be zero
        double imHaij = 0.0;
        for(traj=0; traj<ens.ntraj; traj++){
 
          // Use very general expression:
          // Note: the sign here is not very obvious! Keep in mind:
          // Re(i*z) = -Im(z)  but  Im(i*z) = Re(z)
          // My formula is: P(i->j) = (2*dt/(hbar*|c_i|^2)) * ( -Im(H_ij * c_i^* * c_j) )
          double cij_re = ens.el[traj].q[i] * ens.el[traj].q[j] + ens.el[traj].p[i] * ens.el[traj].p[j];
          double cij_im = ens.el[traj].q[i] * ens.el[traj].p[j] - ens.el[traj].p[i] * ens.el[traj].q[j];
          imHaij += cij_re * ens.ham[traj]->Hvib(i,j).imag() + cij_im * ens.ham[traj]->Hvib(i,j).real(); // Im(H_ij * c_i^* * c_j)     

          // Coherence term:
          double coh = 0.0;
          for(n=0; n<ens.nnucl; n++){  
            coh += (2.0*ens.mol[traj].p[n] * ens.mol[traj].f[n]/p_sigma2[n]);
          }
          coh /= ( (ens.nelec-1) ); // equally distribute coherence over all states
          coh *= (ens.el[traj].q[i] * ens.el[traj].q[i] + ens.el[traj].p[i] * ens.el[traj].p[i]); 

          imHaij -= coh;

        }// for traj

        imHaij /= ((double)ens.ntraj);


        if(a_ii[i]<1e-8){ g_ij = 0.0; }  // avoid division by zero
        else{

          g_ij += -2.0*dt*imHaij/a_ii[i];  // This is a general case
          if(g_ij<0.0){  g_ij = 0.0; }

        }// else

        g->set(i,j,g_ij);
        sum = sum + g_ij;

      }// i!=j

      else{ g->set(i,j,0.0); }

    }// for j

    g->set(i,i,1.0 - sum);

  }// for i

}// compute probabilities_esh




void hop(int ntraj, vector<int>& initstate, vector<Nuclear*>& mol, vector<Hamiltonian*>& ham, 
         vector<double> ksi, vector<MATRIX*>& g, int do_rescaling, int rep, int do_reverse){
// Do actual hop from state initstate
// initstate - state from which we try to hop out - it will also be updated after the hop has happened
// mol - nuclear DOF
// ham - handler of Hamiltonian
// ksi - a random number
// g   - hopping probabilities matrix
// do_rescaling - flag to turn on/off CPA: 0 - no rescaling (CPA), 1 - do rescaling (back-reaction)
// rep - representation:  0 - for diabatic, 1 - for adiabatic

  int nstates = g[0]->num_of_cols;
  double left, right; 
  vector<int> finstate; finstate = initstate;

  for(int traj=0;traj<ntraj++;traj++){

    left = right = 0.0;

    for(int i=0;i<nstates;i++){
      if(i==0){left = 0.0; right = g[traj]->get(initstate[traj],i); }
      else{  left = right; right = right + g[traj]->get(initstate[traj],i); }
 
      if((left<ksi[traj]) && (ksi[traj]<=right)){  finstate[traj] = i;  }

    }// for i

  }// for all copies (trajectories) of the system

  int status = 1; // initial and final are equal
  for(traj=0; traj<ntraj; traj++){  if(finstate[traj]!=initstate[traj]){ status = 0; } }

  if(status==0){  // different multi-particle states

    if(!do_rescaling){ initstate = finstate; }        // CPA-style, no velocity rescaling
    else{                                             // Possibly rescale velocities - normal inclusion of back-electron reaction

      // state is changed or preserved in the function
      if(rep==0){
        //rescale_velocities_diabatic(mol,ham,finstate,initstate); 
      }
      else if(rep==1){
        rescale_velocities_adiabatic(ntraj,mol,ham,finstate,initstate,do_reverse); 
      }

    }// else

  }// finstate!=initstate

  initstate = finstate;

}// hop


vector<int>
hop(int ntraj, vector<int> initstate, vector<Nuclear>& mol, vector<Hamiltonian>& ham, 
    vector<double> ksi, vector<MATRIX>& g, int do_rescaling, int rep, int do_reverse){

/*
  vector<int> res; res = initstate;
  vector<Nuclear*> _mol;
  vector<Hamiltonian*> _ham;
  vector<MATRIX*> _g;

  for(int traj=0; traj<ntraj; traj++){
    _mol[traj] = &mol[traj];
    _ham[traj] = &ham[traj];
    _g[traj] = &g[traj];
  }


  hop(ntraj, res, _mol, _ham, ksi, _g, do_rescaling, rep, do_reverse);

  return res;
*/

  int nstates = g[0].num_of_cols;
  double left, right; 
  vector<int> finstate; finstate = initstate;

  for(int traj=0;traj<ntraj;traj++){

    left = right = 0.0;

    for(int i=0;i<nstates;i++){
      if(i==0){left = 0.0; right = g[traj].get(initstate[traj],i); }
      else{  left = right; right = right + g[traj].get(initstate[traj],i); }
 
      if((left<ksi[traj]) && (ksi[traj]<=right)){  finstate[traj] = i;  }

    }// for i

  }// for all copies (trajectories) of the system

  int status = 1; // initial and final are equal
  for(traj=0; traj<ntraj; traj++){  if(finstate[traj]!=initstate[traj]){ status = 0; } }

  if(status==0){  // different multi-particle states

    if(!do_rescaling){ initstate = finstate; }        // CPA-style, no velocity rescaling
    else{                                             // Possibly rescale velocities - normal inclusion of back-electron reaction

      // state is changed or preserved in the function
      if(rep==0){
        //rescale_velocities_diabatic(mol,ham,finstate,initstate); 
      }
      else if(rep==1){
        rescale_velocities_adiabatic(ntraj,mol,ham,finstate,initstate,do_reverse); 
      }

    }// else

  }// finstate!=initstate


  return finstate;  

}


boost::python::list
hop(int ntraj, boost::python::list initstate, boost::python::list mol, boost::python::list ham, 
    boost::python::list ksi, boost::python::list g, int do_rescaling, int rep, int do_reverse){
/*
  vector<int> res; 
  vector<Nuclear*> _mol;
  vector<Hamiltonian*> _ham;
  vector<MATRIX*> _g;

  for(int traj=0; traj<ntraj; traj++){
     res[traj] = extract<int>(initstate[traj]);
//    _mol[traj] = extract<Nuclear>(mol[traj]);
//    _ham[traj] = extract<Hamiltonian>(ham[traj]);
//    _g[traj] = extract<MATRIX>(g[traj]);
  }

//  hop(ntraj, res, _mol, _ham, ksi, _g, do_rescaling, rep, do_reverse);

//  boost::python::list lres; lres = [];

//  for(int traj=0; traj<ntraj; traj++){
//     lres.append(res[traj]);
//  }

  return lres;
*/   
}




void rescale_velocities_adiabatic(int ntraj, vector<Nuclear*>& mol, vector<Hamiltonian*>& ham,
      vector<int>& new_st, vector<int>& old_st, int do_reverse){

  // Calculate auxiliary variables to determine the case
  // Here i - old state, j - new state
  vector<int> st;

  int status = 1; // new_st == old_st
  for(int traj=0; traj<ntraj; traj++){
    if(new_st[traj]!=old_st[traj]){ status=0; }
  }

  if(status==0){  // different states

    vector<int> final_st; final_st = old_st;  // no hopping by default

    // According to Fabiano
    double a_ij = 0.0;
    double b_ij = 0.0;
    double E_old = 0.0;
    double E_new = 0.0;
    
    for(int traj=0;traj<ntraj;traj++){
      for(int k=0;k<mol[traj]->nnucl;k++){

        double D = ham[traj]->D(old_st[traj],new_st[traj],k).real(); // derivative coupling w.r.t. nuclear DOF k

        a_ij += 0.5*(D*D / mol[traj]->mass[k]); 
        b_ij += (D*mol[traj]->p[k])/mol[traj]->mass[k];

      }// for k

      E_old += ham[traj]->Hvib(old_st[traj],old_st[traj]).real();
      E_new += ham[traj]->Hvib(new_st[traj],new_st[traj]).real();

    }


    double det = b_ij*b_ij + 4.0*a_ij*(E_old - E_new);


    // Calculate the scaling factor and new state
    double gamma_ij = 0.0;

    if(det<0.0){

      if(do_reverse){     gamma_ij = b_ij / a_ij;}
      else{ gamma_ij = 0.0;  }

      final_st = old_st; // # multi-particle hop does not occur - frustrated hop

    }
    else{
      if(b_ij<0){ gamma_ij = 0.5*(b_ij + sqrt(det))/a_ij; }
      else{       gamma_ij = 0.5*(b_ij - sqrt(det))/a_ij; }
      final_st = new_st;
    }

    //Rescale velocities and do the hop
    for(int traj=0; traj<ntraj; traj++){
      for(int k=0;k<mol[traj]->nnucl;k++){ 

        double D = ham[traj]->D(old_st[traj],new_st[traj],k).real(); 
        mol[traj]->p[k] = mol[traj]->p[k] - gamma_ij * D; 

      }// for k
    }// for traj

    st = final_st;

  }
  else{ st = old_st; }

  new_st = st;


} // rescale velocities adiabatic



void rescale_velocities_adiabatic(int ntraj, vector<Nuclear>& mol, vector<Hamiltonian>& ham,
      vector<int>& new_st, vector<int>& old_st, int do_reverse){

  // Calculate auxiliary variables to determine the case
  // Here i - old state, j - new state
  vector<int> st;

  int status = 1; // new_st == old_st
  for(int traj=0; traj<ntraj; traj++){
    if(new_st[traj]!=old_st[traj]){ status=0; }
  }

  if(status==0){  // different states

    vector<int> final_st; final_st = old_st;  // no hopping by default

    // According to Fabiano
    double a_ij = 0.0;
    double b_ij = 0.0;
    double E_old = 0.0;
    double E_new = 0.0;
    
    for(int traj=0;traj<ntraj;traj++){
      for(int k=0;k<mol[traj].nnucl;k++){

        double D = ham[traj].D(old_st[traj],new_st[traj],k).real(); // derivative coupling w.r.t. nuclear DOF k

        a_ij += 0.5*(D*D / mol[traj].mass[k]); 
        b_ij += (D*mol[traj].p[k])/mol[traj].mass[k];

      }// for k

      E_old += ham[traj].Hvib(old_st[traj],old_st[traj]).real();
      E_new += ham[traj].Hvib(new_st[traj],new_st[traj]).real();

    }


    double det = b_ij*b_ij + 4.0*a_ij*(E_old - E_new);


    // Calculate the scaling factor and new state
    double gamma_ij = 0.0;

    if(det<0.0){

      if(do_reverse){     gamma_ij = b_ij / a_ij;}
      else{ gamma_ij = 0.0;  }

      final_st = old_st; // # multi-particle hop does not occur - frustrated hop

    }
    else{
      if(b_ij<0){ gamma_ij = 0.5*(b_ij + sqrt(det))/a_ij; }
      else{       gamma_ij = 0.5*(b_ij - sqrt(det))/a_ij; }
      final_st = new_st;
    }

    //Rescale velocities and do the hop
    for(int traj=0; traj<ntraj; traj++){
      for(int k=0;k<mol[traj].nnucl;k++){ 

        double D = ham[traj].D(old_st[traj],new_st[traj],k).real(); 
        mol[traj].p[k] = mol[traj].p[k] - gamma_ij * D; 

      }// for k
    }// for traj

    st = final_st;

  }
  else{ st = old_st; }

  new_st = st;


} // rescale velocities adiabatic



/*
int rescale_velocities_adiabatic(Nuclear& mol, Hamiltonian& ham, int old_st, int do_reverse){

  int new_st = old_st;
  rescale_velocities_adiabatic(&mol, &ham, new_st, old_st, do_reverse);
  return new_st;

}
*/




}// namespace libdyn

