/*********************************************************************************
* Copyright (C) 2012 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

namespace libdyn{



void compute_hopping_probabilities_fssh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T){

// Assume ham internal processings are done outside
// assume matrix g is allocated


  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;


  // Verify if memory is allocated, if not generate
  if(g->num_of_elems!=el->nstates*el->nstates){cout<<"Matrix g is not allocated\n"; exit(0);  }


  // Now calculate the hopping probabilities
  for(i=0;i<el->nstates;i++){
    sum = 0.0;
    double a_ii = el->q[i]*el->q[i] + el->p[i]*el->p[i]; // c_i^* * c_i

    for(j=0;j<el->nstates;j++){

      if(i!=j){ // according to formula the diagonal probability P(i->i) should be zero

        // Use very general expression:
        // Note: the sign here is not very obvious! Keep in mind:
        // Re(i*z) = -Im(z)  but  Im(i*z) = Re(z)
        // My formula is: P(i->j) = (2*dt/(hbar*|c_i|^2)) * ( -Im(H_ij * c_i^* * c_j) )
        double cij_re = el->q[i]*el->q[j] + el->p[i]*el->p[j];
        double cij_im = el->q[i]*el->p[j] - el->p[i]*el->q[j];
        double imHaij = cij_re * ham->Hvib(i,j).imag() + cij_im * ham->Hvib(i,j).real(); // Im(H_ij * c_i^* * c_j)


        if(a_ii<1e-8){ g_ij = 0.0; }  // avoid division by zero
        else{

          g_ij = -2.0*dt*imHaij/a_ii;  // This is a general case


          if(use_boltz_factor){

            if(ham->H(j,j).real()>ham->H(i,i).real()){
              argg = -(ham->H(j,j).real() - ham->H(i,i).real())/(kb*T);        
              if(argg<500.0){ g_ij = g_ij * std::exp(argg); }
            }

          }// if use_boltz_factor

          if(g_ij<0.0){  g_ij = 0.0; }

        }// else

        g->set(i,j,g_ij);
        sum = sum + g_ij;
      }
      else{ g->set(i,j,0.0); }


    }// for j

    g->set(i,i,1.0 - sum);

  }// for i

}// compute probabilities_fssh

void compute_hopping_probabilities_fssh(Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
                                        double dt, int use_boltz_factor,double T){

  compute_hopping_probabilities_fssh(&mol, &el, &ham, &g, dt, use_boltz_factor, T);

}

void compute_hopping_probabilities_fssh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T){

  compute_hopping_probabilities_fssh(&ens.mol[i], &ens.el[i], ens.ham[i], &g, dt, use_boltz_factor, T);

}





void compute_hopping_probabilities_gfsh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T){


  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;

  // Verify if memory is allocated, if not generate
  if(g->num_of_elems!=el->nstates*el->nstates){cout<<"Matrix g is not allocated\n"; exit(0);  }


  // compute a_kk and a_dot_kk
  vector<double> a(el->nstates,0.0);
  vector<double> a_dot(el->nstates,0.0);
  double norm = 0.0; // normalization factor

  for(i=0;i<el->nstates;i++){
    a[i] = (el->q[i]*el->q[i] + el->p[i]*el->p[i]); // c_j^* * c_j 

    double q_dot = 0.0;
    double p_dot = 0.0;

    for(j=0;j<el->nstates;j++){

      q_dot += ( ham->Hvib(i,j).real()*el->p[j] + ham->Hvib(i,j).imag()*el->q[j]);
      p_dot -= ( ham->Hvib(i,j).real()*el->q[j] - ham->Hvib(i,j).imag()*el->p[j]);

    }// for j

    a_dot[i] = 2.0*(el->q[i]* q_dot + el->p[i] * p_dot);


    if(a_dot[i]<0.0){ norm += a_dot[i]; } // total rate of population decrease in all decaying states

  }// for i


  // Now calculate the hopping probabilities
  for(i=0;i<el->nstates;i++){       
    double sumg = 0.0;

    for(j=0;j<el->nstates;j++){

      if(j!=i){  // off-diagonal = probabilities to hop to other states

        if(a[i]<1e-12){  g->set(i,j,0.0); }  // since the initial population is almost zero, so no need for hops
        else{

          g->set(i,j,  dt*(a_dot[j]/a[i]) * a_dot[i] / norm);  
 
          if(g->get(i,j)<0.0){  // since norm is negative, than this condition means that a_dot[i] and a_dot[j] have same signs
                                // which is bad - so no transitions are assigned
            g->set(i,j,0.0);
          }
          else{  // here we have opposite signs of a_dot[i] and a_dot[j], but this is not enough yet
            if(a_dot[i]<0.0 & a_dot[j]>0.0){ ;; } // this is out transition probability, but it is already computed
            else{  g->set(i,j,0.0); } // wrong transition
          }

        }// a[i]>1e-12

        sumg += g->get(i,j);

      }
    }// for j

    g->set(i,i, 1.0 - sumg);  // probability to stay in state i

  }// for i



}// compute probabilities gfsh


void compute_hopping_probabilities_gfsh(Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
                                        double dt, int use_boltz_factor,double T){

  compute_hopping_probabilities_gfsh(&mol, &el, &ham, &g, dt, use_boltz_factor, T);

}

void compute_hopping_probabilities_gfsh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T){

  compute_hopping_probabilities_gfsh(&ens.mol[i], &ens.el[i], ens.ham[i], &g, dt, use_boltz_factor, T);

}





void compute_hopping_probabilities_mssh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T){


  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;


  // Verify if memory is allocated, if not generate
  if(g->num_of_elems!=el->nstates*el->nstates){cout<<"Matrix g is not allocated\n"; exit(0);  }


  // Now calculate the hopping probabilities
  for(j=0;j<el->nstates;j++){
    double gjj = (el->q[j]*el->q[j] + el->p[j]*el->p[j]); // c_j^* * c_j 

    for(i=0;i<el->nstates;i++){   g->set(i,j,gjj);   }// for i

  }// for j


}// compute probabilities mssh

void compute_hopping_probabilities_mssh(Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
                                        double dt, int use_boltz_factor,double T){

  compute_hopping_probabilities_mssh(&mol, &el, &ham, &g, dt, use_boltz_factor, T);

}

void compute_hopping_probabilities_mssh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T){

  compute_hopping_probabilities_mssh(&ens.mol[i], &ens.el[i], ens.ham[i], &g, dt, use_boltz_factor, T);

}



void hop(int& initstate, Nuclear* mol, Hamiltonian* ham, double ksi, MATRIX* g, int do_rescaling, int rep, int do_reverse){
// Do actual hop from state initstate
// initstate - state from which we try to hop out - it will also be updated after the hop has happened
// mol - nuclear DOF
// ham - handler of Hamiltonian
// ksi - a random number
// g   - hopping probabilities matrix
// do_rescaling - flag to turn on/off CPA: 0 - no rescaling (CPA), 1 - do rescaling (back-reaction)
// rep - representation:  0 - for diabatic, 1 - for adiabatic

  int nstates = g->num_of_cols;
  double left, right; left = right = 0.0;
  int finstate = initstate;

  for(int i=0;i<nstates;i++){
    if(i==0){left = 0.0; right = g->get(initstate,i); }
    else{  left = right; right = right + g->get(initstate,i); }
 
    if((left<ksi) && (ksi<=right)){  finstate = i;  }
  }


  if(finstate!=initstate){

    if(!do_rescaling){ initstate = finstate; }        // CPA-style, no velocity rescaling
    else{                                             // Possibly rescale velocities - normal inclusion of back-electron reaction

      // state is changed or preserved in the function
      if(rep==0){
        rescale_velocities_diabatic(mol,ham,finstate,initstate); 
      }
      else if(rep==1){
        rescale_velocities_adiabatic(mol,ham,finstate,initstate,do_reverse); 
      }

    }// else

  }// finstate!=initstate

  initstate = finstate;

}// hop

int hop(int initstate, Nuclear& mol, Hamiltonian& ham, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse){

  int res = initstate; 
  hop(res, &mol, &ham, ksi, &g, do_rescaling, rep, do_reverse);

  return res;

}

int hop(int initstate, Ensemble& ens, int i, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse){

  int res = initstate; 
  hop(res, &ens.mol[i], ens.ham[i], ksi, &g, do_rescaling, rep, do_reverse);

  return res;

}



void rescale_velocities_adiabatic(Nuclear* mol, Hamiltonian* ham, int& new_st,int& old_st, int do_reverse){
  // Calculate auxiliary variables to determine the case
  // Here i - old state, j - new state
  int st;

  if(new_st!=old_st){

    int final_st = old_st;  // no hopping by default

    // According to Fabiano
    double a_ij = 0.0;
    double b_ij = 0.0;

    for(int k=0;k<mol->nnucl;k++){

      double D = ham->D(old_st,new_st,k).real(); // derivative coupling w.r.t. nuclear DOF k

      a_ij += 0.5*(D*D / mol->mass[k]); 
      b_ij += (D*mol->p[k])/mol->mass[k];
    }
    double det = b_ij*b_ij + 4.0*a_ij*(ham->Hvib(old_st,old_st).real() - ham->Hvib(new_st,new_st).real());


    // Calculate the scaling factor and new state
    double gamma_ij = 0.0;

    if(det<0.0){

      if(do_reverse){     gamma_ij = b_ij / a_ij;}
      else{ gamma_ij = 0.0;  }

      final_st = old_st; // # hop does not occur - frustrated hop

    }
    else{
      if(b_ij<0){ gamma_ij = 0.5*(b_ij + sqrt(det))/a_ij; }
      else{       gamma_ij = 0.5*(b_ij - sqrt(det))/a_ij; }
      final_st = new_st;
    }

    //Rescale velocities and do the hop
    for(int k=0;k<mol->nnucl;k++){ 
      double D = ham->D(old_st,new_st,k).real(); 
      mol->p[k] = mol->p[k] - gamma_ij * D; 
    }

    st = final_st;

  }
  else{ st = old_st; }

  new_st = st;


} // rescale velocities adiabatic

int rescale_velocities_adiabatic(Nuclear& mol, Hamiltonian& ham, int old_st, int do_reverse){

  int new_st = old_st;
  rescale_velocities_adiabatic(&mol, &ham, new_st, old_st, do_reverse);
  return new_st;

}


void rescale_velocities_diabatic(Nuclear* mol, Hamiltonian* ham, int& new_st,int& old_st){

// In this case, derivative couplings are zero by definition of diabatic (position-independent) states
// so one can not rescale velocities along directions of derivative couplings, since there is no such directions
// We just scale velocities uniformly - based on energy conservation principle
// In principle, this rescaling procedure can be applied to surface hopping scheme in adiabatic basis, too

  double T_i = compute_kinetic_energy(mol); // initial kinetic energy
  double E_i = ham->Hvib(old_st,old_st).real();// initial potential energy
  double E_f = ham->Hvib(new_st,new_st).real();// final potential energy
  
  double T_f = T_i + E_i - E_f;             // pedicted final kinetic energy

  if(T_f>0.0){  // hop is possible - accept it
 
    double scl = 0.0; // rescaling factor
     
    if(T_i>0.0){   scl = std::sqrt(T_f/T_i);  }
    else{ scl = 0.0; }

    // Rescale velocities
    for(int i=0;i<mol->nnucl;i++){  mol->p[i] *= scl; }

    // state stays the same

  }
  else{  // hop is not possible - reject it

    new_st = old_st;

  }

} // rescale velocities diabatic

int rescale_velocities_diabatic(Nuclear& mol, Hamiltonian& ham, int old_st){

  int new_st = old_st;
  rescale_velocities_diabatic(&mol, &ham, new_st, old_st);
  return new_st;

}




}// namespace libdyn


