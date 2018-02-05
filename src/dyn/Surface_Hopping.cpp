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
  \file Surface_Hopping.cpp
  \brief The file implements the functions used in surface hopping methods
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{



MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, CMATRIX& Hvib, double dt){
/**

*/

  int nstates = Coeff.n_elts;
  MATRIX g(nstates,nstates);


  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;

  CMATRIX* denmat; denmat = new CMATRIX(nstates, nstates); 
  
  *denmat = (Coeff * Coeff.H() ).conj();

  // Now calculate the hopping probabilities
  for(i=0;i<nstates;i++){
    sum = 0.0;
    double a_ii = denmat->get(i,i).real(); // c_i^* * c_i

    for(j=0;j<nstates;j++){

      if(i!=j){ // according to formula the diagonal probability P(i->i) should be zero
        // Use very general expression:
        // Note: the sign here is not very obvious! Keep in mind:
        // Re(i*z) = -Im(z)  but  Im(i*z) = Re(z)
        // My formula is: P(i->j) = (2*dt/(hbar*|c_i|^2)) * ( -Im(H_ij * c_i^* * c_j) )
        double imHaij = ( Hvib.get(i,j) * denmat->get(i,j) ).imag(); // Im(H_ij * c_i^* * c_j)


        if(a_ii<1e-8){ g_ij = 0.0; }  // avoid division by zero
        else{
//          g_ij = -2.0*dt*imHaij/a_ii;  // This is a general case - wrong sign (opposite of the one in JCC)
            g_ij = 2.0*dt*imHaij/a_ii;  // This is a general case -
          if(g_ij<0.0){  g_ij = 0.0; }

        }// else

        g.set(i,j,g_ij);
        sum = sum + g_ij;
      }
      else{ g.set(i,j,0.0); }

    }// for j

    g.set(i,i,1.0 - sum);

  }// for i

  delete denmat;

  return g;

}// fssh




void compute_hopping_probabilities_fssh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T){
/**
  \brief Compute the FSSH surface hopping probabilities for the trajectory described by mol, el, and ham
  \param[in] mol Describes the nuclear DOF
  \param[in] el Describes electronic DOF
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[out] g The matrix of hopping probabilities
  \param[in] dt Time duration of nuclear propagation step
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Assume ham internal processings are done outside
  Assume matrix g is allocated

  Abbreviation: FSSH - fewest switches surface hopping
  References: 
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061–1071. - the original paper
  (2) Fabiano, E.; Keal, T. W.; Thiel, W. Implementation of Surface Hopping Molecular Dynamics Using Semiempirical Methods. Chem. Phys. 2008, 349, 334–347.
  Here, we generalized the formula, so it works equally well for both diabatic and adiabatic representations

*/



  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;


  // Verify if memory is allocated, if not generate
  if(g->n_elts!=el->nstates*el->nstates){cout<<"Matrix g is not allocated\n"; exit(0);  }


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

//        double cij_re = el->q[i]*el->q[j] + el->p[i]*el->p[j];
//        double cij_im = el->q[i]*el->p[j] - el->p[i]*el->q[j];
//        double imHaij = cij_re * ham->Hvib(i,j).imag() + cij_im * ham->Hvib(i,j).real(); // Im(H_ij * c_i^* * c_j)
        double imHaij = (ham->Hvib(i,j) * el->rho(i,j) ).imag();


        if(a_ii<1e-8){ g_ij = 0.0; }  // avoid division by zero
        else{

//          g_ij = -2.0*dt*imHaij/a_ii;  // This is a general case - wrong sign!
          g_ij = 2.0*dt*imHaij/a_ii;  // This is a general case


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
/**
  \brief Compute the FSSH surface hopping probabilities for the trajectory described by mol, el, and ham - Python-friendly
  \param[in] mol Describes the nuclear DOF
  \param[in] el Describes electronic DOF
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[out] g The matrix of hopping probabilities
  \param[in] dt Time duration of nuclear propagation step
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Assume ham internal processings are done outside
  Assume matrix g is allocated

  Abbreviation: FSSH - fewest switches surface hopping
  References: 
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061–1071. - the original paper
  (2) Fabiano, E.; Keal, T. W.; Thiel, W. Implementation of Surface Hopping Molecular Dynamics Using Semiempirical Methods. Chem. Phys. 2008, 349, 334–347.
  Here, we generalized the formula, so it works equally well for both diabatic and adiabatic representations

*/


  compute_hopping_probabilities_fssh(&mol, &el, &ham, &g, dt, use_boltz_factor, T);

}

void compute_hopping_probabilities_fssh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T){
/**
  \brief Compute the FSSH surface hopping probabilities for a selected trajectory of an ensemble - Python friendly
  \param[in] ens Describes the ensemble of trajectories
  \param[in] i Is the index of the trajectory of interest
  \param[out] g The matrix of hopping probabilities for this tranjectory
  \param[in] dt Time duration of nuclear propagation step
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Assume ham internal processings are done outside
  Assume matrix g is allocated

  Abbreviation: FSSH - fewest switches surface hopping
  References: 
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061–1071. - the original paper
  (2) Fabiano, E.; Keal, T. W.; Thiel, W. Implementation of Surface Hopping Molecular Dynamics Using Semiempirical Methods. Chem. Phys. 2008, 349, 334–347.
  Here, we generalized the formula, so it works equally well for both diabatic and adiabatic representations

*/


  compute_hopping_probabilities_fssh(&ens.mol[i], &ens.el[i], ens.ham[i], &g, dt, use_boltz_factor, T);

}








MATRIX compute_hopping_probabilities_gfsh(CMATRIX& Coeff, CMATRIX& Hvib, double dt){

/**
  \brief Compute the GFSH surface hopping probabilities for the trajectory described by mol, el, and ham
  \param[in] Coeff Wavefunction amplitudes
  \param[in] Hvib vibronic Hamiltonian matrix
  \param[in] dt Time duration of nuclear propagation step

  Abbreviation: GFSH - global flux surface hopping
  References: 
  (1) Wang, L.; Trivedi, D.; Prezhdo, O. V. Global Flux Surface Hopping Approach for Mixed Quantum-Classical Dynamics. J. Chem. Theory Comput. 2014, 10, 3598–3605.

*/

  int nstates = Coeff.n_elts;
  MATRIX g(nstates,nstates);

  CMATRIX* denmat; denmat = new CMATRIX(nstates, nstates);   
  *denmat = (Coeff * Coeff.H() ).conj();

  CMATRIX* denmat_dot; denmat_dot = new CMATRIX(nstates, nstates);   
  *denmat_dot = ((*denmat) *  Hvib.conj() - Hvib * (*denmat)) * complex<double>(0.0, 1.0);


  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;


  // compute a_kk and a_dot_kk
  vector<double> a(nstates,0.0);
  vector<double> a_dot(nstates,0.0);
  double norm = 0.0; // normalization factor

  for(i=0;i<nstates;i++){
    a[i] = denmat->get(i,i).real();
    a_dot[i] = denmat_dot->get(i,i).real();

    if(a_dot[i]<0.0){ norm += a_dot[i]; } // total rate of population decrease in all decaying states

  }// for i


  // Now calculate the hopping probabilities
  for(i=0;i<nstates;i++){       
    double sumg = 0.0;

    for(j=0;j<nstates;j++){

      if(j!=i){  // off-diagonal = probabilities to hop to other states

        if(a[i]<1e-12){  g.set(i,j,0.0); }  // since the initial population is almost zero, so no need for hops
        else{

          g.set(i,j,  dt*(a_dot[j]/a[i]) * a_dot[i] / norm);  
 
          if(g.get(i,j)<0.0){  // since norm is negative, than this condition means that a_dot[i] and a_dot[j] have same signs
                                // which is bad - so no transitions are assigned
            g.set(i,j,0.0);
          }
          else{  // here we have opposite signs of a_dot[i] and a_dot[j], but this is not enough yet
            if(a_dot[i]<0.0 & a_dot[j]>0.0){ ;; } // this is out transition probability, but it is already computed
            else{  g.set(i,j,0.0); } // wrong transition
          }

        }// a[i]>1e-12

        sumg += g.get(i,j);

      }
    }// for j

    g.set(i,i, 1.0 - sumg);  // probability to stay in state i

  }// for i

  delete denmat;
  delete denmat_dot;

}// fssh




void compute_hopping_probabilities_gfsh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T){
/**
  \brief Compute the GFSH surface hopping probabilities for the trajectory described by mol, el, and ham
  \param[in] mol Describes the nuclear DOF
  \param[in] el Describes electronic DOF
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[out] g The matrix of hopping probabilities
  \param[in] dt Time duration of nuclear propagation step
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Assume ham internal processings are done outside
  Assume matrix g is allocated

  Abbreviation: GFSH - global flux surface hopping
  References: 
  (1) Wang, L.; Trivedi, D.; Prezhdo, O. V. Global Flux Surface Hopping Approach for Mixed Quantum-Classical Dynamics. J. Chem. Theory Comput. 2014, 10, 3598–3605.

*/


  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;

  // Verify if memory is allocated, if not generate
  if(g->n_elts!=el->nstates*el->nstates){cout<<"Matrix g is not allocated\n"; exit(0);  }


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
/**
  \brief Compute the GFSH surface hopping probabilities for the trajectory described by mol, el, and ham - Python-friendly
  \param[in] mol Describes the nuclear DOF
  \param[in] el Describes electronic DOF
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[out] g The matrix of hopping probabilities
  \param[in] dt Time duration of nuclear propagation step
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Assume ham internal processings are done outside
  Assume matrix g is allocated

  Abbreviation: GFSH - global flux surface hopping
  References: 
  (1) Wang, L.; Trivedi, D.; Prezhdo, O. V. Global Flux Surface Hopping Approach for Mixed Quantum-Classical Dynamics. J. Chem. Theory Comput. 2014, 10, 3598–3605.

*/


  compute_hopping_probabilities_gfsh(&mol, &el, &ham, &g, dt, use_boltz_factor, T);

}

void compute_hopping_probabilities_gfsh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T){
/**
  \brief Compute the GFSH surface hopping probabilities for a selected trajectory of an ensemble - Python friendly
  \param[in] ens Describes the ensemble of trajectories
  \param[in] i Is the index of the trajectory of interest
  \param[out] g The matrix of hopping probabilities for this tranjectory
  \param[in] dt Time duration of nuclear propagation step
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Assume ham internal processings are done outside
  Assume matrix g is allocated

  Abbreviation: GFSH - global flux surface hopping
  References: 
  (1) Wang, L.; Trivedi, D.; Prezhdo, O. V. Global Flux Surface Hopping Approach for Mixed Quantum-Classical Dynamics. J. Chem. Theory Comput. 2014, 10, 3598–3605.

*/


  compute_hopping_probabilities_gfsh(&ens.mol[i], &ens.el[i], ens.ham[i], &g, dt, use_boltz_factor, T);

}




MATRIX compute_hopping_probabilities_mssh(CMATRIX& Coeff){
/**
  \brief Compute the MSSH surface hopping probabilities for the trajectory described by the coefficients Coeff

  This is the version taking the minimal amount of input information

  \param[in] Coeff The amplitudes of different basis states in the coherent superposition. This matrix is assumed to be
  a column-vector, so of the size N x 1, where N is the number of basis excited states

  The function returns the matrix g with hopping probabilities

  Abbreviation: MSSH - Markov state surface hopping
  References: 
  (1) Akimov, A. V.; Trivedi, D.; Wang, L.; Prezhdo, O. V. Analysis of the Trajectory Surface Hopping Method from the Markov State Model Perspective. J. Phys. Soc. Jpn. 2015, 84, 094002.

*/



  int nst = Coeff.n_elts;
  MATRIX g(nst,nst);

  double norm; norm = (Coeff.H() * Coeff).get(0,0).real();  // <- this is the norm <PSI|PSI>

  // Calculate the hopping probabilities
  for(int j=0;j<nst;j++){

    double gjj = (std::conj(Coeff.get(j)) * Coeff.get(j)).real()/norm; // c_j^* * c_j
      
    for(int i=0;i<nst;i++){   g.set(i,j,gjj);  } // hopping from i to j


  }// for j

  return g;

}




MATRIX compute_hopping_probabilities_mssh(CMATRIX& Coeff, CMATRIX& Hvib, int use_boltz_factor,double T){

  const double kb = 3.166811429e-6; // Hartree/K
  int nst = Coeff.n_elts;
  MATRIX g(nst,nst);
  double norm; norm = (Coeff.H() * Coeff).get(0,0).real();  // <- this is the norm <PSI|PSI>
  double argg;

   // Calculate the hopping probabilities^M
   for(int j=0;j<nst;j++){
       double gjj = (std::conj(Coeff.get(j)) * Coeff.get(j)).real()/norm; // c_j^* * c_j
	 for(int i=0;i<nst;i++){   g.set(i,j,gjj);  } // hopping from i to j
   }

   if(use_boltz_factor){
     for(int i=0;i<nst;i++){
       for(int j=0;j<nst;j++){ 
	 if(Hvib.get(j,j).real()>Hvib.get(i,i).real()){
	   argg = -(Hvib.get(j,j).real() - Hvib.get(i,i).real())/(kb*T);
	   if(argg<500.0){ g.set(i,j, g.get(i,j) * std::exp(argg)); }
	 }
       }
     }
   }// if use_boltz_factor

   return g;
}

void compute_hopping_probabilities_mssh(Nuclear* mol, Electronic* el, Hamiltonian* ham, MATRIX* g,
                                        double dt, int use_boltz_factor,double T){
/**
  \brief Compute the MSSH surface hopping probabilities for the trajectory described by mol, el, and ham
  \param[in] mol Describes the nuclear DOF
  \param[in] el Describes electronic DOF
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[out] g The matrix of hopping probabilities
  \param[in] dt Time duration of nuclear propagation step
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Assume ham internal processings are done outside
  Assume matrix g is allocated

  Abbreviation: MSSH - Markov state surface hopping
  References: 
  (1) Akimov, A. V.; Trivedi, D.; Wang, L.; Prezhdo, O. V. Analysis of the Trajectory Surface Hopping Method from the Markov State Model Perspective. J. Phys. Soc. Jpn. 2015, 84, 094002.

*/



  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;


  // Verify if memory is allocated, if not generate
  if(g->n_elts!=el->nstates*el->nstates){cout<<"Matrix g is not allocated\n"; exit(0);  }


  // Now calculate the hopping probabilities
  for(j=0;j<el->nstates;j++){
    double gjj = (el->q[j]*el->q[j] + el->p[j]*el->p[j]); // c_j^* * c_j 

    for(i=0;i<el->nstates;i++){   g->set(i,j,gjj);   }// for i

  }// for j


}// compute probabilities mssh

void compute_hopping_probabilities_mssh(Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
                                        double dt, int use_boltz_factor,double T){
/**
  \brief Compute the MSSH surface hopping probabilities for the trajectory described by mol, el, and ham - Python-friendly
  \param[in] mol Describes the nuclear DOF
  \param[in] el Describes electronic DOF
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[out] g The matrix of hopping probabilities
  \param[in] dt Time duration of nuclear propagation step
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Assume ham internal processings are done outside
  Assume matrix g is allocated

  Abbreviation: MSSH - Markov state surface hopping
  References: 
  (1) Akimov, A. V.; Trivedi, D.; Wang, L.; Prezhdo, O. V. Analysis of the Trajectory Surface Hopping Method from the Markov State Model Perspective. J. Phys. Soc. Jpn. 2015, 84, 094002.

*/


  compute_hopping_probabilities_mssh(&mol, &el, &ham, &g, dt, use_boltz_factor, T);

}

void compute_hopping_probabilities_mssh(Ensemble& ens, int i, MATRIX& g, double dt, int use_boltz_factor,double T){
/**
  \brief Compute the MSSH surface hopping probabilities for a selected trajectory of an ensemble - Python friendly
  \param[in] ens Describes the ensemble of trajectories
  \param[in] i Is the index of the trajectory of interest
  \param[out] g The matrix of hopping probabilities for this tranjectory
  \param[in] dt Time duration of nuclear propagation step
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Assume ham internal processings are done outside
  Assume matrix g is allocated

  Abbreviation: MSSH - Markov state surface hopping
  References: 
  (1) Akimov, A. V.; Trivedi, D.; Wang, L.; Prezhdo, O. V. Analysis of the Trajectory Surface Hopping Method from the Markov State Model Perspective. J. Phys. Soc. Jpn. 2015, 84, 094002.

*/


  compute_hopping_probabilities_mssh(&ens.mol[i], &ens.el[i], ens.ham[i], &g, dt, use_boltz_factor, T);

}


int hop(int initstate, MATRIX& g, double ksi){
/** 
  \brief Attempts a stochastic hop from the initial state "initstate"
  \param[in] initstate The index of the state from which we try to hop out 
  \param[in] g The hopping probabilities matrix (type MATRIX)
  \param[in] ksi A random number that determines the outcome of the "hop" procedure

  Returned value: the index of the state to which we have hopped
*/

  int nstates = g.n_cols;
  double left, right; left = right = 0.0;
  int finstate = initstate;

  for(int i=0;i<nstates;i++){
    if(i==0){left = 0.0; right = g.get(initstate,i); }
    else{  left = right; right = right + g.get(initstate,i); }
 
    if((left<ksi) && (ksi<=right)){  finstate = i;  }
  }

  return finstate;

}// hop



void hop(int& initstate, Nuclear* mol, Hamiltonian* ham, double ksi, MATRIX* g, int do_rescaling, int rep, int do_reverse){
/** 
  \brief Do actual hop from the state initstate 
  \param[in,out] initstate The state from which we try to hop out - it will also be updated after the hop has happened
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in] ksi A random number that determines the outcome of the "hop" procedure
  \param[in] g The hopping probabilities matrix
  \param[in] do_rescaling The flag to turn on/off CPA: 0 - no velocity rescaling (CPA, no back-reaction),
  in this case one should use Boltzmann factor (consider use_boltz_factor when computing the hopping probability matrix, g)
  1 - do rescaling (back-reaction), in this case it would be wrong to use Boltzmann factor
  \param[in] rep Selects the used representation:  0 - for diabatic, 1 - for adiabatic
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.
*/

  int nstates = g->n_cols;
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
/** 
  \brief Do actual hop from the state initstate - Python-friendly
  \param[in] initstate The state from which we try to hop out - it will also be updated after the hop has happened
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in] ksi A random number that determines the outcome of the "hop" procedure
  \param[in] g The hopping probabilities matrix
  \param[in] do_rescaling The flag to turn on/off CPA: 0 - no velocity rescaling (CPA, no back-reaction),
  in this case one should use Boltzmann factor (consider use_boltz_factor when computing the hopping probability matrix, g)
  1 - do rescaling (back-reaction), in this case it would be wrong to use Boltzmann factor
  \param[in] rep Selects the used representation:  0 - for diabatic, 1 - for adiabatic
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  The function returns the index of the final state (new or old).
*/


  int res = initstate; 
  hop(res, &mol, &ham, ksi, &g, do_rescaling, rep, do_reverse);

  return res;

}

int hop(int initstate, Ensemble& ens, int i, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse){
/** 
  \brief Do actual hop from the state initstate of the i-th trajectory of an ensemble - Python-friendly
  \param[in] initstate The state from which we try to hop out - it will also be updated after the hop has happened
  \param[in] i The index of the trajectory of interest
  \param[in,out] ens Describes the ensemble of trajectories which we propagate (including hops)
  \param[in] ksi A random number that determines the outcome of the "hop" procedure
  \param[in] g The hopping probabilities matrix
  \param[in] do_rescaling The flag to turn on/off CPA: 0 - no velocity rescaling (CPA, no back-reaction),
  in this case one should use Boltzmann factor (consider use_boltz_factor when computing the hopping probability matrix, g)
  1 - do rescaling (back-reaction), in this case it would be wrong to use Boltzmann factor
  \param[in] rep Selects the used representation:  0 - for diabatic, 1 - for adiabatic
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  The function returns the index of the final state (new or old).
*/


  int res = initstate; 
  hop(res, &ens.mol[i], ens.ham[i], ksi, &g, do_rescaling, rep, do_reverse);

  return res;

}



void rescale_velocities_adiabatic(Nuclear* mol, Hamiltonian* ham, int& new_st,int& old_st, int do_reverse){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in,out] new_st The index of the new state: this is a truly new state if the attempted hop was successfull, or just can be 
             an old state, otherwise
  \param[in] old_st The index of the old state (from which we try to hop)
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  This verions implies that the adiabatic representation is used
*/

  int st;

  //cout<<"in rescale_velocities_adiabatic\n";

  if(new_st!=old_st){

    int final_st = old_st;  // no hopping by default

    // According to Fabiano
    double a_ij = 0.0;
    double b_ij = 0.0;

    for(int k=0;k<mol->nnucl;k++){

      double D = ham->D(old_st,new_st,k).real(); // derivative coupling w.r.t. nuclear DOF k
      cout<<"D= "<<D<<endl;

      a_ij += 0.5*(D*D / mol->mass[k]); 
      b_ij += (D*mol->p[k])/mol->mass[k];
    }
    double det = b_ij*b_ij + 4.0*a_ij*(ham->Hvib(old_st,old_st).real() - ham->Hvib(new_st,new_st).real());

    cout<<"det= "<<det<<" a_ij= "<<a_ij<<"  b_ij= "<<b_ij<<endl;

    // Calculate the scaling factor and new state
    double gamma_ij = 0.0;

    if(det<0.0){

      if(do_reverse){     gamma_ij = b_ij / a_ij;}
      else{ gamma_ij = 0.0;  }

      final_st = old_st; // # hop does not occur - frustrated hop

    }
    else{
      cout<<"gamma_ij= "<<gamma_ij<<endl;

      if(b_ij<0){ gamma_ij = 0.5*(b_ij + sqrt(det))/a_ij; }
      else{       gamma_ij = 0.5*(b_ij - sqrt(det))/a_ij; }
      final_st = new_st;
    }

    //Rescale velocities and do the hop
    for(int k=0;k<mol->nnucl;k++){ 
      double D = ham->D(old_st,new_st,k).real(); 
      cout<<"D= "<<D<<endl;
      mol->p[k] = mol->p[k] - gamma_ij * D; 
    }

    st = final_st;

  }
  else{ st = old_st; }

  new_st = st;


} // rescale velocities adiabatic

int rescale_velocities_adiabatic(Nuclear& mol, Hamiltonian& ham, int old_st, int do_reverse){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another - Python-friendly
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in] old_st The index of the old state (from which we try to hop)
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  The function returns the index of the new state: this is a truly new state if the attempted hop was successfull, or just can be 
  an old state, otherwise

  This verions implies that the adiabatic representation is used
*/


  int new_st = old_st;
  rescale_velocities_adiabatic(&mol, &ham, new_st, old_st, do_reverse);
  return new_st;

}


void rescale_velocities_diabatic(Nuclear* mol, Hamiltonian* ham, int& new_st,int& old_st){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in,out] new_st The index of the new state: this is a truly new state if the attempted hop was successfull, or just can be 
             an old state, otherwise
  \param[in] old_st The index of the old state (from which we try to hop)

  This verions implies that the diabatic representation is used

  In this case, derivative couplings are zero by definition of diabatic (position-independent) states
  so one can not rescale velocities along directions of derivative couplings, since there is no such directions
  We just scale velocities uniformly - based on energy conservation principle
  In principle, this rescaling procedure can be applied to surface hopping scheme in adiabatic basis, too
*/


  double T_i = compute_kinetic_energy(mol); // initial kinetic energy
  double E_i = ham->Hvib(old_st,old_st).real();// initial potential energy
  double E_f = ham->Hvib(new_st,new_st).real();// final potential energy
  
  double T_f = T_i + E_i - E_f;             // predicted final kinetic energy

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
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another - Python-friendly
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in] old_st The index of the old state (from which we try to hop)

  The function returns the index of the new state: this is a truly new state if the attempted hop was successfull, or just can be 
  an old state, otherwise

  This verions implies that the diabatic representation is used

  In this case, derivative couplings are zero by definition of diabatic (position-independent) states
  so one can not rescale velocities along directions of derivative couplings, since there is no such directions
  We just scale velocities uniformly - based on energy conservation principle
  In principle, this rescaling procedure can be applied to surface hopping scheme in adiabatic basis, too
*/


  int new_st = old_st;
  rescale_velocities_diabatic(&mol, &ham, new_st, old_st);
  return new_st;

}



int ida(CMATRIX& Coeff, int old_st, int new_st, double E_old, double E_new, double T, double ksi){
/**
  \brief Instantaneous decoherence at attempted hops (ID-A)

  \param[in,out] Coeff The matrix of coefficients (amplitudes of the basis excited states). This matrix is modified at every
  decoherence event
  \param[in] old_st The old state before an attempted hop
  \param[in] new_st The new state after an attempted hop (no matter what the algorithm - FSSH, GFSH, MSSH, or something else)
  \param[in] E_old The energy corresponding to the old state
  \param[in] E_new The energy corresponding to the new state
  \param[in] T is the temperature of the system
  \param[in] ksi a random number from a uniform distribution of the [0,1] interval - needed to decide the hop acceptance/decoherence

  The function returns the index of the new state after hop rejection/decoherence criteria
  The function also modifies the amplitudes of the coherent superposition

*/

  const double kb = 3.166811429e-6;  // Boltzmann constant: Hartree/K
  int istate = old_st;               // by default the resulting state is the old state
  
  if(new_st != old_st){  // attempted hop

    // Now apply energy considerations
    double dE = (E_new - E_old);
    double boltz_f = 1.0;

    if(dE>0.0){
      double argg = dE/(kb*T);
      if(argg > 50.0){ boltz_f = 0.0; }
      else{            boltz_f = exp(-argg); }

      if(ksi<boltz_f){
        istate = new_st;  // accepted hop

        // Collapse the wavefunction onto the new state
        Coeff *= 0.0;
        Coeff.set(new_st, 1.0, 0.0);
      }
      else{
        // Unsuccessful hop - collapse wfc back to the original state
        Coeff *= 0.0;
        Coeff.set(old_st, 1.0, 0.0);
      }
    }// dE > 0
  }// new_st != old_st

  return istate;

}


MATRIX coherence_intervals(CMATRIX& Coeff, MATRIX& rates ){
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

  CMATRIX* denmat; denmat = new CMATRIX(nstates, nstates);   
  *denmat = (Coeff * Coeff.H() ).conj();

  MATRIX tau_m(nstates, 1); tau_m *= 0.0;
  

  for(int i=0;i<nstates;i++){

    double summ = 0.0;
    for(int j=0;j<nstates;j++){

      if(j!=i){
        summ += denmat->get(j,j).real() * rates.get(i,j); 
      }// if

    }// for j

    if(summ>0.0){   tau_m.set(i, 1.0/summ); }
    else        {   tau_m.set(i, 1.0e+100); } // infinite coherence interval
    
     
  }// for i

  delete denmat;

  return tau_m;
}



int dish(Electronic& el, MATRIX& t_m, const MATRIX& tau_m, const CMATRIX& Hvib, int use_boltz_flag, double Ekin, double T, double ksi1, double ksi2){
/**
  \brief Decoherence-induced surface hopping (DISH)

  Reference: Jaeger, H. M.; Fischer, S.; Prezhdo, O. V. Decoherence-Induced Surface Hopping. J. Chem. Phys. 2012, 137, 22A545.

  \param[in,out] el Electronic object, containing the info about electronic amplitudes
  \param[in,out] t_m A matrix N x 1 of the times each state resides in a coherence interval
   (since the last decoherence event)
  \param[in] tau_m A matrix N x 1 of the coherence intervals for each electronic state
  \param[in] Hvib The matrix of vibronic Hamiltonian - we need the energies
  \param[in] use_boltz_flag  if set to 1, the hopping probabilities will be re-scaled by the 
  Boltzmann factor. This is need in neglect of back-reaction approximation (NBRA) calculations, 
  when no 
  \param[in] Ekin Kinetic energy of nuclei
  \param[in] T is the temperature of the system
  \param[in] ksi1 a is random number from a uniform distribution on the [0,1] interval
  \param[in] ksi2 a is random number from a uniform distribution on the [0,1] interval

  The function modifies  el and t_m variables
  Returns: the index of the electronic state after potential hop
*/

  const double kb = 3.166811429e-6; // Hartree/K   
  int i,j; 
  int has_decoherence = 0; /// set to 1 if we have already encountered a decoherence event

  for(i=0;i<el.nstates && !has_decoherence;i++){

    /// The state i has evolved coherently for longer than the coherence interval
    /// so it has to experience a decoherence event 
    if(t_m.get(i) >= tau_m.get(i) ) { 

      /// There are essentially two outcomes when the decoherence takes place:
      /// One: we collapse the wavefunction onto the state i with the probability 
      /// given by the population of that state

      if(ksi1 < el.rho(i,i).real()){ 

        /// Now, lets determine if the hop is possible based on the energy conservation
        /// considerations 



        int can_hop = 0;

        double E_i = Hvib.get(el.istate, el.istate).real();/// initial potential energy
        double E_f = Hvib.get(i,i).real();                 /// proposed potential energy
        double dE = E_f - E_i;

        /// In leu of hop rejection use Boltzmann factors: for NBRA simulations
        if(use_boltz_flag==1){
          double bf = 1.0;
          if(dE>0){  bf = exp(-dE/(kb*T)); }  /// hop to higher energy state is difficult
          if(ksi2<bf){ can_hop = 1; }  /// hop is allowed thermally

        }
        /// Regular energy-conservation based criterion
        else{   
          /// Predicted final kinetic energy is positive - hop is possible
          if(Ekin - dE > 0.0){  can_hop = 1; }
        }
          
        /// Now, decide about decoherence         
        if(can_hop){     el.collapse(i, 1);  } /// Here is the actuall collapse
        else{ el.project_out(i); }

      }

      /// Second: project the system out of that state otherwise
      else{  el.project_out(i);   }


      /// Reset the time axis for state i (only for this state)
      /// other states still reside in a coherent superposition
      t_m.set(i, 0.0);

      /// Set the flag that we have attempted a decoherence event
      /// so we done with DISH at this point in time
      has_decoherence = 1;


    }// t_m[i]>=1.0/tau_m
      
  }// for i

  return el.istate;

} // dish


int dish(Electronic& el, Nuclear& mol, Hamiltonian& ham, MATRIX& t_m, const MATRIX& tau_m, int use_boltz_flag, double T, double ksi1, double ksi2){
/**
  \brief Decoherence-induced surface hopping (DISH) - overloaded version

  Reference: Jaeger, H. M.; Fischer, S.; Prezhdo, O. V. Decoherence-Induced Surface Hopping. J. Chem. Phys. 2012, 137, 22A545.

  \param[in,out] el Electronic object, containing the info about electronic amplitudes
  \param[in,out] t_m A matrix N x 1 of the times each state resides in a coherence interval
   (since the last decoherence event)
  \param[in] tau_m A matrix N x 1 of the coherence intervals for each electronic state
  \param[in] Hvib The matrix of vibronic Hamiltonian - we need the energies
  \param[in] use_boltz_flag  if set to 1, the hopping probabilities will be re-scaled by the 
  Boltzmann factor. This is need in neglect of back-reaction approximation (NBRA) calculations, 
  when no 
  \param[in] Ekin Kinetic energy of nuclei
  \param[in] T is the temperature of the system
  \param[in] ksi1 a is random number from a uniform distribution on the [0,1] interval
  \param[in] ksi2 a is random number from a uniform distribution on the [0,1] interval

  The function modifies  el and t_m variables
  Returns: the index of the electronic state after potential hop
*/

  const double kb = 3.166811429e-6; // Hartree/K   
  int i,j; 
  int has_decoherence = 0; /// set to 1 if we have already encountered a decoherence event

  for(i=0;i<el.nstates && !has_decoherence;i++){

    /// The state i has evolved coherently for longer than the coherence interval
    /// so it has to experience a decoherence event 
    if(t_m.get(i) >= tau_m.get(i) ) { 

      /// There are essentially two outcomes when the decoherence takes place:
      /// One: we collapse the wavefunction onto the state i with the probability 
      /// given by the population of that state

      if(ksi1 < el.rho(i,i).real()){ 

        /// Now, lets determine if the hop is possible based on the energy conservation
        /// considerations 



        int can_hop = 0;

        double E_i = ham.Hvib(el.istate, el.istate).real();/// initial potential energy
        double E_f = ham.Hvib(i,i).real();                 /// proposed potential energy
        double dE = E_f - E_i;

        /// In leu of hop rejection use Boltzmann factors: for NBRA simulations
        if(use_boltz_flag==1){
          double bf = 1.0;
          if(dE>0){  bf = exp(-dE/(kb*T)); }  /// hop to higher energy state is difficult
          if(ksi2<bf){ can_hop = 1; }  /// hop is allowed thermally

        }
        /// Regular energy-conservation based criterion
        else{   
          /// Predicted final kinetic energy is positive - hop is possible
          double Ekin = compute_kinetic_energy(mol); // initial kinetic energy
          if(Ekin - dE > 0.0){  can_hop = 1; }
        }
          
        /// Now, decide about decoherence         
        if(can_hop){     el.collapse(i, 1);  } /// Here is the actuall collapse
        else{ el.project_out(i); }

      }

      /// Second: project the system out of that state otherwise
      else{  el.project_out(i);   }


      /// Reset the time axis for state i (only for this state)
      /// other states still reside in a coherent superposition
      t_m.set(i, 0.0);

      /// Set the flag that we have attempted a decoherence event
      /// so we done with DISH at this point in time
      has_decoherence = 1;


    }// t_m[i]>=1.0/tau_m
      
  }// for i

  return el.istate;

} // dish




}// namespace libdyn
}// liblibra

