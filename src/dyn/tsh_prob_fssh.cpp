/*********************************************************************************
* Copyright (C) 2015-2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file tsh_prob_fssh.cpp
  \brief The file implements the Fewest Switches Surface Hopping  (FSSH) - related algorithms
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;
using namespace libhamiltonian;
namespace bp = boost::python;


/// libdyn namespace
namespace libdyn{




MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, CMATRIX* Hvib, double dt, int use_boltz_factor,double T){
/**
  \brief This function computes the surface hopping probabilities according to Tully's FSSH prescription. 
  The surface-hopping probabilities may be Boltzmann-corrected

  See more details in:
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061–1071. - the original paper
  (2) Fabiano, E.; Keal, T. W.; Thiel, W. Implementation of Surface Hopping Molecular Dynamics Using Semiempirical Methods. Chem. Phys. 2008, 349, 334–347.
  Here, we generalized the formula, so it works equally well for both diabatic and adiabatic representations
  (3) Akimov, A. V. Libra: An Open-Source “Methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626–1649.

  \param[in] Coeff - [ndia x 1] or a [nadi x 1] matrix of electronic basis states amplitudes in a superposition
  \param[in] Hvib - [ndia x ndia] or a [nadi x nadi] vibronic Hamiltonian matrix
  \param[in] dt - the interval of the time step for which we compute the hopping probability
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Returns: A matrix with the hopping probabilities between all pairs of states is returned

*/
  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;

  int nstates = Coeff.n_rows;
  MATRIX g(nstates,nstates);

  CMATRIX denmat(nstates, nstates);   
  
  denmat = Coeff * Coeff.H();

  // Now calculate the hopping probabilities
  for(i=0;i<nstates;i++){
    sum = 0.0;
    double a_ii = denmat.get(i,i).real(); 

    for(j=0;j<nstates;j++){

      if(i!=j){ 
        /**
          dc/dt = -(i/hbar) * Hvib * c
          (dc/dt)^+ = i/hbar * c^+ * Hvib^+

          rho = c * c^+

          Then
          drho/dt = i/hbar * (rho*Hvib - Hvib*rho)

          The diagonal element:
           
          drho_ii/dt = (i/hbar) *sum_a { rho_ia * Hvib_ai - Hvib_ia * rho_ai}

          Then:  P(i->*) = -(drho_ii / rho_ii ) * dt  - probability of leaving state i

          Then:  P(i->a) = - (dt/rho_ii) * Re[ (i/hbar) *sum_a { rho_ia * Hvib_ai - Hvib_ia * rho_ai} ] 

          = (dt/(hbar*rho_ii)) * Im[ rho_ia * Hvib_ai - Hvib_ia * rho_ai ] 

          Or:

          P(i->j) = (dt/(hbar*rho_ii)) * Im[ rho_ij * Hvib_ji - Hvib_ij * rho_ji ] 

        */

        double imHaij = ( denmat.get(i,j) * Hvib->get(j,i) - Hvib->get(i,j) * denmat.get(j,i) ).imag(); 

        if(a_ii<1e-8){ g_ij = 0.0; }  // avoid division by zero
        else{
          g_ij = dt*imHaij/a_ii;  // This is a general case -

          if(use_boltz_factor){

            if(Hvib->get(j,j).real() > Hvib->get(i,i).real()){
              argg = -(Hvib->get(j,j).real() - Hvib->get(i,i).real())/(kb*T);        
              if(argg<500.0){ g_ij = g_ij * std::exp(argg); }
            }

          }// if use_boltz_factor


          if(g_ij<0.0){  g_ij = 0.0; }

        }// else

        g.set(i,j,g_ij);
        sum = sum + g_ij;
      }
      else{ g.set(i,j,0.0); }

    }// for j

    g.set(i,i,1.0 - sum);

  }// for i

  return g;

}// fssh


MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, CMATRIX* Hvib, double dt){
/**
  \brief This function computes the surface hopping probabilities according to Tully's FSSH prescription. 

  See more details in:
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061–1071. - the original paper
  (2) Fabiano, E.; Keal, T. W.; Thiel, W. Implementation of Surface Hopping Molecular Dynamics Using Semiempirical Methods. Chem. Phys. 2008, 349, 334–347.
  Here, we generalized the formula, so it works equally well for both diabatic and adiabatic representations
  (3) Akimov, A. V. Libra: An Open-Source “Methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626–1649.

  \param[in] Coeff - [ndia x 1] or a [nadi x 1] matrix of electronic basis states amplitudes in a superposition
  \param[in] Hvib - [ndia x ndia] or a [nadi x nadi] vibronic Hamiltonian matrix
  \param[in] dt - the interval of the time step for which we compute the hopping probability

  Returns: A matrix with the hopping probabilities between all pairs of states is returned

*/

  return compute_hopping_probabilities_fssh(Coeff, Hvib, dt, 0, 1.0);

}


MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, CMATRIX& Hvib, double dt, int use_boltz_factor,double T){
/**
  \brief See the description of the 
  compute_hopping_probabilities_fssh(CMATRIX& Coeff, CMATRIX* Hvib, double dt, int use_boltz_factor,double T) function
*/

  return compute_hopping_probabilities_fssh(Coeff, &Hvib, dt, use_boltz_factor, T);

}// fssh


MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, CMATRIX& Hvib, double dt){
/**
  \brief See the description of the 
  compute_hopping_probabilities_fssh(CMATRIX& Coeff, CMATRIX* Hvib, double dt) function
*/

  return compute_hopping_probabilities_fssh(Coeff, &Hvib, dt, 0, 1.0);

}


MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, nHamiltonian* ham, int rep, double dt, int use_boltz_factor,double T){
/**
  \brief Compute the FSSH surface hopping probabilities for a single trajectory. May be
  Boltzmann-corrected

  \param[in] Coeff - [ndia x 1] or a [nadi x 1] matrix describing electronic DOFs
  \param[in] ham - Is the nHamiltonian object that organizes energy-related calculations
  \param[in] rep - index selecting the type of representation to be used: 0 - diabatic, 1 - adiabatic
  \param[in] dt Time duration of nuclear propagation step
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Returns: A matrix with the hopping probabilities between all pairs of states is returned
  Abbreviation: FSSH - fewest switches surface hopping

  References: 
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061–1071. - the original paper
  (2) Fabiano, E.; Keal, T. W.; Thiel, W. Implementation of Surface Hopping Molecular Dynamics Using Semiempirical Methods. Chem. Phys. 2008, 349, 334–347.
  Here, we generalized the formula, so it works equally well for both diabatic and adiabatic representations
  (3) Akimov, A. V. Libra: An Open-Source “Methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626–1649.

*/
  if(rep==0){
    return compute_hopping_probabilities_fssh(Coeff, ham->hvib_dia, dt, use_boltz_factor, T);
  }
  else if(rep==1){
    return compute_hopping_probabilities_fssh(Coeff, ham->hvib_adi, dt, use_boltz_factor, T);
  }
  else{
    cout<<"Error in compute_hopping_probabilities_fssh(...): representation rep = "<<rep<<" is not known\n";
    cout<<"Exiting...\n";
    exit(0);
  }

}


MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, nHamiltonian* ham, int rep, double dt){
/**
  \brief See the description of 
  MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, nHamiltonian* ham, int rep, double dt, int use_boltz_factor,double T)
  function
*/

  return compute_hopping_probabilities_fssh(Coeff, ham, rep, dt, 0, 300.0);
}

MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, nHamiltonian& ham, int rep, double dt, int use_boltz_factor,double T){
/**
  \brief See the description of the
  MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, nHamitlonian* ham, int rep, double dt, int use_boltz_factor,double T)
  function

*/
  if(rep==0){
    return compute_hopping_probabilities_fssh(Coeff, ham.hvib_dia, dt, use_boltz_factor, T);
  }
  else if(rep==1){
    return compute_hopping_probabilities_fssh(Coeff, ham.hvib_adi, dt, use_boltz_factor, T);
  }
  else{
    cout<<"Error in compute_hopping_probabilities_fssh(...): representation rep = "<<rep<<" is not known\n";
    cout<<"Exiting...\n";
    exit(0);
  }

}


MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, nHamiltonian& ham, int rep, double dt){
/**
  \brief See the description of 
  MATRIX compute_hopping_probabilities_fssh(CMATRIX& Coeff, nHamiltonian& ham, int rep, double dt, int use_boltz_factor,double T)
  function
*/

  return compute_hopping_probabilities_fssh(Coeff, &ham, rep, dt, 0, 300.0);
}



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






}// namespace libdyn
}// liblibra

