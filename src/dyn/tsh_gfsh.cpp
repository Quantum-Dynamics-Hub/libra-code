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
  \file tsh_gfsh.cpp
  \brief The file implements the Global Flux Surface Hopping  (GFSH) - related algorithms
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

#include "Dynamics.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;
using namespace libhamiltonian;
namespace bp = boost::python;

/// libdyn namespace
namespace libdyn{




MATRIX compute_hopping_probabilities_gfsh(CMATRIX& Coeff, CMATRIX* Hvib, double dt, int use_boltz_factor,double T){

/**
  \brief Compute the GFSH surface hopping probabilities for a single trajectory

  \param[in] Coeff - [ndia x 1] or a [nadi x 1] matrix of electronic basis states amplitudes in a superposition
  \param[in] Hvib - [ndia x ndia] or a [nadi x nadi] vibronic Hamiltonian matrix
  \param[in] dt - the interval of the time step for which we compute the hopping probability
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Returns: A matrix with the hopping probabilities between all pairs of states is returned

  Abbreviation: GFSH - global flux surface hopping
  References: 
  (1) Wang, L.; Trivedi, D.; Prezhdo, O. V. Global Flux Surface Hopping Approach for Mixed Quantum-Classical Dynamics. J. Chem. Theory Comput. 2014, 10, 3598–3605.
  (2) Akimov, A. V. Libra: An Open-Source “Methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626–1649.


*/

  int nstates = Coeff.n_rows;
  MATRIX g(nstates,nstates);

  CMATRIX* denmat; denmat = new CMATRIX(nstates, nstates);   
  *denmat = (Coeff * Coeff.H() ).conj();

  CMATRIX* denmat_dot; denmat_dot = new CMATRIX(nstates, nstates);   
  *denmat_dot = ((*denmat) *  (*Hvib).conj() - (*Hvib) * (*denmat)) * complex<double>(0.0, 1.0);


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

/*
        if(use_boltz_factor){

          if(ham->H(j,j).real() > ham->H(i,i).real()){
            argg = -(ham->H(j,j).real() - ham->H(i,i).real())/(kb*T);        
            if(argg<500.0){ g_ij = g_ij * std::exp(argg); }
          }

        }// if use_boltz_factor

*/

        sumg += g.get(i,j);

      }
    }// for j

    g.set(i,i, 1.0 - sumg);  // probability to stay in state i

  }// for i

  delete denmat;
  delete denmat_dot;

}// gfsh



MATRIX compute_hopping_probabilities_gfsh(CMATRIX& Coeff, CMATRIX* Hvib, double dt){
/**
  \brief Compute the GFSH surface hopping probabilities for a single trajectory 

  \param[in] Coeff - [ndia x 1] or a [nadi x 1] matrix of electronic basis states amplitudes in a superposition
  \param[in] Hvib - [ndia x ndia] or a [nadi x nadi] vibronic Hamiltonian matrix
  \param[in] dt - the interval of the time step for which we compute the hopping probability

  Returns: A matrix with the hopping probabilities between all pairs of states is returned

  Abbreviation: GFSH - global flux surface hopping
  References: 
  (1) Wang, L.; Trivedi, D.; Prezhdo, O. V. Global Flux Surface Hopping Approach for Mixed Quantum-Classical Dynamics. J. Chem. Theory Comput. 2014, 10, 3598–3605.
  (2) Akimov, A. V. Libra: An Open-Source “Methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626–1649.

*/

  return compute_hopping_probabilities_gfsh(Coeff, Hvib, dt, 0, 1.0);

}


MATRIX compute_hopping_probabilities_gfsh(CMATRIX& Coeff, CMATRIX& Hvib, double dt, int use_boltz_factor,double T){
/**
  \brief See the description of the 
  compute_hopping_probabilities_gfsh(CMATRIX& Coeff, CMATRIX* Hvib, double dt, int use_boltz_factor,double T) function
*/

  return compute_hopping_probabilities_gfsh(Coeff, &Hvib, dt, use_boltz_factor, T);

}// gfsh


MATRIX compute_hopping_probabilities_gfsh(CMATRIX& Coeff, CMATRIX& Hvib, double dt){
/**
  \brief See the description of the 
  compute_hopping_probabilities_gfsh(CMATRIX& Coeff, CMATRIX* Hvib, double dt) function
*/

  return compute_hopping_probabilities_gfsh(Coeff, &Hvib, dt, 0, 1.0);

}

MATRIX compute_hopping_probabilities_gfsh(CMATRIX& Coeff, nHamiltonian* ham, int rep, double dt, int use_boltz_factor,double T){
/**
  \brief Compute the GFSH surface hopping probabilities for a single trajectory. May be
  Boltzmann-corrected

  \param[in] Coeff - [ndia x 1] or a [nadi x 1] matrix describing electronic DOFs
  \param[in] ham - Is the nHamiltonian object that organizes energy-related calculations
  \param[in] rep - index selecting the type of representation to be used: 0 - diabatic, 1 - adiabatic
  \param[in] dt Time duration of nuclear propagation step
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Returns: A matrix with the hopping probabilities between all pairs of states is returned
  Abbreviation: GFSH - fewest switches surface hopping

  References: 
  (1) Wang, L.; Trivedi, D.; Prezhdo, O. V. Global Flux Surface Hopping Approach for Mixed Quantum-Classical Dynamics. J. Chem. Theory Comput. 2014, 10, 3598–3605.
  (2) Akimov, A. V. Libra: An Open-Source “Methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626–1649.

*/
  if(rep==0){
    return compute_hopping_probabilities_gfsh(Coeff, ham->hvib_dia, dt, use_boltz_factor, T);
  }
  else if(rep==1){
    return compute_hopping_probabilities_gfsh(Coeff, ham->hvib_adi, dt, use_boltz_factor, T);
  }
  else{
    cout<<"Error in compute_hopping_probabilities_gfsh(...): representation rep = "<<rep<<" is not known\n";
    cout<<"Exiting...\n";
    exit(0);
  }

}

MATRIX compute_hopping_probabilities_gfsh(CMATRIX& Coeff, nHamiltonian& ham, int rep, double dt, int use_boltz_factor,double T){
/**
  \brief See the description of the
  MATRIX compute_hopping_probabilities_gfsh(CMATRIX& Coeff, nHamitlonian* ham, int rep, double dt, int use_boltz_factor,double T)
  function

*/
  if(rep==0){
    return compute_hopping_probabilities_gfsh(Coeff, ham.hvib_dia, dt, use_boltz_factor, T);
  }
  else if(rep==1){
    return compute_hopping_probabilities_gfsh(Coeff, ham.hvib_adi, dt, use_boltz_factor, T);
  }
  else{
    cout<<"Error in compute_hopping_probabilities_gfsh(...): representation rep = "<<rep<<" is not known\n";
    cout<<"Exiting...\n";
    exit(0);
  }

}





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
  (2) Akimov, A. V. Libra: An Open-Source “Methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626–1649.

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




}// namespace libdyn
}// liblibra

