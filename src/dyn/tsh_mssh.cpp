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
  \file tsh_mssh.cpp
  \brief The file implements the Markov State Surface Hopping  (MSSH) - related algorithms
    
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



  int nst = Coeff.n_rows;
  MATRIX g(nst,nst);

  double norm; norm = (Coeff.H() * Coeff).get(0,0).real();  // <- this is the norm <PSI|PSI>

  // Calculate the hopping probabilities
  for(int j=0;j<nst;j++){

    double gjj = (std::conj(Coeff.get(j)) * Coeff.get(j)).real()/norm; // c_j^* * c_j
      
    for(int i=0;i<nst;i++){   g.set(i,j,gjj);  } // hopping from i to j


  }// for j

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






}// namespace libdyn
}// liblibra

