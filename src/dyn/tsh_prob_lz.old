/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file tsh_prob_lz.cpp
  \brief The file implements the Landau-Zener surface hopping probabilities
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;
using namespace libnhamiltonian;
namespace bp = boost::python;


/// libdyn namespace
namespace libdyn{



MATRIX compute_hopping_probabilities_lz(nHamiltonian* ham, int rep, MATRIX& p, const MATRIX& invM, MATRIX& prev_ham_dia){
/**
  \brief This function computes the surface hopping probabilities according to Landau-Zener formula. 

  See more details in:
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061–1071. - the original paper
  (2) Belyaev, A. K.; Lebedev, O. V. Nonadiabatic nuclear dynamics of atomic collisions based on branching classical trajectories
  Phys. Rev. A 2011, 84, 014701

  \param[in] ham Is the nHamiltonian object that organizes energy-related calculations
  \param[in] rep index selecting the type of representation to be used: 0 - diabatic, 1 - adiabatic
  \param[in] p [ndof x 1] or [ndof x ntraj] matrix of nuclear momenta
  \param[in] invM [ndof x 1] matrix of inverse nuclear masses
  \param[in] prev_ham_dia [nstates x nstates] the matrix of previous diabatic Hamiltonian - needed to locate the diabatic crossing 

  Returns: 
      A matrix with the hopping probabilities between all pairs of states is returned
      Convention: P(i,j) is the probability for the i->j transition

*/
  int i,j,k;
  int nstates;
  int ndof = ham->nnucl;

  // Determine the dimensions
  if(rep==0){ nstates = ham->ndia;  } 
  else if(rep==1){  nstates = ham->nadi;  }

  MATRIX g(nstates,nstates);

  // Diabatic representation
  if(rep==0){     
    
    // Get denominator
    MATRIX dham(nstates,nstates);
    MATRIX dHv(nstates,nstates); // |H'_ii - H'_jj| * v
    
    for(k=0;k<ndof;k++){  
    
      dham = ham->get_d1ham_dia(k).real() * p.get(k,0) * invM.get(k,0);
    
      for(i=0;i<nstates;i++){
        for(j=0;j<nstates;j++){
    
          dHv.add(i,j, fabs(dham.get(i,i) - dham.get(j,j)) );
    
        }// for j
      }// for i
    }// for k
    
    
    // Now calculate the hopping probabilities
    MATRIX ham_dia(nstates,nstates);    
    ham_dia = ham->get_ham_dia().real();
    
    // Off-diagonal elements
    for(i=0;i<nstates;i++){
      for(j=i+1;j<nstates;j++){
    
        double dh = ham_dia.get(i,i) - ham_dia.get(j,j);
        double dh_prev = prev_ham_dia.get(i,i) - prev_ham_dia.get(j,j);

        if(dh * dh_prev < 0.0){


          double h_ij = ham_dia.get(i,j);
          double g_ij = exp(-2.0*M_PI*h_ij*h_ij / dHv.get(i,j) );
    
          g.set(i,j,g_ij);
          g.set(j,i,g_ij);
        }
        else{ g.set(i,j, 0.0);  g.set(j,i, 0.0); }
    
      }// for j
    }// for i
    
  }// rep == 0 diabatic

  // Adiabatic representation
  else if(rep==1){

    // Update adiabatic NACs
    ham->compute_nac_adi(p, invM);

    // Now calculate the hopping probabilities
    MATRIX ham_dia(nstates,nstates);    
    MATRIX ham_adi(nstates,nstates);
    MATRIX nac_adi(nstates,nstates);

    ham_dia = ham->get_ham_dia().real();    
    ham_adi = ham->get_ham_adi().real();
    nac_adi = ham->get_nac_adi().real();
    
    // Off-diagonal elements
    for(i=0;i<nstates;i++){
      for(j=i+1;j<nstates;j++){

        double g_ij = 0.0;

        double dh = ham_dia.get(i,i) - ham_dia.get(j,j);
        double dh_prev = prev_ham_dia.get(i,i) - prev_ham_dia.get(j,j);

        if(dh * dh_prev < 0.0){        
    
          double Z_ij = fabs(ham_adi.get(i,i) - ham_adi.get(j,j));        
          double nac_ij = fabs(nac_adi.get(i,j));
          if(nac_ij > 0.0){ g_ij = exp(-0.25*M_PI*Z_ij/nac_ij);     }

        }

        g.set(i,j,g_ij);
        g.set(j,i,g_ij);
        
    
      }// for j
    }// for i
   
  }// rep == 1 adiabatic


  // Diagonal elements - common for both reps
  for(i=0;i<nstates;i++){
    double sum = 0.0;      
  
    for(j=0;j<nstates;j++){
  
      if(j!=i){  sum += g.get(i,j);   }
  
    }// for i
  
    g.set(i,i, 1.0 - sum);
  }


  return g;

}// lz


MATRIX compute_hopping_probabilities_lz(nHamiltonian& ham, int rep, MATRIX& p, const MATRIX& invM, MATRIX& prev_ham_dia){
/**
  \brief This function computes the surface hopping probabilities according to Landau-Zener formula. 

  See more details in:
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061–1071. - the original paper
  (2) Belyaev, A. K.; Lebedev, O. V. Nonadiabatic nuclear dynamics of atomic collisions based on branching classical trajectories
  Phys. Rev. A 2011, 84, 014701

  \param[in] ham Is the nHamiltonian object that organizes energy-related calculations
  \param[in] rep index selecting the type of representation to be used: 0 - diabatic, 1 - adiabatic
  \param[in] p [ndof x 1] or [ndof x ntraj] matrix of nuclear momenta
  \param[in] invM [ndof x 1] matrix of inverse nuclear masses
  \param[in] prev_ham_dia [nstates x nstates] the matrix of previous diabatic Hamiltonian - needed to locate the diabatic crossing 

  Returns: 
      A matrix with the hopping probabilities between all pairs of states is returned
      Convention: P(i,j) is the probability for the i->j transition

*/

  return compute_hopping_probabilities_lz(&ham, rep, p, invM, prev_ham_dia);

}



}// namespace libdyn
}// liblibra

