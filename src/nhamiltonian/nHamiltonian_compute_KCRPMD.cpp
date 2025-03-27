/*********************************************************************************
* Copyright (C) 2018-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian_compute_KCRPMD.cpp
  \brief The file implements the calculations of the KC-RPMD (Kinetically Constrained Ring Polymer Molecular Dynamics)
  terms 
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <stdlib.h>
#endif 

#include "nHamiltonian.h"
#include "../math_meigen/libmeigen.h"

/// liblibra namespace
namespace liblibra{

/// libnhamiltonian namespace 
namespace libnhamiltonian{


using namespace liblinalg;
using namespace libmeigen;


double nHamiltonian::internal_ring_polymer_potential(const MATRIX& q, const MATRIX& invM, double beta){
/**
  Compute the ETHD energy

  q - is a ndof x ntraj matrix of coordinates
  invM - is a ndof x 1 matrix of inverse masses of all DOFs
  beta - the inverse temperature Boltzmann factor in atomic units
*/

  double en = 0.0;

  if(children.size()==0){
    en = 0.0;
  }
  else{ cout<<"Error in internal_ring_polymer_potential(), not implemented for quantum nuclei\n"; exit(0); }

  return en;
}


vector<MATRIX> nHamiltonian::generate_m_matrices(double beta){
/**
  Generate set of M matrices for each trajectory

  beta - the inverse temperature Boltzmann factor in atomic units
*/

  if(ham_dia_mem_status==0){ cout<<"Error in generate_m_matrices(): the diabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(ndia!=2){ cout<<"Error in generate_m_matrices(): implementation only for ndia=2\n"; exit(0); }

  vector<MATRIX> res;

  if(children.size()==0){
    res = vector<MATRIX>(1, MATRIX(2,2));
    double V0 = (this->ham_dia->get(0,0)).real(); 
    double V1 = (this->ham_dia->get(1,1)).real(); 
    double K = abs(this->ham_dia->get(0,1)); 
    res[0].set(0,0, exp(-beta * V0));
    res[0].set(0,1, -beta * K * exp(-beta * V0));
    res[0].set(1,0, -beta * K * exp(-beta * V1));
    res[0].set(1,1, exp(-beta * V1));
  }
  else{
    res = vector<MATRIX>(children.size(), MATRIX(2,2));
    for(int traj=0; traj<children.size(); traj++){
      double V0 = (children[traj]->ham_dia->get(0,0)).real(); 
      double V1 = (children[traj]->ham_dia->get(1,1)).real(); 
      double K = abs(children[traj]->ham_dia->get(0,1));
      res[traj].set(0,0, exp(-beta / children.size() * V0));
      res[traj].set(0,1, -beta / children.size() * K * exp(-beta / children.size() * V0));
      res[traj].set(1,0, -beta / children.size() * K * exp(-beta / children.size() * V1));
      res[traj].set(1,1, exp(-beta / children.size() * V1));
    }
  }
  return res;
}


double nHamiltonian::donor_electronic_potential(double beta){
/**
  Compute the Donor electronic potential contribution to the effective potentials

  beta - the inverse temperature Boltzmann factor in atomic units
*/

  if(ham_dia_mem_status==0){ cout<<"Error in donor_electronic_potential(): the diabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(ndia!=2){ cout<<"Error in donor_electronic_potential(): implementation only for ndia=2\n"; exit(0); }

  double en = 0.0;

  if(children.size()==0){
    double V0 = (this->ham_dia->get(0,0)).real(); 
    en = V0;
  }
  else{ cout<<"Error in donor_electronic_potential() not implemented for quantum nuclei\n"; exit(0); }

  return en;
}


double nHamiltonian::acceptor_electronic_potential(double beta){
/**
  Compute the Acceptor electronic potential contribution to the effective potentials

  beta - the inverse temperature Boltzmann factor in atomic units
*/

  if(ham_dia_mem_status==0){ cout<<"Error in acceptor_electronic_potential(): the diabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(ndia!=2){ cout<<"Error in acceptor_electronic_potential(): implementation only for ndia=2\n"; exit(0); }

  double en = 0.0;

  if(children.size()==0){
    double V1 = (this->ham_dia->get(1,1)).real(); 
    en = V1;
  }
  else{ cout<<"Error in acceptor_electronic_potential() not implemented for quantum nuclei\n"; exit(0); }

  return en;
}


double nHamiltonian::kinked_pair_electronic_potential(double beta){
/**
  Compute the Kinked Pair electronic potential contribution to the effective potentials

  beta - the inverse temperature Boltzmann factor in atomic units
*/

  if(ham_dia_mem_status==0){ cout<<"Error in kinked_pair_electronic_potential(): the diabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(ndia!=2){ cout<<"Error in kinked_pair_electronic_potential(): implementation only for ndia=2\n"; exit(0); }

  double en = 0.0;

  if(children.size()==0){
    double V0 = (this->ham_dia->get(0,0)).real(); 
    double V1 = (this->ham_dia->get(1,1)).real(); 
    double K = abs(this->ham_dia->get(0,1)); 
    if(beta * K > 1e-3){
      double Vg = 0.5 * (V0 + V1) - 0.5 * sqrt((V0 - V1) * (V0 - V1) + 4 * K * K);
      double Ve = 0.5 * (V0 + V1) + 0.5 * sqrt((V0 - V1) * (V0 - V1) + 4 * K * K);
      en = Vg - log(1 + exp(-beta * (Ve - Vg)) - exp(-beta * (V0 - Vg)) - exp(-beta * (V1 - Vg))) / beta;
    }
    else if(beta * abs(V0 - V1) > 1e-7){
      en = 0.5 * (V0 + V1) - log(pow(beta * K, 2) * sinh(0.5 * beta * (V0 - V1)) / (0.5 * beta * (V0 - V1))) / beta;
    }
    else{
      en = 0.5 * (V0 + V1) - log(pow(beta * K, 2)) / beta;
    }
  }
  else{ cout<<"Error in kinked_pair_electronic_potential() not implemented for quantum nuclei\n"; exit(0); }

  return en;
}


double nHamiltonian::kinetic_constraint_potential(double beta, double eta, double a){
/**
  Additional potential constraint added to the kinked-pair electronic contributions as defined in original KC-RPMD paper

  beta - the inverse temperature Boltzmann factor in atomic units
  eta - geometric parameter conserving free energy of kinked pair formation ad defined in second KC-RPMD paper
  a - is the gaussian functional limit parameter
*/
  double en = 0.0;

  if(children.size()==0){
    double V0 = (this->ham_dia->get(0,0)).real(); 
    double V1 = (this->ham_dia->get(1,1)).real(); 
    double K = abs(this->ham_dia->get(0,1)); 
    double w = (V0 - V1) / K;
    en = (a * w * w - log(sqrt(a / 3.1415) * eta)) / beta;
  } 
  else{ cout<<"Error in kinetic_constraint_potential() not implemented for quantum nuclei\n"; exit(0); }

  return en;
}


double nHamiltonian::heavy_side_potential(vector<double>& auxiliary_y, int theta, double beta, double b){
/**
  Heavy side coarse graining potential Vr as defined in original KC-RPMD paper

  auxiliary_y - is the classical electronic coordinate as defined in KC-RPMD
  theta - the coarse graining box center
  beta - the inverse temperature Boltzmann factor in atomic units
  b - is the heavyside functional limit parameter
*/
  double en = 0.0;

  if(abs(auxiliary_y[0] - theta) < 0.5){
    en = -log(1 / (1 + exp(b * (2 * abs(auxiliary_y[0] - theta) - 1)))) / beta;
  }
  else{
    en = (b * (2 * abs(auxiliary_y[0] - theta) - 1) - log(1 / (1 + exp(-b * (2 * abs(auxiliary_y[0] - theta) - 1))))) / beta;
  }

  return en;
}

double nHamiltonian::kcrpmd_effective_potential(vector<double>& auxiliary_y, const MATRIX& q, const MATRIX& invM, double beta, double eta, double a, double b){
/**
  Compute the KC-RPMD effective potential energy

  auxiliary_y - is the classical electronic coordinate as defined in KC-RPMD
  q - is a ndof x ntraj matrix of coordinates
  invM - is a ndof x 1 matrix of inverse masses of all DOFs
  beta - the inverse temperature Boltzmann factor in atomic units
  eta - geometric parameter conserving free energy of kinked pair formation ad defined in second KC-RPMD paper
  a - is the kinetic constraint ad-hoc parameter
  b - is the heavyside functional limit parameter
*/
  double en = 0.0;

  double V0 = donor_electronic_potential(beta) + heavy_side_potential(auxiliary_y,-1,beta,b); 
  double VKP = kinked_pair_electronic_potential(beta) + heavy_side_potential(auxiliary_y,0,beta,b) + kinetic_constraint_potential(beta,eta,a); 
  double V1 = acceptor_electronic_potential(beta) + heavy_side_potential(auxiliary_y,1,beta,b); 
  double Vshift = min({V0, VKP, V1});

  en = internal_ring_polymer_potential(q,invM,beta) + Vshift - log(exp(-beta * (V0 - Vshift)) + exp(-beta * (VKP - Vshift)) + exp(-beta * (V1 - Vshift))) / beta;

  return en;
}





}// namespace libnhamiltonian
}// liblibra

