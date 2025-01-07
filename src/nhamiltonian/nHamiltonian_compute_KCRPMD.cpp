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


double nHamiltonian::kcrpmd_effective_potential(vector<double>& auxiliary_y, double a, double b){
/**
  Compute the KC-RPMD effective potential energy

  auxiliary_y - is the classical electronic coordinate as defined in KC-RPMD
  a - is the kinetic constraint ad-hoc parameter
  b - is the heavyside functional limit parameter
*/
  double en = 0.0;
  en = 0.1 * auxiliary_y[0];
  return en;
}





}// namespace libnhamiltonian
}// liblibra

