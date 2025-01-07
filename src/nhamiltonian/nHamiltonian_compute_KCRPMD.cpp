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


vector<CMATRIX> nHamiltonian::generate_m_matrices(double beta){
/**
  Generate set of M matrices for each trajectory

  beta - the inverse temperature Boltzmann factor in atomic units
*/

  vector<CMATRIX> res; res = vector<CMATRIX>(children.size(), CMATRIX(ndia,ndia));

  for(int traj=0; traj<children.size(); traj++){
    //std::cout << children[traj]->ham_dia->get(0,1) << std::endl;
    res[traj] = *(children[traj]->ham_dia);

  }// traj

  return res;

}


double nHamiltonian::kcrpmd_effective_potential(vector<double>& auxiliary_y){
/**
  Compute the KC-RPMD effective potential energy

  auxiliary_y - is the classical electronic coordinate as defined in KC-RPMD

*/

  double en = 0.0;

  en = 0.1 * auxiliary_y[0];

  return en;


}





}// namespace libnhamiltonian
}// liblibra

