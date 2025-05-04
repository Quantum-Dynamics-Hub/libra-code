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
  \file nHamiltonian_compute_RPMD.cpp
  \brief The file implements the calculations general to RPMD (Ring Polymer Molecular Dynamics)
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


double RPMD_internal_potential(const MATRIX& q, const MATRIX& invM, double beta){
/**
  Compute the ring polymer internal potential Uint

  q - is a ndof x ntraj matrix of coordinates
  invM - is a ndof x 1 matrix of inverse masses of all DOFs
  beta - the inverse temperature Boltzmann factor in atomic units
*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;

  int dof, traj;

  double en;

  if(ntraj==1){
    en = 0.0;
  }
  else{ cout<<"Error in RPMD_internal_potential(), not implemented for quantum nuclei\n"; exit(0); }

  return en;
}


MATRIX RPMD_internal_force(const MATRIX& q, const MATRIX& invM, double beta){
/**
  Compute the ring polymer internal potential Uint

  q - is a ndof x ntraj matrix of coordinates
  invM - is a ndof x 1 matrix of inverse masses of all DOFs
  beta - the inverse temperature Boltzmann factor in atomic units

  Returns:
  f - is a ndof x ntraj matrix that will contain the forces due to Uint

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;

  int dof, traj;

  MATRIX f(ndof, ntraj);

  if(ntraj==1){
    for(dof=0; dof<ndof; dof++){   
      f.set(dof, 0, 0.0);
    }// for dof
  }
  else{ cout<<"Error in RPMD_internal_force(), not implemented for quantum nuclei\n"; exit(0); }

  return f;
}


}// namespace libnhamiltonian
}// liblibra

