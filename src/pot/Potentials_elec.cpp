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
  \file Potentials_elec.cpp
  This file implements the electrostatic interaction potentials
    
*/


#include "Potentials_elec.h"

/// liblibra namespace
namespace liblibra{


namespace libpot{

double Elec_Coulomb(VECTOR& ri,VECTOR& rj,     /*Inputs*/             
                    VECTOR& fi,VECTOR& fj,     /*Outputs*/
                    double qi,double qj,
                    double eps,double delta){  /*Parameters*/
/**
  This is a simple Coulombic interaction potential:

  E = qi*qj/(eps*|rij+delta|) 

  E - the energy in Ha *(atomic energy units)
  qi, qj - the charges in electron charge units
  rij - distance, in atomic units (Bohr)
  delta - a shift (e.g. to avoid charge clashes) in Bohr
  eps - the electric permittivity in units in which the permittivity of vacuum is = 1

  \param[in]  ri Coordinates of the particle i
  \param[in]  rj Coordinates of the particle j
  \param[out] fi Force acting on particle i
  \param[out] fj Force acting on particle j
  \param[in] qi Charge of the particle i
  \param[in] qj Charge of the particle j
  \param[in] e Charge of the particle i
  \param[in] eps Dielectric permittivity ( = 1 for vacuum)
  \param[in] delta Shift in the position

*/

  double energy,r2,r6,r12,d1,d2;
  VECTOR rij = ri - rj;
  d1 = rij.length();
  d2 = d1 + delta;
//  energy = electric*(qi*qj/(eps*d2));
  energy = (qi*qj/(eps*d2));
  fi = (energy/d2)*(rij/d1);
  fj = -fi;
  return energy;
}


}// namespace libpot
}// liblibra

