/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Potentials_elec.h"

using namespace libmmath;
using namespace libmmath::liblinalg;


namespace libpot{

double Elec_Coulomb(VECTOR& ri,VECTOR& rj,     /*Inputs*/             
                    VECTOR& fi,VECTOR& fj,     /*Outputs*/
                    double qi,double qj,
                    double eps,double delta){  /*Parameters*/
//****************** double Lennard-Jones potential **************************
//*                                                                          *
//*       E = electric*qi*qj/(eps*|rij+delta|)                               *
//*                                                                          *
//****************************************************************************
  double energy,r2,r6,r12,d1,d2;
  VECTOR rij = ri - rj;
  d1 = rij.length();
  d2 = d1 + delta;
  energy = electric*(qi*qj/(eps*d2));
  fi = (energy/d2)*(rij/d1);
  fj = -fi;
  return energy;
}


}// namespace libpot


