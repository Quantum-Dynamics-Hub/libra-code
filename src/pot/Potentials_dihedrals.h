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

#ifndef POTENTIALS_DIHEDRALS_H
#define POTENTIALS_DIHEDRALS_H

#include "../mmath/libmmath.h"
using namespace libmmath;
using namespace libmmath::liblinalg;


namespace libpot{




//------------------ Dihedral/Torsion potentials ------------------------------

double Dihedral_General(VECTOR& ri,VECTOR& rj,VECTOR& rk,VECTOR& rl, /*Inputs*/
                        VECTOR& fi,VECTOR& fj,VECTOR& fk,VECTOR& fl, /*Outputs*/
                        double Vphi,double phi0,int n,int opt        /*Parameters*/
                        );
double Dihedral_Fourier(VECTOR& ri,VECTOR& rj,VECTOR& rk,VECTOR& rl,    /*Inputs*/
                        VECTOR& fi,VECTOR& fj,VECTOR& fk,VECTOR& fl,    /*Outputs*/
                        double Vphi1,double Vphi2,double Vphi3,int opt  /*Parameters*/
                        );

}// namespace libpot

#endif //POTENTIALS_DIHEDRALS_H
