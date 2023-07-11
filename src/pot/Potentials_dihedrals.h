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

#ifndef POTENTIALS_DIHEDRALS_H
#define POTENTIALS_DIHEDRALS_H

#include "../math_linalg/liblinalg.h"
#include "../math_specialfunctions/libspecialfunctions.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;


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
}// liblibra

#endif //POTENTIALS_DIHEDRALS_H
