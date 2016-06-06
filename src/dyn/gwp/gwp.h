/*********************************************************************************
* Copyright (C) 2016 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file gwp.h
  \brief This file defines the functions need for Gaussian Wave Packet propagation - including 
  some auxiliary integrals    
*/

#ifndef GWP_H
#define GWP_H


#include "../../mmath/libmmath.h"
using namespace libmmath;

/// libdyn namespace
namespace libdyn{

/// libgwp namespace
namespace libgwp{


complex<double> gwp_overlap(MATRIX& R1, MATRIX& P1, double gamma1, 
                            MATRIX& R2, MATRIX& P2, double gamma2, 
                            double alp, double hbar);

CMATRIX gwp_coupling(MATRIX& R1, MATRIX& P1, double gamma1, 
                     MATRIX& R2, MATRIX& P2, double gamma2, 
                     double alp, double hbar);

complex<double> gwp_kinetic(MATRIX& R1, MATRIX& P1, double gamma1, 
                            MATRIX& R2, MATRIX& P2, double gamma2, 
                            double alp, double hbar);




}// namespace libgwp
}// namespace libdyn

#endif  // GWP_H
