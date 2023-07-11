/*********************************************************************************
* Copyright (C) 2016-2020 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
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


#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libdyn namespace
namespace libdyn{

/// libgwp namespace
namespace libgwp{

///=============== gwp.cpp ===================

int check_dimensions(std::string function_name, MATRIX& R1, MATRIX& P1, MATRIX& R2, MATRIX& P2);

complex<double> gwp_value(MATRIX& r, MATRIX& R, MATRIX& P, double gamma,  double alp, double hbar);


double gwp_product_decomposition(double q1, double p1, double gamma1, double alp1,
                                 double q2, double p2, double gamma2, double alp2,
                                 double& q, double& p, double& gamma, double& alp);

double gwp_product_decomposition(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                                 MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2,
                                 MATRIX& q,  MATRIX& p,  MATRIX& gamma,  MATRIX& alp);


///=============== Overlaps (gwp_overlap.cpp) ===================

complex<double> gwp_overlap(double q1, double p1, double gamma1, double alp1,
                            double q2, double p2, double gamma2, double alp2);

complex<double> gwp_overlap(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                            MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2);

CMATRIX gwp_overlap_matrix(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                           MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2);

CMATRIX gwp_overlap_matrix(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1, vector<int>& state1,
                           MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2, vector<int>& state2);

complex<double> gwp_overlap(MATRIX& R1, MATRIX& P1, double gamma1, 
                            MATRIX& R2, MATRIX& P2, double gamma2, 
                            double alp, double hbar);

complex<double> gwp_overlap(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                            MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2,
                            MATRIX& q3, MATRIX& p3, MATRIX& gamma3, MATRIX& alp3);


///=============== Transition dipole moments (gwp_dipole.cpp) ===================

complex<double> gwp_dipole(double q1, double p1, double gamma1, double alp1,
                           double q2, double p2, double gamma2, double alp2);

CMATRIX gwp_dipole(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                   MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2);

CMATRIX gwp_dipole(MATRIX& R1, MATRIX& P1, double gamma1, 
                   MATRIX& R2, MATRIX& P2, double gamma2, 
                   double alp, double hbar);


///=============== Derivative coupling (gwp_coupling.cpp) ===================

complex<double> gwp_coupling(double q1, double p1, double gamma1, double alp1,
                             double q2, double p2, double gamma2, double alp2);

CMATRIX gwp_coupling(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                     MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2);

CMATRIX gwp_coupling(MATRIX& R1, MATRIX& P1, double gamma1, 
                     MATRIX& R2, MATRIX& P2, double gamma2, 
                     double alp, double hbar);

///=============== Kinetic energy operator (gwp_kinetic.cpp) ===================

complex<double> gwp_kinetic(double q1, double p1, double gamma1, double alp1,
                            double q2, double p2, double gamma2, double alp2);

complex<double> gwp_kinetic(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                            MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2);

complex<double> gwp_kinetic(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                            MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2,
                            MATRIX& iM);

CMATRIX gwp_kinetic_matrix(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                           MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2, 
                           MATRIX& invM);

complex<double> gwp_kinetic(MATRIX& R1, MATRIX& P1, double gamma1, 
                            MATRIX& R2, MATRIX& P2, double gamma2, 
                            double alp, double hbar);




}// namespace libgwp
}// namespace libdyn
}// liblibra

#endif  // GWP_H
