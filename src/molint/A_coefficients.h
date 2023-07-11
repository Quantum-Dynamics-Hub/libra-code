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

#ifndef A_COEFFICIENTS_H
#define A_COEFFICIENTS_H



/// liblibra namespace
namespace liblibra{

namespace libmolint{

// Auxiliary  functions for computation of STO integrals

double A_coefficient_general(int u,int v, int na, int nb, int la, int lb, int m);
void generate_coefficients();
void A_coefficients_fast(int u,int v, int na, int nb, int la, int lb, int m,double** A);

void Aux_F1(double rhoA, double rhoB, double* f,int mu_max);


}// namespace libmolint
}// namespace liblibra

#endif // A_COEFFICIENTS_H
