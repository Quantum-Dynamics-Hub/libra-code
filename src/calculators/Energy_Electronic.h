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

#ifndef ENERGY_ELECTRONIC_H
#define ENERGY_ELECTRONIC_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libcalculators{

double energy_elec(MATRIX* Pao,MATRIX* Hao,MATRIX* Fao);
double energy_elec(MATRIX Pao,MATRIX Hao,MATRIX Fao);

double energy_elec(MATRIX* P_alp, MATRIX* P_bet, 
                   MATRIX* Hao_alp, MATRIX* Hao_bet,
                   MATRIX* Fao_alp, MATRIX* Fao_bet,
                   MATRIX* dFao_alp_dP_alp, MATRIX* dFao_alp_dP_bet,
                   MATRIX* dFao_bet_dP_alp, MATRIX* dFao_bet_dP_bet,
                   MATRIX* temp
                  );
double energy_elec(MATRIX P_alp, MATRIX P_bet, 
                   MATRIX Hao_alp, MATRIX Hao_bet,
                   MATRIX Fao_alp, MATRIX Fao_bet,
                   MATRIX dFao_alp_dP_alp, MATRIX dFao_alp_dP_bet,
                   MATRIX dFao_bet_dP_alp, MATRIX dFao_bet_dP_bet,
                   MATRIX temp
                  );




}// namespace libcalculators

#endif // ENERGY_ELECTRONIC_H
