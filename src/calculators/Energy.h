#ifndef ENERGY_H
#define ENERGY_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libcalculators{

double energy_elec(int Norb,MATRIX* Pao,MATRIX* Hao,MATRIX* Fao);
double energy_elec(int Norb, 
                   MATRIX* P_alp, MATRIX* P_bet, 
                   MATRIX* Hao_alp, MATRIX* Hao_bet,
                   MATRIX* Fao_alp, MATRIX* Fao_bet,
                   MATRIX* dFao_alp_dP_alp, MATRIX* dFao_alp_dP_bet,
                   MATRIX* dFao_bet_dP_alp, MATRIX* dFao_bet_dP_bet,
                   MATRIX* temp
                  );


}// namespace libcalculators

#endif // ENERGY_H
