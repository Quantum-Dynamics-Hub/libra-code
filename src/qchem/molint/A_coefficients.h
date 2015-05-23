#ifndef A_COEFFICIENTS_H
#define A_COEFFICIENTS_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{

// Auxiliary  functions for computation of STO integrals

double A_coefficient_general(int u,int v, int na, int nb, int la, int lb, int m);
void generate_coefficients();
void A_coefficients_fast(int u,int v, int na, int nb, int la, int lb, int m,double** A);

void Aux_F1(double rhoA, double rhoB, double* f,int mu_max);


}// namespace libmolint
}// namespace libqchem

#endif // A_COEFFICIENTS_H
