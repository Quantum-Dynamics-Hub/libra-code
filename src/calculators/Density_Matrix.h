#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libcalculators{

void compute_density_matrix(vector< pair<int,double> >& occ, MATRIX* C, MATRIX* P);

void Fock_to_P(int Norb,int Nocc, int degen, double Nel, std::string eigen_method, int pop_opt, double kT, double etol,
               MATRIX* Fao, MATRIX* Sao, MATRIX* C, MATRIX* E,
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,
               MATRIX* P, int BM, vector<Timer>& bench_t);


}// namespace libcalculators

#endif // DENSITY_MATRIX_H
