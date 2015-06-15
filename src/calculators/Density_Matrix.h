#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libcalculators{

void compute_density_matrix(vector< pair<int,double> >& occ, MATRIX* C, MATRIX* P);
MATRIX compute_density_matrix(boost::python::list occ, MATRIX C);

void Fock_to_P(MATRIX* Fao, MATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               MATRIX* E, MATRIX* C, MATRIX* P,                                              /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,           /*Outputs*/
               int BM, vector<Timer>& bench_t);                                              /*Benchmarking data*/
void Fock_to_P(MATRIX* Fao, MATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               MATRIX* E, MATRIX* C, MATRIX* P,                                                         /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ                       /*Outputs*/
              );
boost::python::list Fock_to_P(MATRIX Fao, MATRIX Sao, double Nel, double degen, double kT, double etol, int pop_opt);




}// namespace libcalculators

#endif // DENSITY_MATRIX_H
