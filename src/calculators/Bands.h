#ifndef BANDS_H
#define BANDS_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libcalculators{


void order_bands(int Norb, MATRIX* E, vector< pair<int,double> >& bands);

void populate_bands(int Nocc, int Norb, int degen, double Nel, int pop_opt, double kT, double etol,
         vector< pair<int,double> >& bands,vector< pair<int,double> >& occ);

void show_bands(int Norb, int Nocc, vector< pair<int,double> >& bands,vector< pair<int,double> >& occ);


}// namespace libcalculators

#endif // BANDS_H
