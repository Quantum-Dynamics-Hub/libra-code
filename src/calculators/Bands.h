#ifndef BANDS_H
#define BANDS_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libcalculators{


void convert_1(boost::python::list bands,  vector< pair<int,double> >& int_bands);
boost::python::list convert_2( vector< pair<int,double> >& bands);

void order_bands(MATRIX* E, vector< pair<int,double> >& bands);
boost::python::list order_bands(MATRIX E);

void populate_bands(double Nel, double degen, double kT, double etol, int pop_opt,
         vector< pair<int,double> >& bands,vector< pair<int,double> >& occ);
boost::python::list populate_bands(double Nel, double degen, double kT, double etol, int pop_opt,
         boost::python::list bands);


void show_bands(int Norb, int Nocc, vector< pair<int,double> >& bands,vector< pair<int,double> >& occ);


}// namespace libcalculators

#endif // BANDS_H
