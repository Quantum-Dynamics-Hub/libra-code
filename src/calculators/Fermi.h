#ifndef FERMI_H
#define FERMI_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libcalculators{


double fermi_population(double e,double ef,double degen, double kT);
double fermi_integral(vector<double>& bnds, double ef, double degen, double kT);
double fermi_energy(vector<double>& bnds,double Nel,double degen, double kT, double etol);


}// namespace libcalculators

#endif // FERMI_H
