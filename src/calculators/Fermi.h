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

#ifndef FERMI_H
#define FERMI_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libcalculators{


double fermi_population(double e,double ef,double degen, double kT);

double fermi_integral(vector<double>& bnds, double ef, double degen, double kT);
double fermi_integral(boost::python::list bnds, double ef, double degen, double kT);

double fermi_energy(vector<double>& bnds,double Nel,double degen, double kT, double etol);
double fermi_energy(boost::python::list bnds,double Nel,double degen, double kT, double etol);



}// namespace libcalculators

#endif // FERMI_H
