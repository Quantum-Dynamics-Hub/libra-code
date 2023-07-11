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
/**
  \file Fermi.h
  \brief The file defines functions for Fermi energy/population calculations
    
*/

#ifndef FERMI_H
#define FERMI_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <vector>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace std;
using namespace boost::python;
using namespace liblinalg;

/// libcalculators namespace
namespace libcalculators{

double fermi_population(double e,double ef,double degen, double kT);  

double fermi_integral(std::vector<double>& bnds, double ef, double degen, double kT);
double fermi_integral(boost::python::list bnds, double ef, double degen, double kT);

double fermi_energy(std::vector<double>& bnds,double Nel,double degen, double kT, double etol);
double fermi_energy(boost::python::list bnds,double Nel,double degen, double kT, double etol); 


/// For FOE: Fermi operator expansion
double p_up(double e, double e_up, double de);
double p_dn(double e, double e_dn, double de);
double p_ef(double e, double ef, double de);
void Chebyshev_coeff(vector<double>& C, double (*f)(double x, double y, double z), double ef, double de, int N);
double Chebyshev_fit(MATRIX& H, MATRIX& P, double (*f)(double _x, double _y, double _z), double ef, double de, int np);


}// namespace libcalculators

}// liblibra

#endif // FERMI_H
