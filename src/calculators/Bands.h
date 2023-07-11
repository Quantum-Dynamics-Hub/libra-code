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
  \file Bands.h
  \brief The file describes functions for ordering, converting, and printing bands (energies and occupations) information
    
*/

#ifndef BANDS_H
#define BANDS_H


#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libcalculators namespace
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
}// liblibra

#endif // BANDS_H
