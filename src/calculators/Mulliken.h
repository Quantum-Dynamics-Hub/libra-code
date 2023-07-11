/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

/**
  \file Mulliken.h
  \brief The file defines functions for Mulliken population analysis    
    
*/

#ifndef MULLIKEN_H
#define MULLIKEN_H

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libcalculators namespace
namespace libcalculators{


void update_Mull_orb_pop(MATRIX*, MATRIX*, vector<double>&, vector<double>&);

void update_Mull_charges(vector<int>&, vector<int>&, vector<vector<int> >&,vector<double>&,
                         vector<double>&, vector<double>&, vector<double>&, vector<double>&);

void update_Mull_charges(vector<int>& ao_to_atom_map, vector<double>& Zeff,
                         vector<double>& Mull_orb_pop_gross, vector<double>& Mull_orb_pop_net,
                         vector<double>& Mull_charges_gross, vector<double>& Mull_charges_net);

boost::python::list update_Mull_orb_pop(MATRIX P, MATRIX S);

boost::python::list update_Mull_charges(
 vector<int>& ao_to_atom_map, vector<double>& Zeff,
 vector<double>& Mull_orb_pop_gross, vector<double>& Mull_orb_pop_net
);


}// namespace libcalculators
}// liblibra

#endif // MULLIKEN_H

