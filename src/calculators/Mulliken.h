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

#ifndef MULLIKEN_H
#define MULLIKEN_H

#include "../mmath/libmmath.h"
using namespace libmmath;


namespace libcalculators{


void update_Mull_orb_pop(MATRIX*, MATRIX*, vector<double>&, vector<double>&);

void update_Mull_charges(vector<int>&, vector<int>&, vector<vector<int> >&,vector<double>&,
                         vector<double>&, vector<double>&, vector<double>&, vector<double>&);

void update_Mull_charges(vector<int>& ao_to_atom_map, vector<double>& Zeff,
                         vector<double>& Mull_orb_pop_gross, vector<double>& Mull_orb_pop_net,
                         vector<double>& Mull_charges_gross, vector<double>& Mull_charges_net);



}// namespace libcalculators

#endif // MULLIKEN_H

