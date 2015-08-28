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

#ifndef DENSITY_OF_STATES_H
#define DENSITY_OF_STATES_H

#include "MOAO.h"
#include "Control_Parameters.h"
#include "Model_Parameters.h"
#include "Electronic.h"
#include "Nuclear.h"
#include "Memory.h"


void compute_dos(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                 vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                 Electronic* el, Memory* mem);


#endif // DENSITY_OF_STATES_H
