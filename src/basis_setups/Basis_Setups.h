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

#ifndef BASIS_SETUPS_H
#define BASIS_SETUPS_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <sstream>
#endif

#include "../math_linalg/liblinalg.h"
#include "../qobjects/libqobjects.h"
#include "../basis/libbasis.h"
#include "../model_parameters/libmodel_parameters.h"

/// liblibra namespace
namespace liblibra{

using namespace std;
using namespace liblinalg;
using namespace libqobjects;
using namespace libbasis;
using namespace libmodel_parameters;


namespace libbasis_setups{




// Basis_STO_3G_DZ.cpp

void set_basis_STO_3G_DZ(vector<std::string>& at_type, vector<VECTOR>& R,  Model_Parameters& modpar, int verb,
  vector<AO>& basis_ao, int& Nelec, int& Norb, 
  vector<vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map);

boost::python::list set_basis_STO_3G_DZ(vector<std::string>& at_type, vector<VECTOR>& R,  Model_Parameters& modpar, int verb);



}// namespace libbasis_setups
}// namespace liblibra



#endif // BASIS_SETUPS_H

