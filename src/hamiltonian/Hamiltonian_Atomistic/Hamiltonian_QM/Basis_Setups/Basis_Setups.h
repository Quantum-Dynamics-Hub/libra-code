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

#ifndef BASIS_SETUPS_H
#define BASIS_SETUPS_H

#include <sstream>
using namespace std;

#include "../../../../mmath/libmmath.h"
using namespace libmmath;

#include "../../../../qchem/qobjects/libqobjects.h"
using namespace libqchem::libqobjects;

#include "../../../../qchem/basis/libbasis.h"
using namespace libqchem::libbasis;

#include "../Model_Parameters/libmodel_parameters.h"




/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_atomistic namespace
namespace libhamiltonian_atomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{

using namespace libmodel_parameters;

namespace libbasis_setups{




// Basis_STO_3G_DZ.cpp

void set_basis_STO_3G_DZ(vector<std::string>& at_type, vector<VECTOR>& R,  Model_Parameters& modpar, int verb,
  vector<AO>& basis_ao, int& Nelec, int& Norb, 
  vector<vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map);

boost::python::list set_basis_STO_3G_DZ(vector<std::string>& at_type, vector<VECTOR>& R,  Model_Parameters& modpar, int verb);



}// namespace libbasis_setups
}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian



#endif // BASIS_SETUPS_H

