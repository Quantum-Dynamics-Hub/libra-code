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
  \file Surface_Hopping.h
  \brief The file describes the functions used in surface hopping methods
    
*/

#ifndef SURFACE_HOPPING_H
#define SURFACE_HOPPING_H


// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../nhamiltonian/libnhamiltonian.h"

// Dynamics classes
#include "nuclear/libnuclear.h"
#include "electronic/libelectronic.h"
//#include "ensemble/libensemble.h"

/// liblibra namespace
namespace liblibra{

using namespace libio;
using namespace libnhamiltonian;
namespace bp = boost::python;


/// libdyn namespace
namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
//using namespace libensemble;

///================  In tsh_prob_lz.cpp  ===================================

//MATRIX compute_hopping_probabilities_lz(nHamiltonian* ham, int rep, MATRIX& p, const MATRIX& invM, MATRIX& prev_ham_dia);
//MATRIX compute_hopping_probabilities_lz(nHamiltonian& ham, int rep, MATRIX& p, const MATRIX& invM, MATRIX& prev_ham_dia);



///================  In tsh_hungarian.cpp  =================================
vector<int> get_permutation(vector<vector<int> >& inp);
vector<int> Munkres_Kuhn_minimize(MATRIX& _X, int verbosity);
vector<int> Munkres_Kuhn_maximize(MATRIX& _X, int verbosity);


}// namespace libdyn
}// liblibra

#endif // SURFACE_HOPPING_H
