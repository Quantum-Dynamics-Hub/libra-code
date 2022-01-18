/*********************************************************************************
* Copyright (C) 2022 Matthew Dutra, Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file qtag.h
  \brief This file defines the functions need for QTAG method

*/

#ifndef QTAG_H
#define QTAG_H


#include "../../math_linalg/liblinalg.h"
#include "../../hamiltonian/libhamiltonian.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libhamiltonian;
namespace bp = boost::python;


/// libdyn namespace
namespace libdyn{

/// libqtag namespace
namespace libqtag{


///=============== (qtag.cpp) ===================

/// Wavefunction
CMATRIX qtag_psi(MATRIX q, MATRIX& q1, MATRIX& p1, MATRIX& alp1, MATRIX& s1, CMATRIX& Coeff);

/// Elementary overlap matrix
CMATRIX qtag_overlap_elementary(MATRIX& q1, MATRIX& p1, MATRIX& alp1, MATRIX& s1,
                                MATRIX& q2, MATRIX& p2, MATRIX& alp2, MATRIX& s2);

/// Elementary kinetic matrix
CMATRIX qtag_kinetic_elementary(MATRIX q, MATRIX& p, MATRIX& alp, MATRIX& s, MATRIX& invM);

/// Global overlap matrix
CMATRIX qtag_overlap(vector<int>& active_states, CMATRIX& ovlp, int nstates);

///Hamiltonian for all trajectories
CMATRIX qtag_hamiltonian(MATRIX q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff,
                         vector<int>& active_states, CMATRIX& ovlp, CMATRIX& kin,
                         MATRIX& invM, nHamiltonian& ham, bp::object compute_ham_funct,  
                         bp::dict& compute_ham_params);
/// QTAG momentum
CMATRIX qtag_momentum(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff);




}// namespace libqtag
}// namespace libdyn
}// liblibra

#endif  // QTAG_H
