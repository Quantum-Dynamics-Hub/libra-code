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
#include "../../nhamiltonian/libnhamiltonian.h"
#include "../dyn_control_params.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libnhamiltonian;
namespace bp = boost::python;


/// libdyn namespace
namespace libdyn{

/// libqtag namespace
namespace libqtag{


///=============== (qtag.cpp) ===================

/// Wavefunction
CMATRIX qtag_psi(MATRIX& q, MATRIX& q1, MATRIX& p1, MATRIX& alp1, MATRIX& s1, CMATRIX& Coeff);

/// Elementary overlap matrix
CMATRIX qtag_overlap_elementary(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s);

/// Elementary kinetic matrix
CMATRIX qtag_kinetic_elementary(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s, MATRIX& invM);

/// Global overlap matrix
CMATRIX qtag_overlap(vector<int>& active_states, CMATRIX& ovlp, int nstates);

/// Bra-Ket Averaged Taylor expnasion Approximation
complex<double> BAT(CMATRIX* Ham1, CMATRIX* Ham2, vector<CMATRIX*>& dHam1, vector<CMATRIX*>& dHam2,
                    MATRIX& q1, MATRIX& p1, MATRIX& s1, MATRIX& alp1, int n1, 
                    MATRIX& q2, MATRIX& p2, MATRIX& s2, MATRIX& alp2, int n2);

complex<double> BATe(int i, int j,
                     MATRIX& q1, MATRIX& p1, MATRIX& s1, MATRIX& alp1, int n1,
                     MATRIX& q2, MATRIX& p2, MATRIX& s2, MATRIX& alp2, int n2,
                     double AA, double BB, double CC, nHamiltonian& ham);

/// Local Harmonic approximation to Hamiltonian
complex<double> LHA(CMATRIX* Ham1, CMATRIX* Ham2, 
                    vector<CMATRIX*>& dHam1, vector<CMATRIX*>& dHam2,
                    vector<CMATRIX*>& d2Ham1, vector<CMATRIX*>& d2Ham2,
                    MATRIX& q1, MATRIX& p1, MATRIX& s1, MATRIX& alp1, int n1, 
                    MATRIX& q2, MATRIX& p2, MATRIX& s2, MATRIX& alp2, int n2);

complex<double> LHAe(int i, int j, 
                     MATRIX& q1, MATRIX& p1, MATRIX& s1, MATRIX& alp1, int n1,
                     MATRIX& q2, MATRIX& p2, MATRIX& s2, MATRIX& alp2, int n2,
                     double AA, double BB, double CC, nHamiltonian& ham);

/// Elementary potential & coupling matrix
CMATRIX qtag_potential(MATRIX& q1, MATRIX& p1, MATRIX& s1, MATRIX& alp1, int n1, vector<int>& traj_on_surf_n1,
                       MATRIX& q2, MATRIX& p2, MATRIX& s2, MATRIX& alp2, int n2, vector<int>& traj_on_surf_n2,
                       nHamiltonian& ham, int method, double AA, double BB, double CC);

/// super-Hamiltonian and super-Overlap for all trajectories
void qtag_hamiltonian_and_overlap(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff,
                                  vector<int>& active_states, MATRIX& invM, 
                                  nHamiltonian& ham, bp::object compute_ham_funct, bp::dict compute_ham_params,
                                  bp::dict& dyn_params,
                                  CMATRIX& super_ovlp, CMATRIX& super_ham);

/// QTAG momentum
CMATRIX qtag_momentum(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff);




}// namespace libqtag
}// namespace libdyn
}// liblibra

#endif  // QTAG_H
