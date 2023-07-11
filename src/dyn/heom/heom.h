/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file heom.h
  \brief The file describes the functions used to setup and perform the HEOM calculations

  Based on the Fortran 90 code of Amber Jain & Joe Subotnik:  https://github.com/subotnikgroup/HEOM_Amber

*/

#ifndef HEOM_H
#define HEOM_H


#include "../../util/libutil.h"
#include "../../math_linalg/liblinalg.h"
//#include "../../math_random/librandom.h"
//#include "../../math_specialfunctions/libspecialfunctions.h"
//#include "../../hamiltonian/libhamiltonian.h"
#include "../../Units.h"


namespace liblibra{

using namespace libutil;
using namespace liblinalg;
//using namespace librandom;
//using namespace libspecialfunctions;
//using namespace libhamiltonian;

/// libdyn namespace
namespace libdyn{

/// libheom namespace
namespace libheom{

///=============== In the heom.cpp ====================

/// General hierarchy setup
int compute_nn_tot(int d, int max_tier);

vector< vector<int> > gen_next_level(vector< vector<int> >& parents);

void gen_hierarchy(int d, int max_tier, int verbosity,
                   vector< vector<int> >& all_vectors, 
                   vector< vector<int> >& vec_plus, 
                   vector< vector<int> >& vec_minus);


/// General calculations
vector<int> filter(vector<CMATRIX>& rho, vector<int>& adm_list, double tolerance, int do_zeroing);

CMATRIX compute_deriv_n(int n, vector<CMATRIX>& rho, CMATRIX& Ham, vector<CMATRIX>& el_phon_coupl,
        double eta, double temperature,
        vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara,
        int truncation_scheme, complex<double> truncation_prefactor, int do_scale, vector<int>& nonzero,
        vector< vector<int> >& nvectors, vector< vector<int> >& vec_plus, vector< vector<int> >& vec_minus        
        );
CMATRIX compute_heom_derivatives(CMATRIX& RHO, bp::dict prms);


/// Bath setup
vector<CMATRIX> initialize_el_phonon_couplings(int nquant);
complex<double> compute_matsubara_sum(vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara, int KK);
void setup_bath(int KK, double eta, double gamma, double temperature,
                vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara);


/// Auxiliary functions
void unpack_mtx(vector<CMATRIX>& rho_unpacked, CMATRIX& RHO);
void pack_mtx(vector<CMATRIX>& rho_unpacked, CMATRIX& RHO);
void scale_rho(vector<CMATRIX>& rho, vector<CMATRIX>& rho_scaled, bp::dict prms);
void inv_scale_rho(vector<CMATRIX>& rho, vector<CMATRIX>& rho_scaled, bp::dict prms);






}// namespace libheom
}// namespace libdyn
}// liblibra

#endif  // HEOM_H
