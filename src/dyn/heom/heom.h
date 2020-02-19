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


#include "../../math_linalg/liblinalg.h"
//#include "../../math_random/librandom.h"
//#include "../../math_specialfunctions/libspecialfunctions.h"
//#include "../../hamiltonian/libhamiltonian.h"
#include "../../Units.h"


namespace liblibra{

using namespace liblinalg;
//using namespace librandom;
//using namespace libspecialfunctions;
//using namespace libhamiltonian;

/// libdyn namespace
namespace libdyn{

/// libheom namespace
namespace libheom{

///=============== In the heom.cpp ====================
vector< vector<int> > gen_next_level(vector<int>& parent);
vector< vector<int> > gen_next_level2(vector< vector<int> >& parents);
int is_equal(vector<int>& vec1, vector<int>& vec2);
int is_included(vector<int>& vec1, vector<vector<int> >& vec);
void gen_hierarchy(int d, int max_tier, int verbosity,
                   vector< vector<int> >& all_vectors, 
                   vector< vector<int> >& vec_plus, 
                   vector< vector<int> >& vec_minus);



int compute_nn_tot(int nquant, int KK, int LL);

vector<int> allocate_1D(int sz1);
vector< vector<int> > allocate_2D(int sz1, int sz2);
vector< vector< vector<int> > > allocate_3D(int sz1, int sz2, int sz3);


void compute_nn(int nquant, int KK, int LL, vector<int>& map_sum, 
                vector< vector< vector<int> > >& nn);

void compute_nn_sum_L(int nquant, int KK, int L, int& n_beg, int& n_end, 
                      vector< vector< vector<int> > >& nn);


vector< vector<int> > index_int2vec(vector< vector< vector<int> > >& nn, int n, int nquant, int KK);
int sum2D(vector< vector<int> >& nvec);
int index_vec2int(vector< vector< vector<int> > >& nn, vector< vector<int> >& nvec, int LL);
void compute_map(int nquant, int KK, int LL, 
                 vector< vector< vector<int> > >& nn,
                 vector< vector< vector<int> > >& map_nplus,
                 vector< vector< vector<int> > >& map_nneg);


vector<int> filter(vector<CMATRIX>& rho, double tolerance);


vector<CMATRIX> initialize_el_phonon_couplings(int nquant);

complex<double> compute_matsubara_sum(vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara, int KK);

CMATRIX compute_deriv_n(int n, vector<CMATRIX>& rho, CMATRIX& Ham, vector<CMATRIX>& projectors,
        double eta, double temperature,
        vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara,
        int do_truncate, int do_scale,
        vector< vector< vector<int> > >& nn, int KK, vector<int>& zero,
        vector< vector< vector<int> > >& map_nplus, vector< vector< vector<int> > >& map_nneg        
        );


void unpack_rho(vector<CMATRIX>& rho_unpacked, CMATRIX& RHO);
void pack_rho(vector<CMATRIX>& rho_unpacked, CMATRIX& RHO);
CMATRIX compute_heom_derivatives(CMATRIX& RHO, bp::dict prms);

void setup_bath(bp::dict params, vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara);





}// namespace libheom
}// namespace libdyn
}// liblibra

#endif  // HEOM_H
