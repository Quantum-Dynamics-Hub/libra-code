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
  \file libheom.cpp
  \brief The file implements Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libheom.h"

/// liblibra namespace
namespace liblibra{

using namespace boost::python;

namespace libdyn{
namespace libheom{




void export_heom_objects(){
/** 
  \brief Exporter of libheom classes and functions
*/
 
  vector< vector<int> > (*expt_gen_next_level_v1)
  (vector<int>& parent) = &gen_next_level;
  def("gen_next_level", expt_gen_next_level_v1);

  vector< vector<int> > (*expt_gen_next_level2_v1)
  (vector< vector<int> >& parents) = &gen_next_level2;
  def("gen_next_level2", expt_gen_next_level2_v1);

  int (*expt_is_equal_v1)
  (vector<int>& vec1, vector<int>& vec2) = &is_equal;
  def("is_equal", expt_is_equal_v1);

  int (*expt_is_included_v1)
  (vector<int>& vec1, vector<vector<int> >& vec) = is_included;
  def("is_included", expt_is_included_v1);

  void (*expt_gen_hierarchy_v1)
  (int d, int max_tier, int verbosity,
   vector< vector<int> >& all_vectors,  
   vector< vector<int> >& vec_plus, 
   vector< vector<int> >& vec_minus) = &gen_hierarchy;
  def("gen_hierarchy", expt_gen_hierarchy_v1);






  int (*expt_compute_nn_tot_v1)(int nquant, int KK, int LL) = &compute_nn_tot;
  def("compute_nn_tot", expt_compute_nn_tot_v1);

  vector<int> (*expt_allocate_1D_v1)
  (int sz1) = &allocate_1D;
  def("allocate_1D", expt_allocate_1D_v1);

  vector< vector<int> > (*expt_allocate_2D_v1)
  (int sz1, int sz2) = &allocate_2D;
  def("allocate_2D", expt_allocate_2D_v1);

  vector< vector< vector<int> > > (*expt_allocate_3D_v1)
  (int sz1, int sz2, int sz3) = &allocate_3D;
  def("allocate_3D", expt_allocate_3D_v1);


  void (*expt_compute_nn_v1)
  (int nquant, int KK, int LL, vector<int>& map_sum, 
   vector< vector< vector<int> > >& nn) = &compute_nn;
  def("compute_nn", expt_compute_nn_v1);


  void (*expt_compute_nn_sum_L_v1)
  (int nquant, int KK, int L, int& n_beg, int& n_end, 
   vector< vector< vector<int> > >& nn) = &compute_nn_sum_L;
  def("compute_nn_sum_L", expt_compute_nn_sum_L_v1);



  vector< vector<int> > (*expt_index_int2vec_v1)
  (vector< vector< vector<int> > >& nn, int n, int nquant, int KK) = &index_int2vec;
  def("index_int2vec", expt_index_int2vec_v1);


  int (*expt_sum2D_v1)
  (vector< vector<int> >& nvec) = &sum2D;
  def("sum2D", expt_sum2D_v1);

  int (*expt_index_vec2int_v1)
  (vector< vector< vector<int> > >& nn, vector< vector<int> >& nvec, int LL) = &index_vec2int;
  def("index_vec2int", expt_index_vec2int_v1);

  void (*expt_compute_map_v1)
  (int nquant, int KK, int LL, vector< vector< vector<int> > >& nn,
   vector< vector< vector<int> > >& map_nplus,
   vector< vector< vector<int> > >& map_nneg) = &compute_map;
  def("compute_map", expt_compute_map_v1);


  vector<int> (*expt_filter_v1)
  (vector<CMATRIX>& rho, double tolerance) = &filter;
  def("filter", expt_filter_v1);


  vector<CMATRIX> (*expt_initialize_el_phonon_couplings_v1)
  (int nquant) = &initialize_el_phonon_couplings;
  def("initialize_el_phonon_couplings", expt_initialize_el_phonon_couplings_v1);

  complex<double> (*expt_compute_matsubara_sum_v1)
  (vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara, int KK) = &compute_matsubara_sum;
  def("compute_matsubara_sum", expt_compute_matsubara_sum_v1);

  CMATRIX (*expt_compute_deriv_n_v1)
  (int n, vector<CMATRIX>& rho, CMATRIX& Ham, vector<CMATRIX>& projectors,
   double eta, double temperature,
   vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara,
   int do_truncate, int do_scale,
   vector< vector< vector<int> > >& nn, int KK, vector<int>& zero,
   vector< vector< vector<int> > >& map_nplus, vector< vector< vector<int> > >& map_nneg        
  ) = &compute_deriv_n;
  def("compute_deriv_n", expt_compute_deriv_n_v1);


  void (*expt_unpack_rho_v1)
  (vector<CMATRIX>& rho_unpacked, CMATRIX& RHO) = &unpack_rho;
  def("unpack_rho", expt_unpack_rho_v1);

  void (*expt_pack_rho_v1)
  (vector<CMATRIX>& rho_unpacked, CMATRIX& RHO) = &pack_rho;
  def("pack_rho", expt_pack_rho_v1);

  CMATRIX (*expt_compute_heom_derivatives_v1)
  (CMATRIX& RHO, bp::dict prms) = &compute_heom_derivatives;
  def("compute_heom_derivatives", expt_compute_heom_derivatives_v1);


  void (*expt_setup_bath_v1)
  (bp::dict params, vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara) = &setup_bath;
  def("setup_bath", expt_setup_bath_v1);


  class_< intList3 >("intList3")
      .def(vector_indexing_suite< intList3 >())
  ;


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygheom){
#else
BOOST_PYTHON_MODULE(libheom){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_heom_objects();

}


}// namespace libheom
}// namespace libdyn
}// liblibra
