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

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

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

  int (*expt_compute_nn_tot_v1)(int d, int max_tier) = &compute_nn_tot;
  def("compute_nn_tot", expt_compute_nn_tot_v1); 

  vector< vector<int> > (*expt_gen_next_level_v1)
  (vector< vector<int> >& parents) = &gen_next_level;
  def("gen_next_level", expt_gen_next_level_v1);

  void (*expt_gen_hierarchy_v1)
  (int d, int max_tier, int verbosity,
   vector< vector<int> >& all_vectors,  
   vector< vector<int> >& vec_plus, 
   vector< vector<int> >& vec_minus) = &gen_hierarchy;
  def("gen_hierarchy", expt_gen_hierarchy_v1);




  vector<int> (*expt_filter_v1)
  (vector<CMATRIX>& rho, vector<int>& adm_list, double tolerance, int do_zeroing) = &filter;
  def("filter", expt_filter_v1);


  CMATRIX (*expt_compute_deriv_n_v1)
  (int n, vector<CMATRIX>& rho, CMATRIX& Ham, vector<CMATRIX>& el_phon_coupl,
   double eta, double temperature,
   vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara,
   int truncation_scheme, complex<double> truncation_prefactor, int do_scale, vector<int>& nonzero,
   vector< vector<int> >& nvectors, vector< vector<int> >& vec_plus, vector< vector<int> >& vec_minus        
  ) = &compute_deriv_n;
  def("compute_deriv_n", expt_compute_deriv_n_v1);

  CMATRIX (*expt_compute_heom_derivatives_v1)
  (CMATRIX& RHO, bp::dict prms) = &compute_heom_derivatives;
  def("compute_heom_derivatives", expt_compute_heom_derivatives_v1);





  vector<CMATRIX> (*expt_initialize_el_phonon_couplings_v1)
  (int nquant) = &initialize_el_phonon_couplings;
  def("initialize_el_phonon_couplings", expt_initialize_el_phonon_couplings_v1);

  complex<double> (*expt_compute_matsubara_sum_v1)
  (vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara, int KK) = &compute_matsubara_sum;
  def("compute_matsubara_sum", expt_compute_matsubara_sum_v1);

  void (*expt_setup_bath_v1)
  (int KK, double eta, double gamma, double temperature, 
   vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara) = &setup_bath;
  def("setup_bath", expt_setup_bath_v1);




  void (*expt_unpack_mtx_v1)
  (vector<CMATRIX>& rho_unpacked, CMATRIX& RHO) = &unpack_mtx;
  def("unpack_mtx", expt_unpack_mtx_v1);

  void (*expt_pack_mtx_v1)
  (vector<CMATRIX>& rho_unpacked, CMATRIX& RHO) = &pack_mtx;
  def("pack_mtx", expt_pack_mtx_v1);

  void (*expt_scale_rho_v1)
  (vector<CMATRIX>& rho, vector<CMATRIX>& rho_scaled, bp::dict prms) = &scale_rho;
  def("scale_rho", expt_scale_rho_v1);

  void (*expt_inv_scale_rho_v1)
  (vector<CMATRIX>& rho, vector<CMATRIX>& rho_scaled, bp::dict prms) = &inv_scale_rho;
  def("inv_scale_rho", expt_inv_scale_rho_v1);



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
