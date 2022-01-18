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
  \file libqtag.cpp
  \brief The file implements Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libqtag.h"


/// liblibra namespace
namespace liblibra{


using namespace boost::python;

/// libdyn namespace
namespace libdyn{

/// libqtag namespace
namespace libqtag{


void export_qtag_objects(){
/** 
  \brief Exporter of libqtag classes and functions

*/

  CMATRIX (*expt_qtag_psi_v1)
  (MATRIX q, MATRIX& q1, MATRIX& p1, MATRIX& alp1, MATRIX& s1, CMATRIX& Coeff) = &qtag_psi;
  def("qtag_psi",  expt_qtag_psi_v1);


  CMATRIX (*expt_qtag_overlap_elementary_v1)
  (MATRIX& q1, MATRIX& p1, MATRIX& alp1, MATRIX& s1,
   MATRIX& q2, MATRIX& p2, MATRIX& alp2, MATRIX& s2) = &qtag_overlap_elementary;
  def("qtag_overlap_elementary",  expt_qtag_overlap_elementary_v1);

  CMATRIX (*expt_qtag_kinetic_elementary_v1)
  (MATRIX q, MATRIX& p, MATRIX& alp, MATRIX& s, MATRIX& invM) = &qtag_kinetic_elementary;
  def("qtag_kinetic_elementary", expt_qtag_kinetic_elementary_v1);

  CMATRIX (*expt_qtag_overlap_v1)
  (vector<int>& active_states, CMATRIX& ovlp, int nstates) = &qtag_overlap;
  def("qtag_overlap", expt_qtag_overlap_v1);


  CMATRIX (*expt_qtag_hamiltonian_v1)
  (MATRIX q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff,
   vector<int>& active_states, CMATRIX& ovlp, CMATRIX& kin,
   MATRIX& invM, nHamiltonian& ham, bp::object compute_ham_funct,  
   bp::dict& compute_ham_params) = &qtag_hamiltonian;
  def("qtag_hamiltonian", expt_qtag_hamiltonian_v1);


  CMATRIX (*expt_qtag_momentum_v1)(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff) = &qtag_momentum;
  def("qtag_momentum",  expt_qtag_momentum_v1);



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygqtag){
#else
BOOST_PYTHON_MODULE(libqtag){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_qtag_objects();

}


}// namespace libqtag
}// namespace libdyn
}// liblibra

