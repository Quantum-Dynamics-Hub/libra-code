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
  \file libsolvers.cpp
  \brief This file implements the exprots of libsolvers objects to Python
        
*/


#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libsolvers.h"

/// liblibra namespace
namespace liblibra{

/// libsolvers namespace
namespace libsolvers{


void export_solvers_objects(){
/** 
  \brief Exporter of libsolvers classes and functions

*/


  //----------------- DIIS.cpp ------------------------------

  void (DIIS::*expt_add_diis_matrices_v1)(MATRIX& X, MATRIX& err) = &DIIS::add_diis_matrices;
//  void (DIIS::*expt_update_diis_coefficients_v1)() = &DIIS::update_diis_coefficients;
  void (DIIS::*expt_extrapolate_matrix_v1)(MATRIX& X) = &DIIS::extrapolate_matrix;

  class_<DIIS>("DIIS",init<int,int>())
      .def("__copy__", &generic__copy__<DIIS>)
      .def("__deepcopy__", &generic__deepcopy__<DIIS>)

      .def("get_diis_X", &DIIS::get_diis_X)
      .def("get_diis_err", &DIIS::get_diis_err)
      .def("get_diis_c", &DIIS::get_diis_c)
      .def("add_diis_matrices",expt_add_diis_matrices_v1)
//      .def("update_diis_coefficients",expt_update_diis_coefficients_v1)
      .def("extrapolate_matrix",expt_extrapolate_matrix_v1)


      .def_readwrite("N_diis_max",&DIIS::N_diis_max)
      .def_readwrite("N_diis",&DIIS::N_diis)
      .def_readwrite("N_diis_eff",&DIIS::N_diis_eff)

  ;


}// export_solvers_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygsolvers){
#else
BOOST_PYTHON_MODULE(libsolvers){
#endif

  export_solvers_objects();

}


}// namespace libsolvers
}// liblibra



