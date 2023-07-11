/*********************************************************************************
* Copyright (C) 2018-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libopt.cpp
  \brief The file implements Python export function
    
*/


#if defined(USING_PCH)
#include "../pch.h"
#else
#include <stdlib.h>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libopt.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libopt namespace
namespace libopt{

using namespace boost::python;
namespace bp = boost::python;



void export_opt_objects(){
/** 
  \brief Exporter of libopt classes and functions

*/


  MATRIX (*expt_grad_descent_v1)
  (bp::object grad_function, MATRIX& dof, bp::object funct_params, 
   double grad_tol, double step_size, int max_steps) = &grad_descent;

  def("grad_descent",expt_grad_descent_v1);


}// export_opt_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygopt){
#else
BOOST_PYTHON_MODULE(libopt){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_opt_objects();

}


}// libopt
}// liblibra

