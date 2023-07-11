/*********************************************************************************
* Copyright (C) 2019-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libintegrators.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libintegrators.h"


/// liblibra namespace
namespace liblibra{


/// libdyn namespace
namespace libintegrators{

namespace bp = boost::python;


void export_integrators_objects(){
/** 
  \brief Exporter of libintegrators classes and functions

*/

  CMATRIX (*expt_RK4_v1)
  (CMATRIX& q, double dt, bp::object compute_derivatives, bp::object function_params) = &RK4; 

  def("RK4", expt_RK4_v1);


}// export_integrators_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygintegrators){
#else
BOOST_PYTHON_MODULE(libintegrators){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_integrators_objects();

}


}// libintegrators
}// liblibra

