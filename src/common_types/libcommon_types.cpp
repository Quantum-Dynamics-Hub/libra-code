/*********************************************************************************
* Copyright (C) 2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libcommon_types.cpp
  \brief The file implements Python export function
    
*/

#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libcommon_types.h"



/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libcommon_types namespace
namespace libcommon_types{


void export_common_types_objects(){
/** 
  \brief Exporter of the libcommon_types classes and functions

*/

  class_<excitation>("excitation",init<>())
      .def(init<int,int,int,int>())
//      .def("__copy__", &generic__copy__<excitation>)
//      .def("__deepcopy__", &generic__deepcopy__<excitation>)

      .def_readwrite("size", &excitation::size)
      .def_readwrite("from_orbit", &excitation::from_orbit)
      .def_readwrite("from_spin", &excitation::from_spin)
      .def_readwrite("to_orbit", &excitation::to_orbit)
      .def_readwrite("to_spin", &excitation::to_spin)

      .def_readwrite("Energy", &excitation::Energy)
      .def_readwrite("f", &excitation::f)

  ;


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygcommon_types){
#else
BOOST_PYTHON_MODULE(libcommon_types){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_common_types_objects();

}



}// namespace libcommon_types
}// namespace liblibra


