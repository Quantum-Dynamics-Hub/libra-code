/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libmmath.cpp
  \brief The file describes Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libmmath.h"

using namespace boost::python;

/// libmmath namespace
namespace libmmath{

using namespace libspecialfunctions;
using namespace liblinalg;
using namespace libgraph;
using namespace liboperators;
using namespace librandom;
using namespace libdata;
using namespace libann;
using namespace libmeigen;
using namespace libsymmetry;

void export_Mathematics_objects(){
/** 
  \brief Exporter of Timer class and other mathematical libraries and their components

*/

  export_linalg_objects();
  export_SpecialFunctions_objects();
  export_GRAPH_objects();
  export_Operators_objects();
  export_Random_objects();
  export_Data_objects();
  export_NeuralNetwork_objects();
  export_mEigen_objects();
  export_symmetry_objects();


  class_<Timer>("Timer",init<>())
      .def("__copy__", &generic__copy__<Timer>)
      .def("__deepcopy__", &generic__deepcopy__<Timer>)

      .def("start", &Timer::start)
      .def("stop", &Timer::stop)
      .def("show", &Timer::show)

  ;




}// export_Mathematics_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmmath){
#else
BOOST_PYTHON_MODULE(libmmath){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
  export_Mathematics_objects();

}


}// libmmath

