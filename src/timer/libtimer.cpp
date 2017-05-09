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
  \file libtimer.cpp
  \brief The file describes Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libtimer.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;


void export_timer_objects(){
/** 
  \brief Exporter of Timer class and other mathematical libraries and their components

*/


  class_<Timer>("Timer",init<>())
//      .def("__copy__", &generic__copy__<Timer>)
//      .def("__deepcopy__", &generic__deepcopy__<Timer>)

      .def("start", &Timer::start)
      .def("stop", &Timer::stop)
      .def("show", &Timer::show)

  ;




}// export_timer_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygtimer){
#else
BOOST_PYTHON_MODULE(libtimer){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
  export_timer_objects();

}


}// liblibra

