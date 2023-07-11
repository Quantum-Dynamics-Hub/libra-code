/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file liboperators.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "liboperators.h"



/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// liboperators namespace
namespace liboperators{


void export_Operators_objects(){
/** 
  \brief Exporter of liboperators classes and functions

*/


  boost::python::list (*expt_rotate1)(double,double,double) = &expt_rotate;
  double (*expt_shift1)(double, double) = &expt_shift;
  double (*expt_scale1)(double, double) = &expt_scale;


  def("rotate", expt_rotate1);
  def("shift", expt_shift1);
  def("scale", expt_scale1);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygoperators){
#else
BOOST_PYTHON_MODULE(liboperators){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Operators_objects();

}

}// liboperators
}// liblibra


