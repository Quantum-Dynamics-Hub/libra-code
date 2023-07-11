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
  \file libsymmetry.cpp
  \brief The file describes Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#endif

#include "libsymmetry.h"



/// liblibra namespace
namespace liblibra{

using namespace boost::python;
using namespace liblinalg;

/// libsymmetry namespace
namespace libsymmetry{


void export_symmetry_objects(){
/** 
  \brief Exporter of libsymmetry classes and functions

*/


  void (*expt_Apply_Symmetry_v1)(std::string space_group_name,VECTOR r,std::vector<VECTOR>& r_equiv) = &Apply_Symmetry;
//  voi (*expt_Apply_Symmetry_v1)(std::string space_group_name,VECTOR r,std::vector<VECTOR>& r_equiv) = &Apply_Symmetry;

  class_<SPACE_GROUP>("SPACE_GROUP",init<>())
      .def(init<std::string>())
      .def("__copy__", &generic__copy__<SPACE_GROUP>) 
      .def("__deepcopy__", &generic__deepcopy__<SPACE_GROUP>)
      .def_readwrite("operators",&SPACE_GROUP::operators)

  ;

  def("Apply_Symmetry",expt_Apply_Symmetry_v1);



}// export_symmetry_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygsymmetry){
#else
BOOST_PYTHON_MODULE(libsymmetry){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_symmetry_objects();

}




}// namespace libsymmetry
}// namespace liblibra




