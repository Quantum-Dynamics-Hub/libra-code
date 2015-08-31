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

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libconverters.h"
using namespace boost::python;


namespace libconverters{

//using namespace libuniverse;
//using namespace libmol;
//using namespace libchemsys;


void export_converters_objects(){

  void (*expt_system_to_nuclear_v1)(System& syst, Nuclear& nucl) = &system_to_nuclear;
  void (*expt_nuclear_to_system_v1)(Nuclear& nucl, System& syst) = &nuclear_to_system;


  def("system_to_nuclear", expt_system_to_nuclear_v1);
  def("nuclear_to_system", expt_nuclear_to_system_v1);

  

  class_<StringDoubleMap >("StringDoubleMap")
    .def(map_indexing_suite <StringDoubleMap >())
    .def("__len__", &StringDoubleMap::size)
    .def("clear", &StringDoubleMap::clear)
    .def("__getitem__", &map_item<StringDoubleMap>::get,  return_value_policy<copy_non_const_reference>())
    .def("__setitem__", &map_item<StringDoubleMap>::set,  with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &map_item<StringDoubleMap>::del)
    .def("in", &map_item<StringDoubleMap>::in)
    .def("items", &items<StringDoubleMap>)
    .def("keys", &keys<StringDoubleMap>)
    .def("values", &values<StringDoubleMap>)
  ;

  class_<StringIntMap >("StringIntMap")
    .def(map_indexing_suite <StringIntMap >())
    .def("__len__", &StringIntMap::size)
    .def("clear", &StringIntMap::clear)
    .def("__getitem__", &map_item<StringIntMap>::get,  return_value_policy<copy_non_const_reference>())
    .def("__setitem__", &map_item<StringIntMap>::set,  with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &map_item<StringIntMap>::del)
    .def("in", &map_item<StringIntMap>::in)
    .def("items", &items<StringIntMap>)
    .def("keys", &keys<StringIntMap>)
    .def("values", &values<StringIntMap>)
  ;

  class_<StringVDoubleMap >("StringVDoubleMap")
    .def(map_indexing_suite <StringVDoubleMap >())
    .def("__len__", &StringVDoubleMap::size)
    .def("clear", &StringVDoubleMap::clear)
    .def("__getitem__", &map_item<StringVDoubleMap>::get,  return_value_policy<copy_non_const_reference>())
    .def("__setitem__", &map_item<StringVDoubleMap>::set,  with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &map_item<StringVDoubleMap>::del)
    .def("in", &map_item<StringVDoubleMap>::in)
    .def("items", &items<StringVDoubleMap>)
    .def("keys", &keys<StringVDoubleMap>)
    .def("values", &values<StringVDoubleMap>)
  ;




}// export_converters_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygconverters){
#else
BOOST_PYTHON_MODULE(libconverters){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
//  export_Mathematics_objects();
  export_converters_objects();

}


}// libconverters


