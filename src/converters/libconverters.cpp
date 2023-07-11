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

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libconverters.h"

/// liblibra namespace
namespace liblibra{


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


  vector<int> (*expt_Py2Cpp_v1)(boost::python::list x) = &Py2Cpp;
  vector<double> (*expt_Py2Cpp_v2)(boost::python::list x) = &Py2Cpp;
  vector<complex<double> > (*expt_Py2Cpp_v3)(boost::python::list x) = &Py2Cpp;
  vector<std::string> (*expt_Py2Cpp_v4)(boost::python::list x) = &Py2Cpp;
  vector<VECTOR> (*expt_Py2Cpp_v5)(boost::python::list x) = &Py2Cpp;
  vector<MATRIX> (*expt_Py2Cpp_v6)(boost::python::list x) = &Py2Cpp;
  vector<CMATRIX> (*expt_Py2Cpp_v7)(boost::python::list x) = &Py2Cpp;
//  vector<listHamiltonian_QM> (*expt_Py2Cpp_v7)(boost::python::list x) = &Py2Cpp;

  def("Py2Cpp_int", expt_Py2Cpp_v1);
  def("Py2Cpp_double", expt_Py2Cpp_v2);
  def("Py2Cpp_complex", expt_Py2Cpp_v3);
  def("Py2Cpp_string", expt_Py2Cpp_v4);
  def("Py2Cpp_VECTOR", expt_Py2Cpp_v5);
  def("Py2Cpp_MATRIX", expt_Py2Cpp_v6);
  def("Py2Cpp_CMATRIX", expt_Py2Cpp_v7);
//  def("Py2Cpp_listHamiltonian_QM", expt_Py2Cpp_v7);

  boost::python::list (*expt_Cpp2Py_v1)(vector<int>& x) = &Cpp2Py;
  boost::python::list (*expt_Cpp2Py_v2)(vector<double>& x) = &Cpp2Py;
  boost::python::list (*expt_Cpp2Py_v3)(vector<complex<double> >& x) = &Cpp2Py;
  boost::python::list (*expt_Cpp2Py_v4)(vector<std::string>& x) = &Cpp2Py;
  boost::python::list (*expt_Cpp2Py_v5)(vector<VECTOR>& x) = &Cpp2Py;
  boost::python::list (*expt_Cpp2Py_v6)(vector<MATRIX>& x) = &Cpp2Py;
  boost::python::list (*expt_Cpp2Py_v7)(vector<CMATRIX>& x) = &Cpp2Py;
//  boost::python::list (*expt_Cpp2Py_v7)(vector<listHamiltonian_QM>& x) = &Cpp2Py;

  def("Cpp2Py", expt_Cpp2Py_v1);
  def("Cpp2Py", expt_Cpp2Py_v2);
  def("Cpp2Py", expt_Cpp2Py_v3);
  def("Cpp2Py", expt_Cpp2Py_v4);
  def("Cpp2Py", expt_Cpp2Py_v5);
  def("Cpp2Py", expt_Cpp2Py_v6);
  def("Cpp2Py", expt_Cpp2Py_v7);


  class_< StringList >("StringList")
    .def(vector_indexing_suite< StringList >())
//    .def("__len__", &StringList::size)
//    .def("clear", &StringList::clear)
//    .def("append", &std_item<StringList>::add,  with_custodian_and_ward<1,2>()) // to let container keep value
//    .def("__getitem__", &std_item<StringList>::get, return_value_policy<copy_non_const_reference>())
//    .def("__setitem__", &std_item<StringList>::set, with_custodian_and_ward<1,2>()) // to let container keep value
//    .def("__delitem__", &std_item<StringList>::del)
//    .def("index", &std_item<StringList>::index)
//    .def("in", &std_item<StringList>::in)
//    .def("__iter__", boost::python::iterator<StringList>())
  ;

  class_< StringMap >("StringMap")
      .def(vector_indexing_suite< StringMap >())
  ;

  

  class_<StringDoubleMap >("StringDoubleMap")
    .def(map_indexing_suite <StringDoubleMap >())
    .def("__len__", &StringDoubleMap::size)
    .def("clear", &StringDoubleMap::clear)
    .def("__getitem__", &map_item<StringDoubleMap>::get,  return_value_policy<copy_non_const_reference>())
    .def("__setitem__", &map_item<StringDoubleMap>::set,  with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &map_item<StringDoubleMap>::del)
    .def("index", &map_item<StringDoubleMap>::index)
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
    .def("index", &map_item<StringIntMap>::index)
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
}// liblibra

