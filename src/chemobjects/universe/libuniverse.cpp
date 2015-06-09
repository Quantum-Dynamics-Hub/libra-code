#define BOOST_PYTHON_MAX_ARITY 30
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libuniverse.h"

using namespace boost::python;
using namespace libmmath;


namespace libchemobjects{
namespace libuniverse{


void export_Universe_objects(){

  class_<Element>("Element",init<>())
      .def_readwrite("Elt_name",&Element::Elt_name)
      .def_readwrite("Elt_id",&Element::Elt_id)
      .def_readwrite("Elt_mass",&Element::Elt_mass)
      .def_readwrite("Elt_number",&Element::Elt_number)
      .def_readwrite("Elt_nucleus_charge",&Element::Elt_nucleus_charge)
      .def_readwrite("Elt_period",&Element::Elt_period)
      .def_readwrite("Elt_group",&Element::Elt_group)
      .def_readwrite("Elt_block",&Element::Elt_block)
      .def_readwrite("Elt_red",&Element::Elt_red)
      .def_readwrite("Elt_green",&Element::Elt_green)
      .def_readwrite("Elt_blue",&Element::Elt_blue)
      .def_readwrite("Elt_bond_length",&Element::Elt_bond_length)
      .def_readwrite("Elt_radius",&Element::Elt_radius)

      .def("set",&Element::set)
      .def("show_info",&Element::show_info)
  ;


  class_<Universe>("Universe",init<>())
      .def("Add_Element_To_Periodic_Table",&Universe::Add_Element_To_Periodic_Table)
      .def("Get_Element",&Universe::Get_Element)
  ;



}// export_


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyguniverse){
#else
BOOST_PYTHON_MODULE(libuniverse){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Universe_objects();

}



}// namespace libuniverse
}// namespace libchemobjects




