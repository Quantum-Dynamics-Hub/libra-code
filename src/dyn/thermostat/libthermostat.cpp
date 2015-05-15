#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libthermostat.h"

using namespace boost::python;


namespace libdyn{
namespace libthermostat{



void export_Thermostat_objects(){

  class_<Thermostat>("Thermostat",init<>())
      .def(init<boost::python::dict>())
      .def_readwrite("s_var",&Thermostat::s_var)
      .def_readwrite("Ps",&Thermostat::Ps)
      .def_readwrite("Q",&Thermostat::Q)
      .def_readwrite("NHC_size",&Thermostat::NHC_size)
      .def_readwrite("nu_therm",&Thermostat::nu_therm)
      .def_readwrite("Temperature",&Thermostat::Temperature)
      .def_readwrite("thermostat_type",&Thermostat::thermostat_type)

      .def("show_info",&Thermostat::show_info)
  ;

}// export_Thermostat_objects



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygthermostat){
#else
BOOST_PYTHON_MODULE(libthermostat){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Thermostat_objects();

}


}// namespace libthermostat
}// namespace libdyn

