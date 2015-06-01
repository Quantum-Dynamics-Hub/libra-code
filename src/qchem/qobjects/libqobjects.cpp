#define BOOST_PYTHON_MAX_ARITY 30
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libqobjects.h"

using namespace boost::python;
using namespace libmmath;


namespace libqchem{
namespace libqobjects{




void export_qobjects_objects(){


  double (*expt_gaussian_overlap_G_v1)
  ( PrimitiveG& GA, PrimitiveG& GB) = &gaussian_overlap;

  double (*expt_gaussian_overlap_G_v2)
  ( PrimitiveG& GA, PrimitiveG& GB,int is_normalize) = &gaussian_overlap;

  boost::python::list (*expt_gaussian_overlap_G_v3)
  ( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs  ) = &gaussian_overlap;



  // ============ Now export functions =============
  def("gaussian_overlap", expt_gaussian_overlap_G_v1);
  def("gaussian_overlap", expt_gaussian_overlap_G_v2);
  def("gaussian_overlap", expt_gaussian_overlap_G_v3);



  class_<PrimitiveG>("PrimitiveG",init<>())
      .def(init<int,int,int,double,VECTOR&>()) 
      .def(init<const PrimitiveG&>())
      .def("__copy__", &generic__copy__<PrimitiveG>) 
      .def("__deepcopy__", &generic__deepcopy__<PrimitiveG>)
      .def_readwrite("x_exp",&PrimitiveG::x_exp)
      .def_readwrite("y_exp",&PrimitiveG::y_exp)
      .def_readwrite("z_exp",&PrimitiveG::z_exp)
      .def_readwrite("alpha",&PrimitiveG::alpha)
      .def_readwrite("R",&PrimitiveG::R)
      .def_readwrite("value",&PrimitiveG::value)

      .def("init",&PrimitiveG::init)
      .def("get_x_exp",&PrimitiveG::get_x_exp)
      .def("get_y_exp",&PrimitiveG::get_y_exp)
      .def("get_z_exp",&PrimitiveG::get_z_exp)
      .def("get_alpha",&PrimitiveG::get_alpha)
      .def("get_R",&PrimitiveG::get_R)
      .def("get_value",&PrimitiveG::get_value)

      .def("set_x_exp",&PrimitiveG::set_x_exp)
      .def("set_y_exp",&PrimitiveG::set_y_exp)
      .def("set_z_exp",&PrimitiveG::set_z_exp)
      .def("set_alpha",&PrimitiveG::set_alpha)
      .def("set_R",&PrimitiveG::set_R)

      .def("compute",&PrimitiveG::compute)
      .def("norm2",&PrimitiveG::norm2)
      .def("norm2",&PrimitiveG::norm1)
      .def("normalization_factor",&PrimitiveG::normalization_factor)
      .def("show_info",&PrimitiveG::show_info)

      .def("shift_position",&PrimitiveG::shift_position)


  ;


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygqobjects){
#else
BOOST_PYTHON_MODULE(libqobjects){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_qobjects_objects();

}


}// namespace libqobject
}// namespace libqchem

