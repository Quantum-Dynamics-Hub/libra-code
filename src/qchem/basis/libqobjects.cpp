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


  double (*expt_gaussian_overlap_AO_v1)
  (AO& AOa, AO& AOb) = &gaussian_overlap;

  double (*expt_gaussian_overlap_AO_v2)
  (AO& AOa, AO& AOb,int is_normalize) = &gaussian_overlap;

  boost::python::list (*expt_gaussian_overlap_AO_v3)
  (AO& AOa, AO& AOb,int is_normalize, int is_derivs) = &gaussian_overlap;





  // ============ Now export functions =============
  def("gaussian_overlap", expt_gaussian_overlap_G_v1);
  def("gaussian_overlap", expt_gaussian_overlap_G_v2);
  def("gaussian_overlap", expt_gaussian_overlap_G_v3);

  def("gaussian_overlap", expt_gaussian_overlap_AO_v1);
  def("gaussian_overlap", expt_gaussian_overlap_AO_v2);
  def("gaussian_overlap", expt_gaussian_overlap_AO_v3);




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
      .def("norm1",&PrimitiveG::norm1)
      .def("normalization_factor",&PrimitiveG::normalization_factor)
      .def("show_info",&PrimitiveG::show_info)

      .def("shift_position",&PrimitiveG::shift_position)

  ;

  class_< PrimitiveGList >("PrimitiveGList")
      .def(vector_indexing_suite< PrimitiveGList >())
  ;



  class_<AO>("AO",init<>())
      .def(init<const AO&>())
      .def("__copy__", &generic__copy__<AO>) 
      .def("__deepcopy__", &generic__deepcopy__<AO>)
      .def_readwrite("element",&AO::element)
      .def_readwrite("ao_shell",&AO::ao_shell)
      .def_readwrite("ao_shell_type",&AO::ao_shell_type)
      .def_readwrite("ao_name",&AO::ao_name)
      .def_readwrite("x_exp",&AO::x_exp)
      .def_readwrite("y_exp",&AO::y_exp)
      .def_readwrite("z_exp",&AO::z_exp)
      .def_readwrite("expansion_size",&AO::expansion_size)

      .def_readwrite("primitives",&AO::primitives)
      .def_readwrite("coefficients",&AO::coefficients)

      .def("clear",&AO::clear)
      .def("add_primitive",&AO::add_primitive)
      .def("show_info",&AO::show_info)

      .def("compute",&AO::compute)
      .def("norm2",&AO::norm2)
      .def("norm1",&AO::norm1)
      .def("normalization_factor",&AO::normalization_factor)
      .def("normalize",&AO::normalize)

      .def("shift_position",&AO::shift_position)

  ;

  class_< AOList >("AOList")
      .def(vector_indexing_suite< AOList >())
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

