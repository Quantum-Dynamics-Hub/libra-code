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

  ///==================  Overlaps ===========================
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


  ///==================  Moments ===========================
  double (*expt_gaussian_moment_G_v1)
  ( PrimitiveG& GA, PrimitiveG& G, PrimitiveG& GB) = &gaussian_moment;
  double (*expt_gaussian_moment_G_v2)
  ( PrimitiveG& GA, PrimitiveG& G, PrimitiveG& GB,int is_normalize) = &gaussian_moment;
  boost::python::list (*expt_gaussian_moment_G_v3)
  ( PrimitiveG& GA, PrimitiveG& G, PrimitiveG& GB,int is_normalize, int is_derivs  ) = &gaussian_moment;

  double (*expt_gaussian_moment_AO_v1)
  (AO& AOa, PrimitiveG& G, AO& AOb) = &gaussian_moment;
  double (*expt_gaussian_moment_AO_v2)
  (AO& AOa, PrimitiveG& G, AO& AOb,int is_normalize) = &gaussian_moment;
  boost::python::list (*expt_gaussian_moment_AO_v3)
  (AO& AOa, PrimitiveG& G, AO& AOb,int is_normalize, int is_derivs) = &gaussian_moment;


  ///==================  Pseudopotentials ===========================
  double (*expt_pseudopot02_G_v1)
  (double C0, double C2, double alp, const VECTOR& R, PrimitiveG& GA, PrimitiveG& GB) = &pseudopot02;
  double (*expt_pseudopot02_G_v2)
  (double C0, double C2, double alp, const VECTOR& R, PrimitiveG& GA, PrimitiveG& GB,int is_normalize) = &pseudopot02;
  boost::python::list (*expt_pseudopot02_G_v3)
  (double C0, double C2, double alp, const VECTOR& R, PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs) = &pseudopot02;


  double (*expt_pseudopot02_AO_v1)
  (double C0, double C2, double alp, const VECTOR& R, AO& AOa, AO& AOb) = &pseudopot02;
  double (*expt_pseudopot02_AO_v2)
  (double C0, double C2, double alp, const VECTOR& R, AO& AOa, AO& AOb,int is_normalize) = &pseudopot02;
  boost::python::list (*expt_pseudopot02_AO_v3)
  (double C0, double C2, double alp, const VECTOR& R, AO& AOa, AO& AOb,int is_normalize, int is_derivs) = &pseudopot02;


  ///==================  Multipoles ===========================
  VECTOR (*expt_transition_dipole_moment_G_v1)
  (PrimitiveG& GA, PrimitiveG& GB) = &transition_dipole_moment;
  VECTOR (*expt_transition_dipole_moment_G_v2)
  (PrimitiveG& GA, PrimitiveG& GB, int is_normalize) = &transition_dipole_moment;
  boost::python::list (*expt_transition_dipole_moment_G_v3)
  (PrimitiveG& GA, PrimitiveG& GB, int is_normalize, int is_derivs) = &transition_dipole_moment;

  VECTOR (*expt_transition_dipole_moment_AO_v1)
  (AO& AOa, AO& AOb) = &transition_dipole_moment;
  VECTOR (*expt_transition_dipole_moment_AO_v2)
  (AO& AOa, AO& AOb, int is_normalize) = &transition_dipole_moment;
  boost::python::list (*expt_transition_dipole_moment_AO_v3)
  (AO& AOa, AO& AOb, int is_normalize, int is_derivs) = &transition_dipole_moment;


  ///==================  Derivative coupling integrals ===========================
  VECTOR (*expt_derivative_coupling_integral_G_v1)
  (PrimitiveG& GA, PrimitiveG& GB) = &derivative_coupling_integral;
  VECTOR (*expt_derivative_coupling_integral_G_v2)
  (PrimitiveG& GA, PrimitiveG& GB, int is_normalize) = &derivative_coupling_integral;
  boost::python::list (*expt_derivative_coupling_integral_G_v3)
  (PrimitiveG& GA, PrimitiveG& GB, int is_normalize, int is_derivs) = &derivative_coupling_integral;

  VECTOR (*expt_derivative_coupling_integral_AO_v1)
  (AO& AOa, AO& AOb) = &derivative_coupling_integral;
  VECTOR (*expt_derivative_coupling_integral_AO_v2)
  (AO& AOa, AO& AOb, int is_normalize) = &derivative_coupling_integral;
  boost::python::list (*expt_derivative_coupling_integral_AO_v3)
  (AO& AOa, AO& AOb, int is_normalize, int is_derivs) = &derivative_coupling_integral;



  ///==================  Kinetic integrals ===========================
  double (*expt_kinetic_integral_G_v1)
  ( PrimitiveG& GA, PrimitiveG& GB) = &kinetic_integral;
  double (*expt_kinetic_integral_G_v2)
  ( PrimitiveG& GA, PrimitiveG& GB,int is_normalize) = &kinetic_integral;
  boost::python::list (*expt_kinetic_integral_G_v3)
  ( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs  ) = &kinetic_integral;

  double (*expt_kinetic_integral_AO_v1)
  (AO& AOa, AO& AOb) = &kinetic_integral;
  double (*expt_kinetic_integral_AO_v2)
  (AO& AOa, AO& AOb,int is_normalize) = &kinetic_integral;
  boost::python::list (*expt_kinetic_integral_AO_v3)
  (AO& AOa, AO& AOb,int is_normalize, int is_derivs) = &kinetic_integral;


  ///==================  NAI ===========================
  double (*expt_nuclear_attraction_integral_G_v1)
  ( PrimitiveG& GA, PrimitiveG& GB, VECTOR& Rc  ) = &nuclear_attraction_integral;
  double (*expt_nuclear_attraction_integral_G_v2)
  ( PrimitiveG& GA, PrimitiveG& GB, VECTOR& Rc, int is_normalize ) = &nuclear_attraction_integral;
  boost::python::list (*expt_nuclear_attraction_integral_G_v3)
  ( PrimitiveG& GA, PrimitiveG& GB, VECTOR& Rc, int is_normalize, int is_derivs  ) = &nuclear_attraction_integral;

  double (*expt_nuclear_attraction_integral_AO_v1)
  ( AO& AOa, AO& AOb, VECTOR& Rc  ) = &nuclear_attraction_integral;
  double (*expt_nuclear_attraction_integral_AO_v2)
  ( AO& AOa, AO& AOb, VECTOR& Rc,  int is_normalize ) = &nuclear_attraction_integral;
  boost::python::list (*expt_nuclear_attraction_integral_AO_v3)
  ( AO& AOa, AO& AOb, VECTOR& Rc,  int is_normalize, int is_derivs  ) = &nuclear_attraction_integral;



  ///==================  ERI ===========================
  double (*expt_electron_repulsion_integral_G_v1)
  ( PrimitiveG& GA, PrimitiveG& GB, PrimitiveG& GC, PrimitiveG& GD
  ) = &electron_repulsion_integral;
  double (*expt_electron_repulsion_integral_G_v2)
  ( PrimitiveG& GA, PrimitiveG& GB, PrimitiveG& GC, PrimitiveG& GD,
    int is_normalize ) = &electron_repulsion_integral;
  boost::python::list (*expt_electron_repulsion_integral_G_v3)
  ( PrimitiveG& GA, PrimitiveG& GB, PrimitiveG& GC, PrimitiveG& GD,
    int is_normalize, int is_derivs  ) = &electron_repulsion_integral;

  double (*expt_electron_repulsion_integral_AO_v1)
  ( AO& AOa, AO& AOb, AO& AOc, AO& AOd  ) = &electron_repulsion_integral;
  double (*expt_electron_repulsion_integral_AO_v2)
  ( AO& AOa, AO& AOb, AO& AOc, AO& AOd, int is_normalize ) = &electron_repulsion_integral;
  boost::python::list (*expt_electron_repulsion_integral_AO_v3)
  ( AO& AOa, AO& AOb, AO& AOc, AO& AOd, int is_normalize, int is_derivs) = &electron_repulsion_integral;



  // ============ Now export functions =============
  // Overlaps
  def("gaussian_overlap", expt_gaussian_overlap_G_v1);
  def("gaussian_overlap", expt_gaussian_overlap_G_v2);
  def("gaussian_overlap", expt_gaussian_overlap_G_v3);

  def("gaussian_overlap", expt_gaussian_overlap_AO_v1);
  def("gaussian_overlap", expt_gaussian_overlap_AO_v2);
  def("gaussian_overlap", expt_gaussian_overlap_AO_v3);

  // Moments
  def("gaussian_moment", expt_gaussian_moment_G_v1);
  def("gaussian_moment", expt_gaussian_moment_G_v2);
  def("gaussian_moment", expt_gaussian_moment_G_v3);

  def("gaussian_moment", expt_gaussian_moment_AO_v1);
  def("gaussian_moment", expt_gaussian_moment_AO_v2);
  def("gaussian_moment", expt_gaussian_moment_AO_v3);

  // Pseudopotentials
  def("pseudopot02", expt_pseudopot02_G_v1);
  def("pseudopot02", expt_pseudopot02_G_v2);
  def("pseudopot02", expt_pseudopot02_G_v3);

  def("pseudopot02", expt_pseudopot02_AO_v1);
  def("pseudopot02", expt_pseudopot02_AO_v2);
  def("pseudopot02", expt_pseudopot02_AO_v3);

  // Multipoles
  def("transition_dipole_moment", expt_transition_dipole_moment_G_v1);
  def("transition_dipole_moment", expt_transition_dipole_moment_G_v2);
  def("transition_dipole_moment", expt_transition_dipole_moment_G_v3);

  def("transition_dipole_moment", expt_transition_dipole_moment_AO_v1);
  def("transition_dipole_moment", expt_transition_dipole_moment_AO_v2);
  def("transition_dipole_moment", expt_transition_dipole_moment_AO_v3);


  // Derivative coupling integrals
  def("derivative_coupling_integral", expt_derivative_coupling_integral_G_v1);
  def("derivative_coupling_integral", expt_derivative_coupling_integral_G_v2);
  def("derivative_coupling_integral", expt_derivative_coupling_integral_G_v3);

  def("derivative_coupling_integral", expt_derivative_coupling_integral_AO_v1);
  def("derivative_coupling_integral", expt_derivative_coupling_integral_AO_v2);
  def("derivative_coupling_integral", expt_derivative_coupling_integral_AO_v3);


  // Kinetic integrals
  def("kinetic_integral", expt_kinetic_integral_G_v1);
  def("kinetic_integral", expt_kinetic_integral_G_v2);
  def("kinetic_integral", expt_kinetic_integral_G_v3);

  def("kinetic_integral", expt_kinetic_integral_AO_v1);
  def("kinetic_integral", expt_kinetic_integral_AO_v2);
  def("kinetic_integral", expt_kinetic_integral_AO_v3);


  // NAIs
  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_G_v1);
  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_G_v2);
  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_G_v3);

  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_AO_v1);
  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_AO_v2);
  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_AO_v3);



  // ERIs
  def("electron_repulsion_integral", expt_electron_repulsion_integral_G_v1);
  def("electron_repulsion_integral", expt_electron_repulsion_integral_G_v2);
  def("electron_repulsion_integral", expt_electron_repulsion_integral_G_v3);

  def("electron_repulsion_integral", expt_electron_repulsion_integral_AO_v1);
  def("electron_repulsion_integral", expt_electron_repulsion_integral_AO_v2);
  def("electron_repulsion_integral", expt_electron_repulsion_integral_AO_v3);



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
      .def("set_position",&AO::set_position)

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

