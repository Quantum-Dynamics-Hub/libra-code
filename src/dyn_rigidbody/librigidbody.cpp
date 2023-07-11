/*********************************************************************************
* Copyright (C) 2015-2021 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file librigidbody.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

#include "librigidbody.h"


/// liblibra namespace
namespace liblibra{

/// librigidbody namespace
namespace librigidbody{



void export_RigidBody_objects(){
/** 
  \brief Exporter of librigidbody classes and functions

*/


  int (RigidBody::*init1)(int,double*,VECTOR*) = &RigidBody::init;
  int (RigidBody::*init2)(int sz,boost::python::list masses,boost::python::list positions) = &RigidBody::init;


  int (RigidBody::*set_orientation1)(const MATRIX3x3&)                            = &RigidBody::set_orientation;
  int (RigidBody::*set_orientation2)(const VECTOR&, const VECTOR&, const VECTOR&) = &RigidBody::set_orientation;
  int (RigidBody::*set_orientation3)(const QUATERNION&)                           = &RigidBody::set_orientation;
  int (RigidBody::*set_angular_momentum1)(const VECTOR&)                                = &RigidBody::set_angular_momentum;
  int (RigidBody::*set_angular_momentum2)(const double&, const double&, const double&)  = &RigidBody::set_angular_momentum;
  int (RigidBody::*set_angular_velocity1)(const VECTOR&)                                = &RigidBody::set_angular_velocity;
  int (RigidBody::*set_angular_velocity2)(const double&, const double&, const double&)  = &RigidBody::set_angular_velocity;
  int (RigidBody::*scale_angular_1)(double)                 = &RigidBody::scale_angular_;
  int (RigidBody::*scale_angular_2)(double,double,double)   = &RigidBody::scale_angular_;
  int (RigidBody::*scale_angular_3)(const MATRIX3x3&)       = &RigidBody::scale_angular_;

  int (RigidBody::*scale_linear_1)(double)                  = &RigidBody::scale_linear_;
  int (RigidBody::*scale_linear_2)(double,double,double)    = &RigidBody::scale_linear_;
  int (RigidBody::*scale_linear_3)(const MATRIX3x3&)        = &RigidBody::scale_linear_;
  int (RigidBody::*scale_position1)(double)                 = &RigidBody::scale_position;
  int (RigidBody::*scale_position2)(double,double,double)   = &RigidBody::scale_position;
  int (RigidBody::*scale_position3)(const MATRIX3x3&)       = &RigidBody::scale_position;
  int (RigidBody::*apply_force1)(double)      = &RigidBody::apply_force;
  int (RigidBody::*apply_force2)(MATRIX3x3&)  = &RigidBody::apply_force;

  void (RigidBody::*Rotate1)(const MATRIX3x3& rot) = &RigidBody::Rotate;
  void (RigidBody::*Rotate2)(const MATRIX3x3& rot, const VECTOR& pivot) = &RigidBody::Rotate;
  void (RigidBody::*Rotate3)(const QUATERNION& rot) = &RigidBody::Rotate;
  void (RigidBody::*Rotate4)(const QUATERNION& rot, const VECTOR& pivot) = &RigidBody::Rotate;
  void (RigidBody::*Rotate5)(double phi, const VECTOR& dir) = &RigidBody::Rotate;
  void (RigidBody::*Rotate6)(double phi, const VECTOR& dir, const VECTOR& pivot) = &RigidBody::Rotate;

  void (RigidBody::*expt_Rotate_I_v1)(double, const VECTOR& dir)   = &RigidBody::Rotate_I;


  void (RigidBody::*expt_propagate_dlml_v1)(double t,double&) = &RigidBody::propagate_dlml;
  double (RigidBody::*expt_propagate_dlml_v2)(double t) = &RigidBody::propagate_dlml;




  class_<RigidBody>("RigidBody",init<>())
      .def(init<const RigidBody&>())
      .def("__copy__", &generic__copy__<RigidBody>)
      .def("__deepcopy__", &generic__deepcopy__<RigidBody>)
       // Topology
      .def_readwrite("rb_centers",&RigidBody::rb_centers)

      .def_readwrite("rb_mass",&RigidBody::rb_mass)
      .def_readwrite("rb_iM",&RigidBody::rb_iM)
      .def_readwrite("rb_cm",&RigidBody::rb_cm)
      .def_readwrite("rb_p",&RigidBody::rb_p)
      .def_readwrite("rb_v",&RigidBody::rb_v)
      .def_readwrite("rb_force",&RigidBody::rb_force)
      .def_readwrite("rb_I_I",&RigidBody::rb_I_I)
      .def_readwrite("rb_I_e",&RigidBody::rb_I_e)
      .def_readwrite("rb_invI_I",&RigidBody::rb_invI_I)
      .def_readwrite("rb_invI_e",&RigidBody::rb_invI_e)
      .def_readwrite("rb_A",&RigidBody::rb_A)
      .def_readwrite("rb_B",&RigidBody::rb_B)
      .def_readwrite("rb_C",&RigidBody::rb_C)

      .def_readwrite("rb_A_I_to_e",&RigidBody::rb_A_I_to_e)
      .def_readwrite("rb_A_I_to_e_T",&RigidBody::rb_A_I_to_e_T)
      .def_readwrite("rb_L",&RigidBody::rb_L)
      .def_readwrite("rb_p_r",&RigidBody::rb_p_r)
      .def_readwrite("rb_l_e",&RigidBody::rb_l_e)
      .def_readwrite("rb_w_e",&RigidBody::rb_w_e)
      .def_readwrite("rb_torque_e",&RigidBody::rb_torque_e)

      .def_readwrite("is_fixed_translation",&RigidBody::is_fixed_translation)
      .def_readwrite("is_fixed_rotation",&RigidBody::is_fixed_rotation)

      .def_readwrite("set_orientation_option", &RigidBody::set_orientation_option)
      .def_readwrite("is_set_orientation_option", &RigidBody::is_set_orientation_option)


      .def("set",&RigidBody::set)
      .def("show_info",&RigidBody::show_info)


      .def("init",init1)
      .def("init",init2)
      .def("set_mass",&RigidBody::set_mass)
      .def("set_position",&RigidBody::set_position)
      .def("set_momentum",&RigidBody::set_momentum)
      .def("set_velocity",&RigidBody::set_velocity)
      .def("set_force",&RigidBody::set_force)
      .def("set_torque",&RigidBody::set_torque)
      .def("set_forces_and_torques",&RigidBody::set_forces_and_torques)
      .def("set_inertia",&RigidBody::set_inertia)
      .def("set_orientation",set_orientation1)
      .def("set_orientation",set_orientation2)
      .def("set_orientation",set_orientation3)
      .def("set_angular_momentum",set_angular_momentum1)
      .def("set_angular_momentum",set_angular_momentum2)
      .def("set_angular_velocity",set_angular_velocity1)
      .def("set_angular_velocity",set_angular_velocity2)
      .def("set_quaternion_momentum",&RigidBody::set_quaternion_momentum)
      .def("scale_angular_",scale_angular_1)
      .def("scale_angular_",scale_angular_2)
      .def("scale_angular_",scale_angular_3)

      .def("scale_linear_",scale_linear_1)
      .def("scale_linear_",scale_linear_2)
      .def("scale_linear_",scale_linear_3)
      .def("scale_position",scale_position1)
      .def("scale_position",scale_position2)
      .def("scale_position",scale_position3)
      .def("shift_angular_momentum",&RigidBody::shift_angular_momentum)
      .def("shift_angular_velocity",&RigidBody::shift_angular_velocity)
      .def("shift_linear_momentum",&RigidBody::shift_linear_momentum)
      .def("shift_linear_velocity",&RigidBody::shift_linear_velocity)
      .def("shift_position",&RigidBody::shift_position)
      .def("fix_translation",&RigidBody::fix_translation)
      .def("fix_rotation",&RigidBody::fix_rotation)
      .def("unfix_translation",&RigidBody::unfix_translation)
      .def("unfix_rotation",&RigidBody::unfix_rotation)
      .def("apply_torque",&RigidBody::apply_torque)
      .def("apply_force",apply_force1)
      .def("apply_force",apply_force2)


      .def("get_Nf_t",&RigidBody::get_Nf_t, "returns the number of translational degrees of freedom for this RB")      
      .def("get_Nf_r",&RigidBody::get_Nf_r, "returns the number of rotational degrees of freedom for this RB")
      .def("get_center_in_global_frame",&RigidBody::get_center_in_global_frame)  
      .def("get_center_in_lab_frame",&RigidBody::get_center_in_lab_frame)
      .def("get_center_in_body_frame",&RigidBody::get_center_in_body_frame)
      .def("body_frame_to_lab_frame",&RigidBody::body_frame_to_lab_frame, "convert vector from body coordinate system to the lab(external) coordinate system")
      .def("lab_frame_to_body_frame",&RigidBody::lab_frame_to_body_frame)
      .def("ekin_rot",&RigidBody::ekin_rot,"returns rotational kinetic energy") 
      .def("ekin_tr",&RigidBody::ekin_tr,"returns translational kinetic energy")
      .def("Rotate_I_x",&RigidBody::Rotate_I_x,"rotate arond x axis in global(Cartesian) coordinate system")
      .def("Rotate_I_y",&RigidBody::Rotate_I_y,"rotate arond y axis in global(Cartesian) coordinate system")
      .def("Rotate_I_z",&RigidBody::Rotate_I_z,"rotate arond z axis in global(Cartesian) coordinate system")
      .def("Rotate_e_x",&RigidBody::Rotate_e_x,"rotate arond x axis in rigid body coordinate system")
      .def("Rotate_e_y",&RigidBody::Rotate_e_y,"rotate arond y axis in rigid body coordinate system")
      .def("Rotate_e_z",&RigidBody::Rotate_e_z,"rotate arond z axis in rigid body coordinate system")
      .def("Rotate",Rotate1,"arbitrary rotation of RB, rotation is defined by the transformation matrix")
      .def("Rotate",Rotate2,"arbitrary rotation of RB, rotation is defined by the transformation matrix")
      .def("Rotate",Rotate3,"arbitrary rotation of RB, rotation is defined by the quaternion")
      .def("Rotate",Rotate4,"arbitrary rotation of RB, rotation is defined by the quaternion")
      .def("Rotate",Rotate5,"")
      .def("Rotate",Rotate6,"")


      .def("Rotate_I",expt_Rotate_I_v1, 
           "rotate the RB around the specified vector in the external (moving) frame \
            for a magnitude defined by the provided argument (in radians)"
          )


      .def("initialize_exact_rb",&RigidBody::initialize_exact_rb,"initialization of the R van Zon integrator")
      .def("propagate_exact_rb",&RigidBody::propagate_exact_rb,"exact solution for free RB problem, based on Jacoby method as discussed by R van Zon")
      .def("propagate_dlml", expt_propagate_dlml_v1,"DLML = Dullweber, Leimkuhler, McLachlan")
      .def("propagate_dlml", expt_propagate_dlml_v2,"DLML = Dullweber, Leimkuhler, McLachlan")

      .def("propagate_no_squish",&RigidBody::propagate_no_squish,"Symplectic quaternion scheme of Miller and co-workers")
      .def("initialize_terec",&RigidBody::initialize_terec,"inialization of the Terec integrator, see JCTC paper by Akimov & Kolomeisky")
      .def("propagate_terec",&RigidBody::propagate_terec,"Terec integrator with matrices, see JCTC paper by Akimov & Kolomeisky")
      .def("propagate_qterec",&RigidBody::propagate_qterec,"Terec integrator with quaternions, see JCTC paper by Akimov & Kolomeisky")
      .def("propagate_kln",&RigidBody::propagate_kln,"4-rotation integrator of Kamberaj, Low, Neal")
      .def("propagate_omelyan",&RigidBody::propagate_omelyan,"Omelyan RB integrator")




  ;



//  class_<std::vector<RigidBody> >("RigidBodyList")
//      .def(vector_indexing_suite<std::vector<RigidBody> >())
//  ;

}



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygrigidbody){
#else
BOOST_PYTHON_MODULE(librigidbody){
#endif


//  export_Mathematics_objects();  // also register mmath python functions!
  export_RigidBody_objects();

}

}// namespace librigidbody

}// liblibra



