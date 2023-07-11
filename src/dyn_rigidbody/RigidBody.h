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
  \file RigidBody.h
  \brief The file describes RigidBody class and related functions
    
*/

#ifndef RigidBody_H
#define RigidBody_H

#include "../io/libio.h"
#include "../math_linalg/liblinalg.h"
#include "../math_specialfunctions/libspecialfunctions.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;
using namespace libio;


/// librigidbody namespace
namespace librigidbody{

class RigidBody{
/** 
  \brief the RigidBody class to represent the coordinates of a rigid body

  See reference  "Rigid_Body_Mechancs"
  Defenitions of coordinate system names.   
   i - laboratory coordinate system          
   I - movable    coordinate system (parallel to i)
   e - internal   coordinate system (rigidly connected with the body)

*/


public:

  vector<VECTOR> rb_centers;  int rb_centers_size;

  //------ Basic dynamical variables and constants --------
  // Total mass
  double rb_mass;             int is_rb_mass;  ///< the mass of the RB and the flag showing its status (defined or not)
  double rb_iM;               int is_rb_iM;    ///< the inverse mass of the RB and the flag showing its status (defined or not)

  // Center of mass
  VECTOR rb_cm;               int is_rb_cm;    ///< the position of the center of mass of the RB and the flag showing its status (defined or not)

  // Linear momentum and velocity of the center of mass
  VECTOR rb_p;                int is_rb_p;     ///< the RB momentum and the flag showing its status (defined or not)
  VECTOR rb_v;                int is_rb_v;     ///< the RB velocity and the flag showing its status (defined or not)

  // Total force acting on the center of mass 
  VECTOR rb_force;            int is_rb_force; ///< the RB force and the flag showing its status (defined or not)

  // Inertia moments
  MATRIX3x3 rb_I_I;           int is_rb_I_I;   ///< the RB intertia tensor in movable coordinate system and the flag showing its status (defined or not)
  MATRIX3x3 rb_I_e;           int is_rb_I_e;   ///< the RB intertia tensor (diagonal) in internal coordinate system and the flag showing its status (defined or not)

  // Inverse inertia moments
  MATRIX3x3 rb_invI_I;        int is_rb_invI_I;///< the inverse RB intertia tensor in movable coordinate system and the flag showing its status (defined or not)            
  MATRIX3x3 rb_invI_e;        int is_rb_invI_e;///< the inverse RB intertia tensor (diagonal) in internal coordinate system and the flag showing its status (defined or not)

  // Rotational constants
  double rb_A;                int is_rb_A;  ///< the RB rotational constant the flag showing its status (defined or not)            
  double rb_B;                int is_rb_B;  ///< the RB rotational constant the flag showing its status (defined or not)            
  double rb_C;                int is_rb_C;  ///< the RB rotational constant the flag showing its status (defined or not)            

  //  Attitude matrix and its transpose
  MATRIX3x3 rb_A_I_to_e;      int is_rb_A_I_to_e;   ///< the RB transformation matrix: from the movable to internal and the flag showing its status (defined or not)            
                                                    ///< this is essentially the orientation variable (matrix)
  MATRIX3x3 rb_A_I_to_e_T;    int is_rb_A_I_to_e_T; ///< the transpose of the RB transformation matrix: from the movable to internal and the flag showing its status (defined or not)  
                                                    ///< this is essentially inverse of the rb_A_I_to_e matrix (inverse rotation)

  // Orientation quaternion and its conjugate momentum 
  QUATERNION rb_L;            int is_rb_L;   ///< The orientation variables in the quaternion form and the the flag showing its status (defined or not)   
  QUATERNION rb_p_r;          int is_rb_p_r; ///< The momentum conjugate to quaternion variables and the flag showing its status (defined or not)            

  // Angular momentum and angular velocity (in body frame) 
  VECTOR rb_l_e;              int is_rb_l_e; ///< The angular momentum in the body frame and the the flag showing its status (defined or not)   
  VECTOR rb_w_e;              int is_rb_w_e; ///< The angular velocity in the body frame and the the flag showing its status (defined or not)   

  // Total torque (in body-fixed coordinate system)
  VECTOR rb_torque_e;         int is_rb_torque_e; ///< The torque (angular generalized force) in the body frame and the the flag showing its status (defined or not)   

  // Fixation (constraint) flags
  int is_fixed_translation;  ///< The flag showing whether the translation of this RB is forbidden (if = 1) or not (if = 0)
  int is_fixed_rotation;     ///< The flag showing whether the rotation of this RB is forbidden (if = 1) or not (if = 0)

  // Orientation update option - controls how set_orientation function works
  int set_orientation_option;  ///< if 1 - use the previous approach that is we convert the updated orientation matrix to the quaternions for the later use
                               ///< if 0 - forget about updating quaternions
  int is_set_orientation_option; /// the flag of the setting up


private:

  //---------- Auxiliary internal variables ------------

  // Constants      
  double MACHPREC;      ///< "Machine precision". Defined to 1e-15 during the initialization. 
  int IntN;             ///< Elliptic integrals evaluation parameter (???): 10 
  double tol;           ///< Tolerance for some (???) calculations: 1e-15
  double IEPS;          ///< another accurcy parameter (???): 1e-3
  double BIG;           ///< When to set rotational constants to zero. If they are larger than this: 1e+10
  double MAX_NO;        ///< Stopping criteria for Jacobi orthogonalization method: 1e-15
  double minDet;        ///< Minimal value of determinant to treat matrix as non-degenerate: 1e-10
    
  MATRIX3x3 U[7];       ///< Permutation matrices
  int permutindx;       ///< Index of permutation matrix necessary for Jacobi ordering
  int invpermutindx;    ///< Index of permutation matrix inverse for that with permutindx;
  int orderflag;
  int SERIES_EXPANSION; ///< For TEREC algorithm, the number of terms in expansion. Default: 10

  // Internal variables used for exact solution of free rigid-body problem

  VECTOR top_l;         ///< initial angular momentum
  VECTOR top_w;         ///< initial angular velocity
  VECTOR top_wm;        ///< angular velocity amplitudes
  MATRIX3x3 Iint;       ///< Tensor of inertia in principal axes (diagonal)
  MATRIX3x3 invI;       ///< inversion of the inertia moment tensor (diagonal)

  double I1, I2, I3;    ///< Jacobi-ordered principal inertia moments
  double m;             ///< elliptic parameter
  double eps;
  double K, Kcompl;     ///< quater period and complimentary quater period K = F(1|m), Kcompl = F(1|1-m)
  double q;             ///< nome q = exp(-pi*Kcompl/K)
  double A1,A2;
  int NT;
  double* cr;           ///< Real expansion coefficients
  double* ci;           ///< Imaginary expansion coefficient
  double wp;            ///< precession frequency

  ///< Internal variables for no_squish integrator
  MATRIX *P1, *P2, *P3;

  ///< Internal variables for TEREC integrator
  double* Coeffs;

  //--------- Auxiliary internal functions -------------
  // Defined in RigidBody.cpp
  void init_permutations(); ///< Initializes permutation matrices U
  void init_variables(int);///< Initializes dynamical variables
  void copy_content(const RigidBody&); ///< Copies the content which is defined

  // The following functions are the sub-routines of 
  // public method init(int,double*,VECTOR*) which must
  // be performed in this order!!!
  // Defined in RigidBody_methods.cpp
  void calc_mass(int, double*);
  void calc_center_of_mass(int, double*, VECTOR*);
  void calc_inertia_tensors(int, double*, VECTOR*);
  void calc_orientations(int, double*, VECTOR*);
  void calc_inverse_tensors(int);
  void calc_rot_constants();
  void calc_rb_centers(int,VECTOR*);

  // Aux function for no_squish  
  void init_S_matrix(MATRIX&);       // <-- Defined in RigidBody_methods1.cpp
  void rotate_no_squish(int,double); // <-- Defined in RigidBody_methods5_3.cpp
  


public:

  // Basic class methods
  // Defined in RigidBody.cpp
  RigidBody();                  ///< constructor
  RigidBody(const RigidBody&);  ///< copy-constructor
 ~RigidBody();                  ///< destructor

  RigidBody& operator=(const RigidBody&); ///< assignment operator
  void show_info();
  void set(object);
  void save(boost::property_tree::ptree& pt,std::string path);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);


  // Class interface - initialization and setters
  // Defined in RigidBody_methods1.cpp
  int init(int,double*,VECTOR*);
  int init(int sz,boost::python::list masses,boost::python::list positions);

  int set_mass(const double&);
  int set_position(const VECTOR&);
  int set_momentum(const VECTOR&);
  int set_velocity(const VECTOR&);
  int set_force(const VECTOR&);
  int set_torque(const VECTOR&);
  int set_forces_and_torques(int,VECTOR*,VECTOR*);
  int set_inertia(const MATRIX3x3&); 
  int set_orientation(const MATRIX3x3&);
  int set_orientation(const VECTOR&, const VECTOR&, const VECTOR&);
  int set_orientation(const QUATERNION&);
  int set_angular_momentum(const VECTOR&);
  int set_angular_momentum(const double&, const double&, const double&);
  int set_angular_velocity(const VECTOR&);
  int set_angular_velocity(const double&, const double&, const double&);
  int set_quaternion_momentum(const QUATERNION&);

  // Class interface - modifiers, internal variables transformations/propagation
  // Defined in RigidBody_methods2.cpp
  int scale_angular_(double);
  int scale_angular_(double,double,double);
  int scale_angular_(const MATRIX3x3&);
  int scale_linear_(double);
  int scale_linear_(double,double,double);
  int scale_linear_(const MATRIX3x3&);
  int scale_position(double);
  int scale_position(double,double,double);
  int scale_position(const MATRIX3x3&);
  int shift_angular_momentum(const VECTOR&);
  int shift_angular_velocity(const VECTOR&);
  int shift_linear_momentum(const VECTOR&);
  int shift_linear_velocity(const VECTOR&);
  int shift_position(const VECTOR&);
  int fix_translation();
  int fix_rotation();
  int unfix_translation();
  int unfix_rotation();
  int apply_torque(double);
  int apply_force(double);
  int apply_force(MATRIX3x3&);


  // Class interface - getters, properties, external variable transformations
  // Defined in RigidBody_methods3.cpp
  int get_Nf_t();
  int get_Nf_r();
  VECTOR get_center_in_global_frame(int);
  VECTOR get_center_in_lab_frame(int);
  VECTOR get_center_in_body_frame(int);
  void body_frame_to_lab_frame(const VECTOR&,VECTOR&);
  void lab_frame_to_body_frame(const VECTOR&,VECTOR&);
  double ekin_rot();
  double ekin_tr();

  // Rotations
  // Defined in RigidBody_methods4.cpp
  void Rotate_I_x(double);
  void Rotate_I_y(double);
  void Rotate_I_z(double);

  void Rotate_e_x(double);
  void Rotate_e_y(double);
  void Rotate_e_z(double);

  void Rotate(const MATRIX3x3& rot);
  void Rotate(const MATRIX3x3& rot, const VECTOR& pivot);
  void Rotate(const QUATERNION& rot);
  void Rotate(const QUATERNION& quat, const VECTOR& pivot);
  void Rotate(double phi, const VECTOR& dir);
  void Rotate(double phi, const VECTOR& dir, const VECTOR& pivot);


  void Rotate_I(double, const VECTOR&);


  // Integrators
  // Defined in RigidBody_methods5_*.cpp
  // Defined in RigidBody_methods5_1.cpp
  void initialize_exact_rb();
  void propagate_exact_rb(double t);
  // Defined in RigidBody_methods5_2.cpp
  void propagate_dlml(double t,double&);
  double propagate_dlml(double t);
  // Defined in RigidBody_methods5_3.cpp
  void propagate_no_squish(double t);
  // Defined in RigidBody_methods5_4.cpp
  void initialize_terec(int);
  void propagate_terec(double);
  void propagate_qterec(double);
  // Defined in RigidBody_methods5_5.cpp
  void propagate_kln(double);
  // Defined in RigidBody_methods5_6.cpp
  void propagate_omelyan(double);


};

void save(boost::property_tree::ptree& pt,std::string path,vector<RigidBody>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<RigidBody>& vt,int& status);


}// namespace librigidbody

}// liblibra

#endif // RigidBody_H
