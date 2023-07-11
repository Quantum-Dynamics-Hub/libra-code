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
  \file libmolint.cpp
  \brief The file describes Python export function
    
*/

#define BOOST_PYTHON_MAX_ARITY 30

//#if defined(USING_PCH)
//#include "../pch.h"
//#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
//#endif 

#include "libmolint.h"

#include "../math_specialfunctions/libspecialfunctions.h"
#include "../math_linalg/liblinalg.h"


/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;
using namespace boost::python;


/// libmolint namespace
namespace libmolint{



void export_molint_objects(){
/** 
  \brief Exporter of libmolint classes and functions

*/


  // Basic overlaps
  // 1D Gaussians
  double (*expt_gaussian_overlap_ref)(int,double,double, int,double,double) = &gaussian_overlap_ref;
  double (*expt_gaussian_overlap_1D_v1)(int,double,double, int,double,double) = &gaussian_overlap;
  double (*expt_gaussian_overlap_1D_v2)(int,double,double, int,double,double, int) = &gaussian_overlap;
  boost::python::list (*expt_gaussian_overlap_1D_v3)(int,double,double, int,double,double, int, int) = &gaussian_overlap;


  // 3D Gaussians
  double (*expt_gaussian_overlap_3D_v1)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&  ) = &gaussian_overlap;
  double (*expt_gaussian_overlap_3D_v2)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&, int ) = &gaussian_overlap;
  boost::python::list (*expt_gaussian_overlap_3D_v3)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&, int, int ) = &gaussian_overlap;


  // Moments 
  // 1D Moments
  double (*expt_gaussian_moment_ref)(int,double,double, int,double, double, int,double, double ) = &gaussian_moment;
  double (*expt_gaussian_moment_1D_v1)(int,double,double, int,double,double, int,double,double ) = &gaussian_moment;
  double (*expt_gaussian_moment_1D_v2)(int,double,double, int,double,double, int,double,double, int ) = &gaussian_moment;
  boost::python::list (*expt_gaussian_moment_1D_v3)(int,double,double, int,double,double, int,double,double, int, int ) = &gaussian_moment;
 

  // 3D Moments
  double (*expt_gaussian_moment_3D_v1)(int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&          ) = &gaussian_moment;
  double (*expt_gaussian_moment_3D_v2)(int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&, int     ) = &gaussian_moment;
  boost::python::list (*expt_gaussian_moment_3D_v3)(int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&, int, int) = &gaussian_moment;



  // Pseudopotentials
  // 3D pseudopotentials
  double (*expt_pseudopot02_v1)(double,double,double,const VECTOR&,
                                int,int,int,double,const VECTOR&,
                                int,int,int,double,const VECTOR&                 ) = &pseudopot02;
  boost::python::list (*expt_pseudopot02_v2)(double, double, double, const VECTOR&,
                                int,int,int,double, const VECTOR&,
                                int,int,int,double, const VECTOR&,  int, int      ) = &pseudopot02;

  // Multipole moments
  // 1D - Transition dipole moments
  double (*expt_transition_dipole_moment_1D_v1)
  ( int nxa,double alp_a, double Xa, 
    int nxb,double alp_b, double Xb
  ) = &transition_dipole_moment;
  double (*expt_transition_dipole_moment_1D_v2)
  ( int nxa,double alp_a, double Xa, 
    int nxb,double alp_b, double Xb,
    int is_normalize
  ) = &transition_dipole_moment;
  boost::python::list (*expt_transition_dipole_moment_1D_v3)
  ( int nxa,double alp_a, double Xa, 
    int nxb,double alp_b, double Xb,
    int is_normalize, int is_derivs
  ) = &transition_dipole_moment;

  // 3D - Transition dipole moments
  VECTOR (*expt_transition_dipole_moment_3D_v1)
  ( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
    int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb
  ) = &transition_dipole_moment;
  VECTOR (*expt_transition_dipole_moment_3D_v2)
  ( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
    int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
    int is_normalize
  ) = &transition_dipole_moment;
  boost::python::list (*expt_transition_dipole_moment_3D_v3)
  ( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
    int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
    int is_normalize,int is_derivs
  ) = &transition_dipole_moment;




  // Kinetic energy integrals (KEI)
  // 1D KEI
  double (*expt_kinetic_integral_1D_v1)(int,double, double, int, double, double)      = &kinetic_integral;
  double (*expt_kinetic_integral_1D_v2)(int,double, double, int, double, double, int) = &kinetic_integral;
  boost::python::list (*expt_kinetic_integral_1D_v3)(int,double, double, int, double, double, int, int) = &kinetic_integral;

  // 3D KEI
  double (*expt_kinetic_integral_3D_v1)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&          ) = &kinetic_integral;
  double (*expt_kinetic_integral_3D_v2)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&, int     ) = &kinetic_integral;
  boost::python::list (*expt_kinetic_integral_3D_v3)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&, int, int) = &kinetic_integral;


  // Nuclear attraction integrals (3D)
  double (*expt_nuclear_attraction_integral_v1)(int,int, int, double, VECTOR&,
                        int,int,int,double,VECTOR&,VECTOR&                  ) = &nuclear_attraction_integral;
  double (*expt_nuclear_attraction_integral_v2)(int,int, int, double, VECTOR&,
                        int,int,int,double,VECTOR&,VECTOR&, int             ) = &nuclear_attraction_integral;
  boost::python::list (*expt_nuclear_attraction_integral_v3)(int,int, int, double, VECTOR&,
                        int,int,int,double,VECTOR&,VECTOR&, int, int        ) = &nuclear_attraction_integral;



  // Electron repulsion integrals
  double (*expt_electron_repulsion_integral_v1)
  (
   int,int,int,double,VECTOR&,   int,int,int,double,VECTOR&,
   int,int,int,double,VECTOR&,   int,int,int,double,VECTOR&
  ) = &electron_repulsion_integral;

  double (*expt_electron_repulsion_integral_v2)
  (
   int,int,int,double,VECTOR&,   int,int,int,double,VECTOR&,
   int,int,int,double,VECTOR&,   int,int,int,double,VECTOR&, int
  ) = &electron_repulsion_integral;

  boost::python::list (*expt_electron_repulsion_integral_v3)
  (
   int,int,int,double,VECTOR&,   int,int,int,double,VECTOR&,
   int,int,int,double,VECTOR&,   int,int,int,double,VECTOR&, int, int
  ) = &electron_repulsion_integral;



  // Derivative couplings
  double (*expt_derivative_coupling_integral_1D_v1)
  ( int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb
  ) = &derivative_coupling_integral;

  double (*expt_derivative_coupling_integral_1D_v2)
  ( int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
    int is_normalize
  ) = &derivative_coupling_integral;

  boost::python::list (*expt_derivative_coupling_integral_1D_v3)
  ( int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
    int is_normalize, int is_derivs 
  ) = &derivative_coupling_integral;


  VECTOR (*expt_derivative_coupling_integral_3D_v1)
  ( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
    int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
  ) = &derivative_coupling_integral;

  VECTOR (*expt_derivative_coupling_integral_3D_v2)
  ( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
    int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
    int is_normalize
  ) = &derivative_coupling_integral;

  boost::python::list (*expt_derivative_coupling_integral_3D_v3)
  ( int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
    int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
    int is_normalize, int is_derivs
  ) = &derivative_coupling_integral;


  // Approximate and misc integrals
  double (*expt_Coulomb_Integral_v1)
  (  double R, int n_i, double Jii, double ksi_i, std::string type_i, double q_i,
               int n_j, double Jjj, double ksi_j, std::string type_j, double q_j,
               double epsilon, int mode
  ) = &Coulomb_Integral;








  // ============ Now export functions =============
  // ==== Overlaps ====
  def("gaussian_overlap_ref", expt_gaussian_overlap_ref);
  def("gaussian_overlap", expt_gaussian_overlap_1D_v1);
  def("gaussian_overlap", expt_gaussian_overlap_1D_v2);
  def("gaussian_overlap", expt_gaussian_overlap_1D_v3);

  def("gaussian_overlap", expt_gaussian_overlap_3D_v1);
  def("gaussian_overlap", expt_gaussian_overlap_3D_v2);
  def("gaussian_overlap", expt_gaussian_overlap_3D_v3);

  def("A_coefficient_general", A_coefficient_general); 
  def("generate_coefficients", generate_coefficients);

  def("sto_norm", sto_norm);
  def("sto_overlap", sto_overlap);
  def("sto_overlap_fast", sto_overlap_fast);


  // ==== Moments ====                                     
  def("gaussian_moment_ref", expt_gaussian_moment_ref);
  def("gaussian_moment", expt_gaussian_moment_1D_v1);
  def("gaussian_moment", expt_gaussian_moment_1D_v2);
  def("gaussian_moment", expt_gaussian_moment_1D_v3);

  def("gaussian_moment", expt_gaussian_moment_3D_v1);
  def("gaussian_moment", expt_gaussian_moment_3D_v2);
  def("gaussian_moment", expt_gaussian_moment_3D_v3);


  // ==== Pseudopotentials ====
  def("pseudopot02", expt_pseudopot02_v1);
  def("pseudopot02", expt_pseudopot02_v2);


  // ==== Multipole moments ====
  def("transition_dipole_moment", expt_transition_dipole_moment_1D_v1);
  def("transition_dipole_moment", expt_transition_dipole_moment_1D_v2);
  def("transition_dipole_moment", expt_transition_dipole_moment_1D_v3);

  def("transition_dipole_moment", expt_transition_dipole_moment_3D_v1);
  def("transition_dipole_moment", expt_transition_dipole_moment_3D_v2);
  def("transition_dipole_moment", expt_transition_dipole_moment_3D_v3);




  // ==== Kinetic energy integrals ====
  def("kinetic_integral", expt_kinetic_integral_1D_v1);
  def("kinetic_integral", expt_kinetic_integral_1D_v2);
  def("kinetic_integral", expt_kinetic_integral_1D_v3);

  def("kinetic_integral", expt_kinetic_integral_3D_v1);
  def("kinetic_integral", expt_kinetic_integral_3D_v2);
  def("kinetic_integral", expt_kinetic_integral_3D_v3);


  // ==== Nuclear attraction integral ====
  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_v1);
  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_v2);
  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_v3);


  // ==== Electron repulsion integral ====
  def("electron_repulsion_integral", expt_electron_repulsion_integral_v1);
  def("electron_repulsion_integral", expt_electron_repulsion_integral_v2);
  def("electron_repulsion_integral", expt_electron_repulsion_integral_v3);


  // ==== Derivative couplings =====
  def("derivative_coupling_integral", expt_derivative_coupling_integral_1D_v1);
  def("derivative_coupling_integral", expt_derivative_coupling_integral_1D_v2);
  def("derivative_coupling_integral", expt_derivative_coupling_integral_1D_v3);

  def("derivative_coupling_integral", expt_derivative_coupling_integral_3D_v1);
  def("derivative_coupling_integral", expt_derivative_coupling_integral_3D_v2);
  def("derivative_coupling_integral", expt_derivative_coupling_integral_3D_v3);


  // ==== Approximate integrals ====
  def("Coulomb_Integral", expt_Coulomb_Integral_v1);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmolint){
#else
BOOST_PYTHON_MODULE(libmolint){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_molint_objects();

}


}// namespace libmolint
}// namespace liblibra

