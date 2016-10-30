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
#include "libspecialfunctions.h"

using namespace boost::python;


/// libmmath namespace
namespace libmmath{

/// libspecialfunctions namespace
namespace libspecialfunctions{


void export_SpecialFunctions_objects(){

  boost::python::list (*expt_binomial_expansion)(int, int, double, double, int) = &binomial_expansion;

  // Now introduce normal functions:
  def("FAST_POW", FAST_POW);

  def("sinh_",sinh_);  // sinh(x)/x
  def("sin_", sin_);   // sin(x)/x
  def("ERF",ERF);      // error function
  def("ERFC",ERFC);    // complementary error function
  def("gamma_lower", gamma_lower);  // lower gamma function divided by the power
  def("Fn", Fn);
  def("gaussian_int", gaussian_int);  
  def("gaussian_norm2", gaussian_norm2);
  def("gaussian_norm1", gaussian_norm1);
  def("gaussian_normalization_factor", gaussian_normalization_factor);


  def("FACTORIAL",FACTORIAL); // n!
  def("DFACTORIAL",DFACTORIAL);  // n!!
  def("BINOM",BINOM);  // C_n^i
  def("binomial_expansion", expt_binomial_expansion);

  def("LEGENDRE",LEGENDRE);
  def("CHEBYSHEV",CHEBYSHEV);
  def("LAGUERRE",LAGUERRE);
  def("HERMITE",HERMITE); 

  def("Ellipe",Ellipe);
  def("Ellipe2",Ellipe2);
  def("Jacobi_Elliptic", Jacobi_Elliptic);
  def("Km", Km);
  def("Ellint",Ellint);

  def("randperm", randperm);

  MATRIX (*expt_exp__v1)(MATRIX&, double) = &exp_;
  MATRIX3x3 (*expt_exp__v2)(MATRIX3x3&, double) = &exp_;
  MATRIX (*expt_exp1__v1)(MATRIX&, double) = &exp1_;
  MATRIX3x3 (*expt_exp1__v2)(MATRIX3x3&, double) = &exp1_;

  def("exp_", expt_exp__v1);
  def("exp_", expt_exp__v2);
  def("exp1_", expt_exp1__v1);
  def("exp1_", expt_exp1__v2);


  void (*expt_MATRIX_TO_QUATERNION_v1)(MATRIX&,QUATERNION&) = &MATRIX_TO_QUATERNION;
  void (*expt_MATRIX_TO_QUATERNION_v2)(MATRIX3x3&,QUATERNION&) = &MATRIX_TO_QUATERNION;

  void (*expt_QUATERNION_TO_MATRIX_v1)(QUATERNION&,MATRIX&) = &QUATERNION_TO_MATRIX;
  void (*expt_QUATERNION_TO_MATRIX_v2)(QUATERNION&,MATRIX3x3&) = &QUATERNION_TO_MATRIX;

  def("MATRIX_TO_QUATERNION", expt_MATRIX_TO_QUATERNION_v1);
  def("MATRIX_TO_QUATERNION", expt_MATRIX_TO_QUATERNION_v2);
  def("QUATERNION_TO_MATRIX", expt_QUATERNION_TO_MATRIX_v1);
  def("QUATERNION_TO_MATRIX", expt_QUATERNION_TO_MATRIX_v2);

  void (*expt_solve_linsys_v1)(MATRIX&,MATRIX&, MATRIX&,double,int) = &solve_linsys;

  def("solve_linsys", expt_solve_linsys_v1); 



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygspecialfunctions){
#else
BOOST_PYTHON_MODULE(libspecialfunctions){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
  export_SpecialFunctions_objects();

}

}// libspecialfunctions
}// libmmath
