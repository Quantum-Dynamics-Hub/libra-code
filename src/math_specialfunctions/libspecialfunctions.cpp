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

#include "libspecialfunctions.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;

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

  void (*expt_sample_v1)(MATRIX& x, MATRIX& mean_x, MATRIX& sigma_x, Random& rnd) = &sample;
  def("sample", expt_sample_v1);

  def("set_random_state", set_random_state);

  MATRIX3x3 (*expt_exp__v1)(MATRIX3x3&, double) = &exp_;
  MATRIX (*expt_exp__v2)(MATRIX&, double) = &exp_;
  CMATRIX (*expt_exp__v3)(CMATRIX&, complex<double>) = &exp_;

  def("exp_", expt_exp__v1);
  def("exp_", expt_exp__v2);
  def("exp_", expt_exp__v3);


  CMATRIX (*expt_exp_2_v1)(CMATRIX& x, complex<double> dt, int nterms, double max_tol) = &exp_2;
  CMATRIX (*expt_exp_2_v2)(CMATRIX& x, complex<double> dt, int nterms) = &exp_2;
  CMATRIX (*expt_exp_2_v3)(CMATRIX& x, complex<double> dt) = &exp_2;

  MATRIX (*expt_exp_2_v4)(MATRIX& x, double dt, int nterms, double max_tol) = &exp_2;
  MATRIX (*expt_exp_2_v5)(MATRIX& x, double dt, int nterms) = &exp_2;
  MATRIX (*expt_exp_2_v6)(MATRIX& x, double dt) = &exp_2;

  def("exp_2", expt_exp_2_v1);
  def("exp_2", expt_exp_2_v2);
  def("exp_2", expt_exp_2_v3);
  def("exp_2", expt_exp_2_v4);
  def("exp_2", expt_exp_2_v5);
  def("exp_2", expt_exp_2_v6);


  MATRIX (*expt_exp1__v1)(MATRIX&, double) = &exp1_;
  MATRIX3x3 (*expt_exp1__v2)(MATRIX3x3&, double) = &exp1_;

  def("exp1_", expt_exp1__v1);
  def("exp1_", expt_exp1__v2);

  boost::python::list (*expt_merge_sort_v1)(boost::python::list) = &merge_sort;
  int (*expt_merge_sort_v2)(vector<double>&, vector<double>&) = &merge_sort;
  def("merge_sort", expt_merge_sort_v1); 
  def("merge_sort", expt_merge_sort_v2); 


  MATRIX (*expt_mean_v1)(MATRIX& X) = &mean;
  CMATRIX (*expt_mean_v2)(CMATRIX& X) = &mean;
  def("mean", expt_mean_v1); 
  def("mean", expt_mean_v2); 

  MATRIX (*expt_deviation_v1)(MATRIX& X) = &deviation;
  CMATRIX (*expt_deviation_v2)(CMATRIX& X) = &deviation;
  def("deviation", expt_deviation_v1); 
  def("deviation", expt_deviation_v2); 

  MATRIX (*expt_variance_v1)(MATRIX& X, int opt) = &variance;
  MATRIX (*expt_variance_v2)(CMATRIX& X, int opt) = &variance;
  def("variance", expt_variance_v1);
  def("variance", expt_variance_v2);

  MATRIX (*expt_std_dev_v1)(MATRIX& X, int opt) = &std_dev;
  MATRIX (*expt_std_dev_v2)(CMATRIX& X, int opt) = &std_dev;
  def("std_dev", expt_std_dev_v1);
  def("std_dev", expt_std_dev_v2);


  MATRIX (*expt_covariance_v1)(MATRIX& X) = &covariance;
  MATRIX (*expt_covariance_v2)(MATRIX& X, MATRIX& Y) = &covariance;
  CMATRIX (*expt_covariance_v3)(CMATRIX& X) = &covariance;
  CMATRIX (*expt_covariance_v4)(CMATRIX& X, CMATRIX& Y) = &covariance;
  def("covariance", expt_covariance_v1); 
  def("covariance", expt_covariance_v2); 
  def("covariance", expt_covariance_v3); 
  def("covariance", expt_covariance_v4); 

  def("permutations_reiteration", permutations_reiteration);
  def("compute_all_permutations", compute_all_permutations);


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
}// liblibra
