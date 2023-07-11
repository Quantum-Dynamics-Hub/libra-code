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

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

#include "libmeigen.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;


/// libmeigen namespace
namespace libmeigen{


void export_mEigen_objects(){


//boost::python::list (*expt_rotate1)(double,double,double) = &expt_rotate;
//double (*expt_shift1)(double, double) = &expt_shift;
//double (*expt_scale1)(double, double) = &expt_scale;
//  def("rotate", expt_rotate1);
//  def("shift", expt_shift1);
//  def("scale", expt_scale1);

  double (*expt_det_v1)(MATRIX&) = &det;
  complex<double> (*expt_det_v2)(CMATRIX&) = &det;

  double (*expt_FullPivLU_det_v1)(MATRIX&) = &FullPivLU_det;
  complex<double> (*expt_FullPivLU_det_v2)(CMATRIX&) = &FullPivLU_det;

  def("det", expt_det_v1);
  def("det", expt_det_v2);

  def("FullPivLU_det", expt_FullPivLU_det_v1);
  def("FullPivLU_det", expt_FullPivLU_det_v2);


  void (*expt_solve_eigen_v1)(CMATRIX& H, CMATRIX& S, CMATRIX& E, CMATRIX& C, int symm) = &solve_eigen;
  void (*expt_solve_eigen_v2)(MATRIX& H, MATRIX& S, CMATRIX& E, CMATRIX& C, int symm) = &solve_eigen;
  void (*expt_solve_eigen_v3)(MATRIX& H, MATRIX& S, MATRIX& E, MATRIX& C, int symm) = &solve_eigen;
  def("solve_eigen", expt_solve_eigen_v1);
  def("solve_eigen", expt_solve_eigen_v2);
  def("solve_eigen", expt_solve_eigen_v3);

  void (*expt_solve_eigen_v4)(CMATRIX& H, CMATRIX& E, CMATRIX& C, int symm) = &solve_eigen;
  void (*expt_solve_eigen_v5)(MATRIX& H, CMATRIX& E, CMATRIX& C, int symm) = &solve_eigen;
  void (*expt_solve_eigen_v6)(MATRIX& H, MATRIX& E, MATRIX& C, int symm) = &solve_eigen;
  def("solve_eigen", expt_solve_eigen_v4);
  def("solve_eigen", expt_solve_eigen_v5);
  def("solve_eigen", expt_solve_eigen_v6);

  void (*expt_solve_eigen_nosort_v1)(CMATRIX& H, CMATRIX& E, CMATRIX& C, int symm) = &solve_eigen_nosort;
  void (*expt_solve_eigen_nosort_v2)(MATRIX& H, CMATRIX& E, CMATRIX& C, int symm) = &solve_eigen_nosort;
  void (*expt_solve_eigen_nosort_v3)(MATRIX& H, MATRIX& E, MATRIX& C, int symm) = &solve_eigen_nosort;
  def("solve_eigen_nosort", expt_solve_eigen_nosort_v1);
  def("solve_eigen_nosort", expt_solve_eigen_nosort_v2);
  def("solve_eigen_nosort", expt_solve_eigen_nosort_v3);





  void (*expt_sqrt_matrix_v1)(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half, double thresh, int do_phase_correction) = &sqrt_matrix;
  void (*expt_sqrt_matrix_v2)(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half, double thresh) = &sqrt_matrix;
  void (*expt_sqrt_matrix_v3)(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half) = &sqrt_matrix;

  def("sqrt_matrix", expt_sqrt_matrix_v1);
  def("sqrt_matrix", expt_sqrt_matrix_v2);
  def("sqrt_matrix", expt_sqrt_matrix_v3);


  void (*expt_inv_matrix_v1)(CMATRIX& S, CMATRIX& S_inv, double thresh, int do_phase_correction) = &inv_matrix;
  void (*expt_inv_matrix_v2)(CMATRIX& S, CMATRIX& S_inv, double thresh) = &inv_matrix;
  void (*expt_inv_matrix_v3)(CMATRIX& S, CMATRIX& S_inv) = &inv_matrix;
  void (*expt_inv_matrix_v4)(MATRIX& S, MATRIX& S_inv, double thresh, int do_phase_correction) = &inv_matrix;
  void (*expt_inv_matrix_v5)(MATRIX& S, MATRIX& S_inv, double thresh) = &inv_matrix;
  void (*expt_inv_matrix_v6)(MATRIX& S, MATRIX& S_inv) = &inv_matrix;

  def("inv_matrix", expt_inv_matrix_v1);
  def("inv_matrix", expt_inv_matrix_v2);
  def("inv_matrix", expt_inv_matrix_v3);
  def("inv_matrix", expt_inv_matrix_v4);
  def("inv_matrix", expt_inv_matrix_v5);
  def("inv_matrix", expt_inv_matrix_v6);


  void (*expt_exp_matrix_v1)(CMATRIX& res, CMATRIX& S, complex<double> dt, int do_phase_correction) = &exp_matrix;
  void (*expt_exp_matrix_v2)(CMATRIX& res, CMATRIX& S, complex<double> dt) = &exp_matrix;

  def("exp_matrix", expt_exp_matrix_v1);
  def("exp_matrix", expt_exp_matrix_v2);


  void (*expt_FullPivLU_rank_invertible_v1)(MATRIX& A, int& rank, int& is_inver) = &FullPivLU_rank_invertible;
  void (*expt_FullPivLU_rank_invertible_v2)(CMATRIX& A, int& rank, int& is_inver) = &FullPivLU_rank_invertible;
  boost::python::list (*expt_FullPivLU_rank_invertible_v3)(MATRIX& A) = &FullPivLU_rank_invertible;
  boost::python::list (*expt_FullPivLU_rank_invertible_v4)(CMATRIX& A) = &FullPivLU_rank_invertible;

  def("FullPivLU_rank_invertible", expt_FullPivLU_rank_invertible_v1);
  def("FullPivLU_rank_invertible", expt_FullPivLU_rank_invertible_v2);
  def("FullPivLU_rank_invertible", expt_FullPivLU_rank_invertible_v3);
  def("FullPivLU_rank_invertible", expt_FullPivLU_rank_invertible_v4);


  void (*expt_FullPivLU_decomposition_v1)(MATRIX& A, MATRIX& P, MATRIX& L, MATRIX& U, MATRIX& Q) = &FullPivLU_decomposition;
  void (*expt_FullPivLU_decomposition_v2)(CMATRIX& A, CMATRIX& P, CMATRIX& L, CMATRIX& U, CMATRIX& Q) = &FullPivLU_decomposition;

  def("FullPivLU_decomposition", expt_FullPivLU_decomposition_v1);
  def("FullPivLU_decomposition", expt_FullPivLU_decomposition_v2);


  void (*expt_FullPivLU_inverse_v1)(MATRIX& A, MATRIX& invA) = &FullPivLU_inverse;
  void (*expt_FullPivLU_inverse_v2)(CMATRIX& A, CMATRIX& invA) = &FullPivLU_inverse;

  def("FullPivLU_inverse", expt_FullPivLU_inverse_v1);
  def("FullPivLU_inverse", expt_FullPivLU_inverse_v2);



  void (*expt_JacobiSVD_decomposition_v1)(CMATRIX& A, CMATRIX& U, CMATRIX& S, CMATRIX& V) = &JacobiSVD_decomposition;
  void (*expt_BDCSVD_decomposition_v1)(CMATRIX& A, CMATRIX& U, CMATRIX& S, CMATRIX& V) = &BDCSVD_decomposition;

  def("JacobiSVD_decomposition", expt_JacobiSVD_decomposition_v1);
  def("BDCSVD_decomposition", expt_BDCSVD_decomposition_v1);


  bool (*expt_linsys_solver_v1)(const MATRIX& A, MATRIX& X, const MATRIX& B, const double NormThreshold) = &linsys_solver;
  def("linsys_solver", expt_linsys_solver_v1);

  void (*expt_solve_linsys_v1)(MATRIX& C,MATRIX& D, MATRIX& X,double eps, int maxiter) = &solve_linsys;
  def("solve_linsys", expt_solve_linsys_v1);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmeigen){
#else
BOOST_PYTHON_MODULE(libmeigen){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_mEigen_objects();

}


}//namespace libmeigen
}//namespace liblibra




