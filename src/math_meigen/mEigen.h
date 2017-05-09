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

#ifndef MMATH_EIGEN_H
#define MMATH_EIGEN_H

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace boost::python;
using namespace liblinalg; 

/// libmeigen namespace
namespace libmeigen{


///=========== Look in: mEigen_determinant.cpp =================
///< Several versions of the determinants
double det(MATRIX&);
double FullPivLU_det(MATRIX&);
complex<double> det(CMATRIX&);
complex<double> FullPivLU_det(CMATRIX& A);


///=========== Look in: mEigen_eigensolve.cpp ==================
///< Solving the generalized eigenvalue problem: H * C = S * C * E
void solve_eigen(MATRIX* H, MATRIX* S, MATRIX* E, MATRIX* C, int symm);
void solve_eigen(MATRIX& H, MATRIX& S, MATRIX& E, MATRIX& C, int symm);
void solve_eigen(MATRIX* H, MATRIX* S, CMATRIX* E, CMATRIX* C, int symm);
void solve_eigen(MATRIX& H, MATRIX& S, CMATRIX& E, CMATRIX& C, int symm);
void solve_eigen(CMATRIX* H, CMATRIX* S, CMATRIX* E, CMATRIX* C, int symm);
void solve_eigen(CMATRIX& H, CMATRIX& S, CMATRIX& E, CMATRIX& C, int symm);

///< Solving the eigenvalue problem: H * C = C * E
void solve_eigen(MATRIX* H, MATRIX* E, MATRIX* C, int symm);
void solve_eigen(MATRIX& H, MATRIX& E, MATRIX& C, int symm);
void solve_eigen(MATRIX* H, CMATRIX* E, CMATRIX* C, int symm);
void solve_eigen(MATRIX& H, CMATRIX& E, CMATRIX& C, int symm);
void solve_eigen(CMATRIX* H, CMATRIX* E, CMATRIX* C, int symm);
void solve_eigen(CMATRIX& H, CMATRIX& E, CMATRIX& C, int symm);


///=========== Look in: mEigen_inverse.cpp ==================
///< Invertability check
void FullPivLU_rank_invertible(MATRIX& A, int& rank, int& is_inver);
void FullPivLU_rank_invertible(CMATRIX& A, int& rank, int& is_inver);
boost::python::list FullPivLU_rank_invertible(MATRIX& A);
boost::python::list FullPivLU_rank_invertible(CMATRIX& A);

///< Matrix inversion via an LU decomposition
void FullPivLU_inverse(MATRIX& A, MATRIX& invA);
void FullPivLU_inverse(CMATRIX& A, CMATRIX& invA);

///< Matrix inversion via solving the eigenvalue problem
void inv_matrix(MATRIX& S, MATRIX& S_inv, double thresh);
void inv_matrix(MATRIX& S, MATRIX& S_inv);

void inv_matrix(CMATRIX& S, CMATRIX& S_inv, double thresh);
void inv_matrix(CMATRIX& S, CMATRIX& S_inv);


///=========== Look in: mEigen_matrix_functions.cpp ==================
///< Square root and inverse of a matrix
void sqrt_matrix(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half, double thresh);
void sqrt_matrix(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half);
void exp_matrix(CMATRIX& res, CMATRIX& S, complex<double> dt);


///=========== Look in: mEigen_decompositions.cpp ==================
///< LU decomposition
void FullPivLU_decomposition(MATRIX& A, MATRIX& P, MATRIX& L, MATRIX& U, MATRIX& Q);
void FullPivLU_decomposition(CMATRIX& A, CMATRIX& P, CMATRIX& L, CMATRIX& U, CMATRIX& Q);


///=========== Look in: mEigen_linsolve.cpp ==================
///< Solver for a system of linear equations (iterative schemes)
bool linsys_solver(const MATRIX& A, MATRIX& X, const MATRIX& B, const double NormThreshold);
void solve_linsys(MATRIX&,MATRIX&, MATRIX&,double,int);

}// namespace libmeigen
}// namespace liblibra


#endif
