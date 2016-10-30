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
using namespace boost::python;


#include "../linalg/liblinalg.h"

/// libmmath namespace
namespace libmmath{

using namespace liblinalg; 

/// libmeigen namespace
namespace libmeigen{


double det(MATRIX&);
double FullPivLU_det(MATRIX&);
complex<double> det(CMATRIX&);
complex<double> FullPivLU_det(CMATRIX& A);

void solve_eigen(int, MATRIX*, MATRIX*, MATRIX*, MATRIX*);
void solve_eigen(int, MATRIX&, MATRIX&, MATRIX&, MATRIX&);
void solve_eigen_gen(int, MATRIX*, MATRIX*, MATRIX*, MATRIX*);
void solve_eigen_gen(int, MATRIX&, MATRIX&, MATRIX&, MATRIX&);


void solve_eigen(int Norb, MATRIX* H, MATRIX* S, CMATRIX* E, CMATRIX* C);
void solve_eigen(int Norb, MATRIX& H, MATRIX& S, CMATRIX& E, CMATRIX& C);
void solve_eigen_gen(int Norb, MATRIX* H, MATRIX* S, CMATRIX* E, CMATRIX* C);
void solve_eigen_gen(int Norb, MATRIX& H, MATRIX& S, CMATRIX& E, CMATRIX& C);


void solve_eigen(int Norb, CMATRIX* H, CMATRIX* S, CMATRIX* E, CMATRIX* C);
void solve_eigen(int Norb, CMATRIX& H, CMATRIX& S, CMATRIX& E, CMATRIX& C);
void solve_eigen_gen(int Norb, CMATRIX* H, CMATRIX* S, CMATRIX* E, CMATRIX* C);
void solve_eigen_gen(int Norb, CMATRIX& H, CMATRIX& S, CMATRIX& E, CMATRIX& C);


void solve_eigen(int, MATRIX*, MATRIX*, MATRIX*);


void sqrt_matrix(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half, double thresh);
void sqrt_matrix(CMATRIX& S, CMATRIX& S_half, CMATRIX& S_i_half);
void inv_matrix(CMATRIX& S, CMATRIX& S_inv, double thresh);
void inv_matrix(CMATRIX& S, CMATRIX& S_inv);


void FullPivLU_rank_invertible(MATRIX& A, int& rank, int& is_inver);
void FullPivLU_rank_invertible(CMATRIX& A, int& rank, int& is_inver);
boost::python::list FullPivLU_rank_invertible(MATRIX& A);
boost::python::list FullPivLU_rank_invertible(CMATRIX& A);

void FullPivLU_decomposition(MATRIX& A, MATRIX& P, MATRIX& L, MATRIX& U, MATRIX& Q);
void FullPivLU_decomposition(CMATRIX& A, CMATRIX& P, CMATRIX& L, CMATRIX& U, CMATRIX& Q);

void FullPivLU_inverse(MATRIX& A, MATRIX& invA);
void FullPivLU_inverse(CMATRIX& A, CMATRIX& invA);


}// namespace libmeigen
}// namespace libmmath


#endif
