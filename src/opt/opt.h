/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file opt.h
  \brief The file describes various optimization procedures
    
*/

#ifndef OPT_H
#define OPT_H


// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"



/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;

/// libopt namespace 
namespace libopt{

namespace bp = boost::python;

// In matrix_gradients.cpp
void normalize_transition_matrix(MATRIX& flux, MATRIX& num_coeff, MATRIX& denom_coeff);
void hyperbolic(double x, double alpha, double& val, double& deriv);
MATRIX transform_to_hyperbolic(MATRIX& A, double alp);
double derivs0(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA);
double derivs1(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA, double& L, double& gradL, int opt, int opt2);
double derivs1_new(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA, double& L, double& gradL, MATRIX& num_coeff, MATRIX& denom_coeff, int opt);
double derivs2(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA, MATRIX& mu, MATRIX& gradmu, int opt, int opt2);
double derivs3(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA, double& L, double& gradL, int opt, int opt2, MATRIX& s);
double derivs4(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA, MATRIX& mu, MATRIX& gradmu, int opt);


// In grad_descent.cpp
MATRIX grad_descent(bp::object grad_function, MATRIX& dof, bp::object funct_params, double grad_tol, double step_size, int max_steps);

void grad_step(MATRIX& A, MATRIX& x, MATRIX& x_new, double dt, double& L, MATRIX& mu, int opt,
               double& Lag, double& L0, double& L1, double& L2, double& L3, double dT, int approach_option);
MATRIX run_opt(MATRIX& x, MATRIX& x_new, double dt, double nsteps, double err_tol, int opt, vector<double>& err, double dT, int approach_option);


}// namespace libopt
}// liblibra

#endif // OPT_H
