/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file integrators.cpp
  \brief The file implements some ODE integrators

*/
#include "integrators.h"


/// liblibra namespace
namespace liblibra{


/// libdyn namespace
namespace libintegrators{

namespace bp = boost::python;


CMATRIX RK4(CMATRIX& q, double dt, bp::object compute_derivatives, bp::object function_params){
/**
   This function solves a general multidimensional ODE:  dq/dt = f(t,q) for a 1 timestep dt

   It propagates the q variable (complex matrix of arb dimensions by the 4-th order Runge-Kutta algorithm)

   `compute_derivatives` represents the RHS function f(t,q) and is assumed to return CMATRIX object of the same size as the input q

  rho=rho+dt/6.d0*(k1+2*k2+2*k3+k4)
 
*/
    int nrows, ncols;
    nrows = q.n_rows;
    ncols = q.n_cols;

    CMATRIX der1(nrows, ncols);
    CMATRIX der2(nrows, ncols);
    CMATRIX der3(nrows, ncols);
    CMATRIX der4(nrows, ncols);
    CMATRIX tmp(nrows, ncols);
    CMATRIX res(nrows, ncols);

    // Call the Python function with such arguments
    tmp = q;
    der1 = extract<CMATRIX>(compute_derivatives(tmp, function_params));  

    tmp = q + 0.5*dt*der1;
    der2 = extract<CMATRIX>(compute_derivatives(tmp, function_params));  

    tmp = q + 0.5*dt*der2;
    der3 = extract<CMATRIX>(compute_derivatives(tmp, function_params));  

    tmp = q + dt*der3;
    der4 = extract<CMATRIX>(compute_derivatives(tmp, function_params));  

    res = q + (dt/6.0)*(der1 + 2.0*der2 + 2.0*der3 + der4);
 
    return res; 

}



}// namespace libdyn
}// libintegrators


