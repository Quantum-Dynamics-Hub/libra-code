/*********************************************************************************
* Copyright (C) 2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file grad_descent.cpp
  \brief The file implements the gradient descent minimization algorithm
    
*/

#include "opt.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;

/// libopt namespace 
namespace libopt{

namespace bp = boost::python;

MATRIX grad_descent(bp::object grad_function, MATRIX& dof, bp::object funct_params, double grad_tol, double step_size, int max_steps){
/**
  The Python function should correspond to the following C++ signature:

  MATRIX dof_grad = py_funct(MATRIX& dof, bp::object params);

  \param[in] grad_function - the Python function that computes the gradient of the
   function which we minimize
  \param[in] dof - Degrees of Freedom to be varied
  \param[in] funct_params - the parameters of the function
  \param[in] grad_tol - tolerance of the gradient (stopping critetion)
  \param[in] step_size - algorithm parameter
*/


  int ncols = dof.n_cols; 
  int nrows = dof.n_rows; 
  
  MATRIX res(nrows, ncols);  // the result
  MATRIX grd(nrows, ncols);  // the gradient

  res = MATRIX(dof);
  double err = 2.0*grad_tol;
  int iter = 0;

  while(err>grad_tol && iter<max_steps){

      // Call the Python function with such arguments
      grd = bp::extract<MATRIX>(grad_function(res, funct_params));  

      // Update the coordinates
      res = res - grd * step_size;

      // Compute the error
      err = grd.max_elt();

      iter++;
  }

  return res;

}



}// namespace libopt
}// liblibra


