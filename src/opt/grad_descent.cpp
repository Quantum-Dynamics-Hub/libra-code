/*********************************************************************************
* Copyright (C) 2017-2024 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
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


void grad_step(MATRIX& A, MATRIX& x, MATRIX& x_new, double dt, double& L, MATRIX& mu, int opt,
               double& Lag, double& L0, double& L1, double& L2, double& L3, double dT, int approach_option){
    /**
    Optimization of |Ax - y|^2 subject to constraints
                                                                                       Evolve:
    opt = 0 - unconstrained optimization  L0 term                                      A
    opt = 1 -


    opt = 1 - add only lambda-constraint  L0, L1 terms                                 A, L
    opt = 2 - add only mu-constraint      L0, L2 terms                                 A, mu
    opt = 3 - add both constraints        L0, L1, L2 terms                             A, L, mu
    opt = 4 - the constraint on the corrected fluxes            L0, L3                 A, L
    opt = 5 - same as 4 + sum of row should be close to 1       L0, L3, L4             A, L, mu
    opt = 6 - same as 4, but scale fluxes by J_out (see below)  L0, L3 (scaled)        A, L


    */
    int i, j, k;
    int n = A.n_cols;

    //*****************************************
    int opt2 = 0; // 0 - linear, 1 - quadratic
    //*****************************************

    MATRIX y(n, 1);
    MATRIX num_coeff(n,n);
    MATRIX denom_coeff(n,n);

    if(approach_option==0){ y = x_new; }

    vector<double> J_out(n, 0.0); // total outfluxes from all states
    if(approach_option==1 or approach_option==2){
      y = (x_new-x)/dT;

      for(j=0;j<n;j++){
        double prob = 0.0;
        if(x.get(j,0)>0.0){ prob = -dT*y.get(j,0)/x.get(j,0); }
        else{ prob = 0.0; }

        if(prob<0.0){ prob = 0.0; }
        if(prob>1.0){ prob = 1.0; }
        J_out[j] = prob;

        for(i=0;i<n;i++){
          if(i==j){  num_coeff.set(i, j, 0.0);  denom_coeff.set(i,j, 0.0); }
          else{  num_coeff.set(i, j, prob);  denom_coeff.set(i,j, 1.0);}
        }
      }
    }

//    exit(0);

    MATRIX id(n,n); id.identity();

    // Total ones
    Lag = 0.0;
    MATRIX gradA(n,n);   gradA = 0.0;
    MATRIX gradmu(n,1);  gradmu = 0.0;
    double gradL;        gradL = 0.0;

    // Temporary ones
    MATRIX _gradA(n,n);
    MATRIX _gradmu(n,1);
    double _gradL;

    L0 = 0.0; L1 = 0.0; L2 = 0.0; L3 = 0.0;
                                                   
    // For all methods, we have the |Ax - y|^2 -> min condition
    if(opt>=0){
      _gradA = 0.0;
      L0 = derivs0(x, y, A, _gradA);
      gradA += _gradA;
      Lag += L0;
    }

    if(opt==1){
      L1 = derivs1_new(x, y, A, _gradA, L, _gradL, num_coeff, denom_coeff, 0);
      gradA += _gradA; gradL+= _gradL; Lag += L1;
    }
    if(opt==101){
      L1 = derivs1_new(x, y, A, _gradA, L, _gradL, num_coeff, denom_coeff, 1);
      gradA += _gradA; gradL+= _gradL; Lag += L1;
    }
/*
    if(opt==1 or opt==3){
      _gradA = 0.0; _gradL = 0.0;
      L1 = derivs1(x, y, A, _gradA, L, _gradL, opt2, 0); // Last argument: 0 - no diagonal elements, 1 - with diagonal elts.
      gradA += _gradA; gradL += _gradL;
      Lag += L1;
    }
*/
    if(opt==2 or opt==3){
      _gradA = 0.0; _gradmu = 0.0;
      L2 = derivs2(x, y, A, _gradA, mu, _gradmu, opt2, 0); // Last argument: 0 - constraint on columns, 1 - constraint on rows
      gradA += _gradA; gradmu += _gradmu;
      Lag += L2;
    }

    if(opt==4 or opt==5 or opt==6){

      MATRIX s(n,n);
      if(opt==4){   for(i=0;i<n;i++){  for(j=0;j<n;j++){  s.set(i,j, 1.0); }  }      }
      else if(opt==6){

        for(i=0;i<n;i++){  for(j=0;j<n;j++){  s.set(i,j, 1.0); }  }

        if(approach_option==1 or approach_option==2){
          for(i=0;i<n;i++){
            double prob = 0.0; // total outflux from the state i
            if(x.get(i,0)>0.0){ prob = -dT*y.get(i,0)/x.get(i,0); }
            if(prob<0.0){ prob = 0.0; }
            if(prob>1.0){ prob = 1.0; }

            for(j=0;j<n;j++){
              if(j==i){ s.set(j, j, 0.0); }
              else{  s.set(j, i, prob); }
            }// for j
          }// for i
        }// if
      }

      _gradA = 0.0; _gradL = 0.0;
      L3 = derivs3(x, y, A, _gradA, L, _gradL, opt2, 1, s); // second to the last argument: 0 - no diagonal elements, 1 - with diagonal elts.
      gradA += _gradA; gradL += _gradL;
      Lag += L3;
    }

    if(opt==5){
      _gradA = 0.0; _gradmu = 0.0;
      L2 = derivs4(x, y, A, _gradA, mu, _gradmu, opt2); //
      gradA += _gradA; gradmu += _gradmu;
      Lag += L2;

    }// opt == 5


  A = A - dt * gradA;
  if(opt==0 or opt==1 or opt==101 or opt==3 or opt==4 or opt==5 or opt==6){
    L = L - dt * gradL;
  }else{ L = 0.0; }

  if(opt==2 or opt==3 or opt==5){
    mu = mu - dt* gradmu;
  }else{ mu = 0.0; }

}

MATRIX run_opt(MATRIX& x, MATRIX& y, double dt, double nsteps, double err_tol, int opt, vector<double>& err, double dT, int approach_option){
/**
  Looks for a solution `A` of for the problem ||Ax - y||^2 -> min, where x and y are specified inputs

  \param[in] x MATRIX(N, 1) - the x(t) part of the Ax = y equation
  \param[in] y MATRIX(N, 1) - the x(t+dT) part of the Ax = y equation
  \param[in] dt - the "time-step" for the steepest descent steps
  \param[in] nsteps - the maximal number of steps of the steepest descent
  \param[in] err_tol - the maximal allowed value of the error ||Ax - y||^2
  \param[in] opt - the parameter determining what kinds of constraints to also include during the optimization, see the 
    `grad_step` function for more details. 
  \param[out] err 5-element array that will keep the errors due to different terms of the overall Lagrangian
     the first element is always the value of the total Lagrangian; the second element is the ||Ax - y ||^2 value
  \param[in] dT - is the external time-step separating x and y
  \param[in] approach_option - defines what kind of data constitutes the vectors x and y, it also defines the meaning of
     the resulting matrix A (returned by this function)
     approach_option = 0, y = rho(t+dT) and x = rho(t), so A is considered a "transition probability" matrix 
                          and is initialized as the identity matrix
     approach_option = 1 or 2, y = (rho(t+dT) - rho(t))/dT and x = rho(t), so A is considered the matrix of rate constants and 
                          is initialized to be zero matrix.

  Returns: the matrix A of type MATRIX(N,N) that minimizes the identified Lagrangian

*/

    int i,j;
    int n = x.n_rows;
    MATRIX A(n,n);  A = 0.0;
    if(approach_option==0){  A.identity(); }
    else if(approach_option==1){
      for(i=0;i<n;i++){
        for(j=0;j<n;j++){  A.set(i,j, 0.0 ); }
      }
    }
    else if(approach_option==2){
      for(i=0;i<n;i++){
        for(j=0;j<n;j++){  
          if(j-i==1){  A.set(i,j, 1.0); }
          else if(j-i==-1){ A.set(i,j, -1.0); }
          else{         A.set(i,j, 0.0 );  }
        }
      }
    }

    double Lag, L0, L1, L2, L3;
    double L = 1.0;
    MATRIX mu(n,1);
    for(i=0; i<n; i++){ mu.set(i, 0, 1.0); }

    Lag = 2.0*err_tol;
    while(i<nsteps and fabs(Lag) > err_tol){
      grad_step(A, x, y, dt, L, mu, opt, Lag, L0, L1, L2, L3, dT, approach_option);
      i++;
    }

    if( fabs(Lag) > 1.0){ cout<<"WARNING: run_opt did not converge. Error: "<<fabs(Lag)<<"\n";  }

    err[0] = Lag;
    err[1] = L0;
    err[2] = L1;
    err[3] = L2;
    err[4] = L3;

   return A;

}// run_opt





}// namespace libopt
}// liblibra


