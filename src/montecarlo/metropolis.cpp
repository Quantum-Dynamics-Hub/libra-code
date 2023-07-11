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
  \file metropolis.cpp
  \brief The file implements the Metropolis MC sampling algorithm
    
*/

#include "montecarlo.h"




/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;
using namespace librandom;

/// libmontecarlo namespace 
namespace libmontecarlo{

namespace bp = boost::python;

vector<MATRIX> metropolis_gau(Random& rnd, bp::object target_distribution, MATRIX& dof, bp::object distribution_params, 
                              int sample_size, int start_sampling, double gau_var){
/**
  The Python function should correspond to the following C++ signature:

  double target_distribution(MATRIX& dof, bp::object params);

  \param[in] target_distribution - the Python function that computes the probability distribution function
  that depends on multiple DOFs
  \param[in] dof - Degrees of Freedom (DOF) on which the probability distribution function depends
  \param[in] distribution_params - the parameters of the distribution function
  \param[in] sample_size - how many multidimensional points to sample from the distribution
  \param[in] start_sampling - how many first moves to disregard
  \param[in] ksi - Gaussian variance - which is used to sample hop events
*/


  int ncols = dof.n_cols; 
  int nrows = dof.n_rows; 
  int ndof = ncols * nrows;

  vector<MATRIX> res;
  MATRIX s_old(nrows, ncols);  // the resuls
  MATRIX s_new(nrows, ncols);  // the resuls

  int act_sample = 0;  // actual number of sampled points
  int acc_count = 0;   // number of accepted counts

  // Initialization
  s_old = dof;
  double p_old = bp::extract<double>( target_distribution(s_old, distribution_params) );
  

  while(act_sample<sample_size){

      // Attempted move
      for(int i=0;i<ndof;i++){
          double si = s_old.get(i) + gau_var * rnd.normal();
          s_new.set(i, si);
      }
      
      // New probability
      double p_new = bp::extract<double>( target_distribution(s_new, distribution_params) );

      // Compute the acceptance probability
      double acc_prob = p_new /  p_old;
      if(acc_prob>1.0){  acc_prob = 1.0; }

      // Attempt the move
      double ksi = rnd.uniform(0.0, 1.0);
      if(ksi<acc_prob){

          acc_count++;

          // Successful move:
          if(acc_count>start_sampling){

              act_sample++;
              res.push_back(s_new);
              p_old = p_new;
              s_old = s_new;

          }// if 

      }// if 


  }// while 

  return res;

}





}// namespace libmontecarlo
}// liblibra


