/*********************************************************************************
* Copyright (C) 2018-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libmontecarlo.cpp
  \brief The file implements Python export function
    
*/


#if defined(USING_PCH)
#include "../pch.h"
#else
#include <stdlib.h>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libmontecarlo.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libmontecarlo namespace
namespace libmontecarlo{

using namespace boost::python;
namespace bp = boost::python;

void export_montecarlo_objects(){
/** 
  \brief Exporter of libmontecarlo classes and functions

*/

  vector<MATRIX> (*expt_metropolis_gau_v1)
                 (Random& rnd, bp::object target_distribution, MATRIX& dof, bp::object distribution_params, 
                  int sample_size, int start_sampling, double gau_var) = &metropolis_gau;

  def("metropolis_gau",expt_metropolis_gau_v1);


}// export_montecarlo_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmontecarlo){
#else
BOOST_PYTHON_MODULE(libmontecarlo){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_montecarlo_objects();

}


}// libmontecarlo
}// liblibra

