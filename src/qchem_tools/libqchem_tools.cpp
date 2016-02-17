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
/**
  \file libqchem_tools.cpp
  \brief The file implements Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libqchem_tools.h"
using namespace boost::python;

/// libqchem_tools namespace
namespace libqchem_tools{

void export_qchem_tools_objects(){
/** 
  \brief Exporter of libqchem_tools classes and functions

*/

  void (*expt_charge_density_v1)
  ( Electronic_Structure& el, System& syst, vector<AO>& basis_ao, Control_Parameters& prms) = & charge_density;

  def("charge_density", expt_charge_density_v1);



}// export_qchem_tools_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygqchem_tools){
#else
BOOST_PYTHON_MODULE(libqchem_tools){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_qchem_tools_objects();

}


}// libqchem_tools

