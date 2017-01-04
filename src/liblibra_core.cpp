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
  \file liblibra_core.cpp
  \brief The file exports all core functionality of the Libra code
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "liblibra_core.h"

using namespace boost::python;

using namespace libcalculators;
using namespace libcell;
using namespace libchemobjects;
using namespace libcontext;
using namespace libconverters;
using namespace libdyn;
using namespace libhamiltonian;
using namespace libio;
using namespace libmmath;
using namespace libpot;
using namespace libqchem;
using namespace libqchem_tools;
using namespace libscripts;
using namespace libsolvers;
using namespace libutil;



void export_libra_core_objects(){
/** 
  \brief Exporter of Timer class and other mathematical libraries and their components

*/

  export_calculators_objects();
  export_context_objects();
  export_io_objects();
  export_Mathematics_objects();

  export_Cell_objects();
  export_chemobjects_objects();

  export_converters_objects();
  export_Dyn_objects();
  export_Hamiltonian_objects();


  export_Pot_objects();
  export_Qchem_objects();
  export_qchem_tools_objects();
  export_scripts_objects();
  export_solvers_objects();

  export_util_objects();

}// export_libra_core_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyglibra_core){
#else
BOOST_PYTHON_MODULE(liblibra_core){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
  export_libra_core_objects();

}

