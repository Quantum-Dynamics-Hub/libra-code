/*********************************************************************************
* Copyright (C) 2015-2021 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
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

/// liblibra namespace
namespace liblibra{


using namespace boost::python;

using namespace libutil;
using namespace libio;
using namespace libcontext;
using namespace libcommon_types;

using namespace libann;
using namespace libdata;
using namespace libgraph;
using namespace liblinalg;
using namespace libmeigen;
using namespace liboperators;
using namespace librandom;
using namespace libspecialfunctions;
using namespace libsymmetry;

using namespace libmolint;
using namespace libqobjects;
using namespace libbasis;
using namespace libbasis_setups;
using namespace liblibint2_wrappers;

using namespace libcalculators;

using namespace librigidbody;

using namespace libchemobjects;
using namespace libcell;
using namespace libpot;
using namespace libforcefield;

using namespace libcontrol_parameters;
using namespace libnhamiltonian;

using namespace libatomistic;
using namespace libmodels;

using namespace libfgr;
using namespace libivr;
using namespace libdyn;

using namespace libconverters;
using namespace libscripts;
using namespace libqchem_tools;
using namespace libsolvers;
using namespace libintegrators;

using namespace libmontecarlo;
using namespace libopt;


void export_libra_core_objects(){
/** 
  \brief Exporter of Timer class and other mathematical libraries and their components

*/

  export_util_objects();
  export_io_objects();
  export_context_objects();
  export_common_types_objects();

  export_NeuralNetwork_objects();
  export_Data_objects();
  export_GRAPH_objects();
  export_linalg_objects();
  export_mEigen_objects();
  export_Operators_objects();
  export_Random_objects();
  export_SpecialFunctions_objects();
  export_symmetry_objects();
  export_timer_objects();

  export_molint_objects();
  export_qobjects_objects();
  export_basis_objects();
  export_libint2_wrappers_objects();


  export_calculators_objects();

  export_models_objects();

  export_RigidBody_objects();
  export_chemobjects_objects();
  export_Cell_objects();

  export_Pot_objects();
  export_forcefield_objects();


  export_Control_Parameters_objects();
  export_Model_Parameters_objects();
  export_basis_setups_objects();

  export_atomistic_objects();

  export_nhamiltonian_objects();

  export_fgr_objects();
  export_ivr_objects();

  export_Dyn_objects();

  export_converters_objects();
  export_scripts_objects();
  export_qchem_tools_objects();

  export_solvers_objects();

  export_integrators_objects();


  export_montecarlo_objects();
  export_opt_objects();


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

}// liblibra
