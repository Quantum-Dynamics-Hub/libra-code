/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libhamiltonian.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libhamiltonian.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libhamiltonian namespace
namespace libhamiltonian{

using namespace libhamiltonian_generic;
//using namespace libhamiltonian_model;
//using namespace libhamiltonian_atomistic;
//using namespace libhamiltonian_extern;

void export_Hamiltonian_objects(){
/** 
  \brief Exporter of the libhamiltonian classes and functions

*/


//  export_hamiltonian_generic_objects();
//  export_hamiltonian_model_objects();
//  export_hamiltonian_atomistic_objects();
//  export_hamiltonian_extern_objects();

  export_nhamiltonian_generic_objects();

}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian){
#else
BOOST_PYTHON_MODULE(libhamiltonian){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Hamiltonian_objects();

}


}// namespace libhamiltonian
}// liblibra

