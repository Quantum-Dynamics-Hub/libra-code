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
  \file libqchem.cpp
  \brief The file implements Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libqchem.h"
using namespace boost::python;

/// libqchem namespace
namespace libqchem{

using namespace libmolint;
using namespace libqobjects;
using namespace libbasis;


void export_Qchem_objects(){
/** 
  \brief Exporter of libqchem classes and functions

*/

  export_molint_objects();
  export_qobjects_objects();
  export_basis_objects();

}// export_Qchem_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygqchem){
#else
BOOST_PYTHON_MODULE(libqchem){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Qchem_objects();

}


}// libqchem

