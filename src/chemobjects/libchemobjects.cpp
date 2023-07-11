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
  \file libchemobjects.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libchemobjects.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libchemobjects
namespace libchemobjects{

using namespace libuniverse;
using namespace libmol;
using namespace libchemsys;


void export_chemobjects_objects(){
/** 
  \brief Exporter of libchemobjects classes and functions

*/

  export_Universe_objects();
  export_Mol_objects();
  export_Chemsys_objects();

}// export_chemobjects_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygchemobjects){
#else
BOOST_PYTHON_MODULE(libchemobjects){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
//  export_Mathematics_objects();
  export_chemobjects_objects();

}


}// libchemobjects
}// liblibra

