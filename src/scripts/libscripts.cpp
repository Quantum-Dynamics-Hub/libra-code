/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libscripts.h"

/// liblibra namespace
namespace liblibra{

using namespace boost::python;


namespace libscripts{

using namespace libstate;


void export_scripts_objects(){

  export_state_objects();

}// export_scripts_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygscripts){
#else
BOOST_PYTHON_MODULE(libscripts){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
//  export_Mathematics_objects();
  export_scripts_objects();

}


}// libscripts
}// liblibra


