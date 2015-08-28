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

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libconverters.h"
using namespace boost::python;


namespace libconverters{

//using namespace libuniverse;
//using namespace libmol;
//using namespace libchemsys;


void export_converters_objects(){

  void (*expt_system_to_nuclear_v1)(System& syst, Nuclear& nucl) = &system_to_nuclear;
  void (*expt_nuclear_to_system_v1)(Nuclear& nucl, System& syst) = &nuclear_to_system;


  def("system_to_nuclear", expt_system_to_nuclear_v1);
  def("nuclear_to_system", expt_nuclear_to_system_v1);

}// export_converters_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygconverters){
#else
BOOST_PYTHON_MODULE(libconverters){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
//  export_Mathematics_objects();
  export_converters_objects();

}


}// libconverters


