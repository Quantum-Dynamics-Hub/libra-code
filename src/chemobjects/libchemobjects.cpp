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
#include "libchemobjects.h"
using namespace boost::python;


namespace libchemobjects{

using namespace libuniverse;
using namespace libmol;
using namespace libchemsys;


void export_chemobjects_objects(){

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


