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
#include "libio.h"


namespace libio{


void export_io_objects(){

// Nothing to export to Python


}// export_io_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygio){
#else
BOOST_PYTHON_MODULE(libio){
#endif

  export_io_objects();

}



}// namespace libio




