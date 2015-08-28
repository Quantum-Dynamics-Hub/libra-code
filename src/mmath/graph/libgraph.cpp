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
#include "libgraph.h"

using namespace boost::python;


namespace libmmath{
namespace libgraph{

void export_GRAPH_objects(){
}



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyggraph){
#else
BOOST_PYTHON_MODULE(libgraph){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
  export_GRAPH_objects();

}



}// libgraph
}// libmmath

