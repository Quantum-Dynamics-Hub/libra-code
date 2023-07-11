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

#if defined(USING_PCH)
#include "../pch.h"
#else

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#endif 

#include "libgraph.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libgraph namespace
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
}// liblibra

