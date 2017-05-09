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
#include "libcontext.h"

/// liblibra namespace
namespace liblibra{

namespace libcontext{


void export_context_objects(){

  export_Context_objects();
  export_ctx_Control_Parameters_objects();

}// export_context_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygcontext){
#else
BOOST_PYTHON_MODULE(libcontext){
#endif

  export_context_objects();

}



}// namespace libcontext
}// liblibra



