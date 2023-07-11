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
  \file libio.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libio.h"


/// liblibra 
namespace liblibra{

/// libio namespace
namespace libio{


void export_io_objects(){
/** 
  \brief Exporter of libio classes and functions

  It is empty - so far no functions/classes are exported to Python
  mostly, they are for C++ utilization

*/

}// export_io_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygio){
#else
BOOST_PYTHON_MODULE(libio){
#endif

  export_io_objects();

}



}// namespace libio
}// liblibra



