/*********************************************************************************
* Copyright (C) 2022 Matthew Dutra, Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libqtag.cpp
  \brief The file implements Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libqtag.h"


/// liblibra namespace
namespace liblibra{


using namespace boost::python;

/// libdyn namespace
namespace libdyn{

/// libqtag namespace
namespace libqtag{


void export_qtag_objects(){
/** 
  \brief Exporter of libqtag classes and functions

*/

  double (*expt_qtag_momentum_v1)(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff) = &qtag_momentum;

  def("qtag_momentum",  expt_qtag_momentum_v1);



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygqtag){
#else
BOOST_PYTHON_MODULE(libqtag){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_qtag_objects();

}


}// namespace libqtag
}// namespace libdyn
}// liblibra

