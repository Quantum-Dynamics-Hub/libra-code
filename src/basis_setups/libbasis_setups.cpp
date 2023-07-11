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

#define BOOST_PYTHON_MAX_ARITY 30

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libbasis_setups.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;
using namespace liblibra;

namespace libbasis_setups{


void export_basis_setups_objects(){


  // Basis_STO_3G_DZ.cpp

  void (*expt_set_basis_STO_3G_DZ_v1)(vector<std::string>& at_type, vector<VECTOR>& R,  Model_Parameters& modpar, int verb,
  vector<AO>& basis_ao, int& Nelec, int& Norb, 
  vector<vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map) = &set_basis_STO_3G_DZ;

  boost::python::list (*expt_set_basis_STO_3G_DZ_v2)
  (vector<std::string>& at_type, vector<VECTOR>& R,  Model_Parameters& modpar, int verb) = &set_basis_STO_3G_DZ;

//  def("set_basis_STO_3G_DZ", expt_set_basis_STO_3G_DZ_v1);
  def("set_basis_STO_3G_DZ", expt_set_basis_STO_3G_DZ_v2);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygbasis_setups){
#else
BOOST_PYTHON_MODULE(libbasis_setups){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_basis_setups_objects();

}


}// namespace libbasis_setups
}// namespace liblibra

