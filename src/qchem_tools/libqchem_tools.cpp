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
  \file libqchem_tools.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libqchem_tools.h"

/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libqchem_tools namespace
namespace libqchem_tools{

void export_qchem_tools_objects(){
/** 
  \brief Exporter of libqchem_tools classes and functions

*/

  void (*expt_charge_density_v1)
  ( Electronic_Structure& el, System& syst, vector<AO>& basis_ao, Control_Parameters& prms) = &charge_density;

  void (*expt_charge_density_v2)
  (MATRIX& C, vector<listHamiltonian_QM>& ham, System& syst, vector<vector<int> >& active_orb, Control_Parameters& prms) = &charge_density;

  void (*expt_charge_density_v3)
  (MATRIX& C, boost::python::list ham, System& syst, boost::python::list active_orb, Control_Parameters& prms) = &charge_density;


  def("charge_density", expt_charge_density_v1);
  def("charge_density", expt_charge_density_v2);
  def("charge_density", expt_charge_density_v3);



  void (*expt_compute_dos_v1)
  ( Electronic_Structure& el, vector<AO>& basis_ao, Control_Parameters& prms,
    vector<int>& fragment, vector< vector<int> >& atom_to_ao_map
  ) = &compute_dos;

  void (*expt_compute_dos_v2)
  ( Electronic_Structure& el, vector<AO>& basis_ao, Control_Parameters& prms,
    boost::python::list fragment, vector< vector<int> >& atom_to_ao_map
  ) = &compute_dos;

  void (*expt_compute_dos_v3)
  ( Electronic_Structure& el, System& syst, vector<AO>& basis_ao, 
    Control_Parameters& prms, vector< vector<int> >& atom_to_ao_map
  ) = &compute_dos;


  def("compute_dos", expt_compute_dos_v1);
  def("compute_dos", expt_compute_dos_v2);
  def("compute_dos", expt_compute_dos_v3);



}// export_qchem_tools_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygqchem_tools){
#else
BOOST_PYTHON_MODULE(libqchem_tools){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_qchem_tools_objects();

}


}// libqchem_tools
}// liblibra

