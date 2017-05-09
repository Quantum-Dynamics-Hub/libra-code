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
/**
  \file libgwp.cpp
  \brief The file implements Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libgwp.h"


/// liblibra namespace
namespace liblibra{


using namespace boost::python;

/// libdyn namespace
namespace libdyn{

/// libgwp namespace
namespace libgwp{


void export_gwp_objects(){
/** 
  \brief Exporter of libgwp classes and functions

*/

  complex<double> (*expt_gwp_value_v1)(MATRIX& r, MATRIX& R, MATRIX& P, double gamma,  double alp, double hbar) = &gwp_value;

  complex<double> (*expt_gwp_overlap_v1)
  (MATRIX& R1, MATRIX& P1, double gamma1, MATRIX& R2, MATRIX& P2, double gamma2, double alp, double hbar) = &gwp_overlap;

  CMATRIX (*expt_gwp_dipole_v1)
  (MATRIX& R1, MATRIX& P1, double gamma1, MATRIX& R2, MATRIX& P2, double gamma2, double alp, double hbar) = &gwp_dipole; 

  CMATRIX (*expt_gwp_coupling_v1)
  (MATRIX& R1, MATRIX& P1, double gamma1, MATRIX& R2, MATRIX& P2, double gamma2, double alp, double hbar) = &gwp_coupling; 

  complex<double> (*expt_gwp_kinetic_v1)
  (MATRIX& R1, MATRIX& P1, double gamma1, MATRIX& R2, MATRIX& P2, double gamma2, double alp, double hbar) = &gwp_kinetic;


  def("gwp_value",  expt_gwp_value_v1);
  def("gwp_overlap",  expt_gwp_overlap_v1);
  def("gwp_dipole", expt_gwp_dipole_v1);
  def("gwp_coupling", expt_gwp_coupling_v1);
  def("gwp_kinetic",  expt_gwp_kinetic_v1);



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyggwp){
#else
BOOST_PYTHON_MODULE(libgwp){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_gwp_objects();

}


}// namespace libgwp
}// namespace libdyn
}// liblibra

