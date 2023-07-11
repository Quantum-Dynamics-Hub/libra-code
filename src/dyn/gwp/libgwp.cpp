/*********************************************************************************
* Copyright (C) 2015-2020 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libgwp.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

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

  def("gwp_value",  expt_gwp_value_v1);


  double (*expt_gwp_product_decomposition_v1)
  (double q1, double p1, double gamma1, double alp1,
   double q2, double p2, double gamma2, double alp2,
   double& q, double& p, double& gamma, double& alp) = &gwp_product_decomposition;

  double (*expt_gwp_product_decomposition_v2)
  (MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
   MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2,
   MATRIX& q,  MATRIX& p,  MATRIX& gamma,  MATRIX& alp) = &gwp_product_decomposition;

  def("gwp_product_decomposition",  expt_gwp_product_decomposition_v1);
  def("gwp_product_decomposition",  expt_gwp_product_decomposition_v2);


  ///=============== Overlaps ===================

  complex<double> (*expt_gwp_overlap_v1)
  (double q1, double p1, double gamma1, double alp1,
   double q2, double p2, double gamma2, double alp2) = &gwp_overlap;

  complex<double> (*expt_gwp_overlap_v2)
  (MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
   MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2) = &gwp_overlap;

  complex<double> (*expt_gwp_overlap_v3)
  (MATRIX& R1, MATRIX& P1, double gamma1, 
   MATRIX& R2, MATRIX& P2, double gamma2, double alp, double hbar) = &gwp_overlap;

  complex<double> (*expt_gwp_overlap_v4)
  (MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
   MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2,
   MATRIX& q3, MATRIX& p3, MATRIX& gamma3, MATRIX& alp3) = &gwp_overlap;

  def("gwp_overlap",  expt_gwp_overlap_v1);
  def("gwp_overlap",  expt_gwp_overlap_v2);
  def("gwp_overlap",  expt_gwp_overlap_v3);
  def("gwp_overlap",  expt_gwp_overlap_v4);


  CMATRIX (*expt_gwp_overlap_matrix_v1)
  (MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
   MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2) = &gwp_overlap_matrix;

  CMATRIX (*expt_gwp_overlap_matrix_v2)
  (MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1, vector<int>& state1,
   MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2, vector<int>& state2) = &gwp_overlap_matrix;

  def("gwp_overlap_matrix",  expt_gwp_overlap_matrix_v1);
  def("gwp_overlap_matrix",  expt_gwp_overlap_matrix_v2);


  ///=============== Transition dipole moments ===================
  complex<double> (*expt_gwp_dipole_v1)
  (double q1, double p1, double gamma1, double alp1,
   double q2, double p2, double gamma2, double alp2) = &gwp_dipole;

  CMATRIX (*expt_gwp_dipole_v2)
  (MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
   MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2) = &gwp_dipole;

  CMATRIX (*expt_gwp_dipole_v3)
  (MATRIX& R1, MATRIX& P1, double gamma1, MATRIX& R2, MATRIX& P2, double gamma2, double alp, double hbar) = &gwp_dipole; 


  def("gwp_dipole", expt_gwp_dipole_v1);
  def("gwp_dipole", expt_gwp_dipole_v2);
  def("gwp_dipole", expt_gwp_dipole_v3);


  ///=============== Derivative couplings ===================

  complex<double> (*expt_gwp_coupling_v1)
  (double q1, double p1, double gamma1, double alp1,
   double q2, double p2, double gamma2, double alp2) = &gwp_coupling;

  CMATRIX (*expt_gwp_coupling_v2)
  (MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
   MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2) = &gwp_coupling;

  CMATRIX (*expt_gwp_coupling_v3)
  (MATRIX& R1, MATRIX& P1, double gamma1, MATRIX& R2, MATRIX& P2, double gamma2, double alp, double hbar) = &gwp_coupling; 


  def("gwp_coupling", expt_gwp_coupling_v1);
  def("gwp_coupling", expt_gwp_coupling_v2);
  def("gwp_coupling", expt_gwp_coupling_v3);

  ///=============== Kinetic energy ===================

  complex<double> (*expt_gwp_kinetic_v1)
  (double q1, double p1, double gamma1, double alp1,
   double q2, double p2, double gamma2, double alp2) = &gwp_kinetic;

  complex<double> (*expt_gwp_kinetic_v2)
  (MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
   MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2) = &gwp_kinetic;

  complex<double> (*expt_gwp_kinetic_v3)
  (MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
   MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2,
   MATRIX& iM) = &gwp_kinetic;

  complex<double> (*expt_gwp_kinetic_v4)
  (MATRIX& R1, MATRIX& P1, double gamma1, MATRIX& R2, MATRIX& P2, double gamma2, double alp, double hbar) = &gwp_kinetic;

  def("gwp_kinetic",  expt_gwp_kinetic_v1);
  def("gwp_kinetic",  expt_gwp_kinetic_v2);
  def("gwp_kinetic",  expt_gwp_kinetic_v3);
  def("gwp_kinetic",  expt_gwp_kinetic_v4);


  CMATRIX (*expt_gwp_kinetic_matrix_v1)
  (MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
   MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2, 
   MATRIX& invM) = &gwp_kinetic_matrix;

  def("gwp_kinetic_matrix",  expt_gwp_kinetic_matrix_v1);

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

