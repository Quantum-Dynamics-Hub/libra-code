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

/**
  \file libcalculators.cpp
  \brief The file that implements the exprots of libcalculator objects to Python
        
*/


#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libcalculators.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libcalculators namespace
namespace libcalculators{


void export_calculators_objects(){
/** 
  \brief Exporter of libcalculators classes and functions

*/


  //----------------- Fermi.cpp ------------------------------

  double (*expt_fermi_population_v1)
  (double e,double ef,double degen, double kT) = &fermi_population;

  double (*expt_fermi_integral_v1)
  (vector<double>& bnds, double ef, double degen, double kT) = &fermi_integral;

  double (*expt_fermi_integral_v2)
  (boost::python::list bnds, double ef, double degen, double kT) = &fermi_integral;

  double (*expt_fermi_energy_v1)
  (vector<double>& bnds,double Nel,double degen, double kT, double etol) = &fermi_energy;

  double (*expt_fermi_energy_v2)
  (boost::python::list bnds,double Nel,double degen, double kT, double etol) = &fermi_energy;

  //----------------- Bands.cpp --------------------------------

  boost::python::list (*expt_order_bands_v1)(MATRIX E) = &order_bands;

  boost::python::list (*expt_populate_bands_v1)
  (double Nel, double degen, double kT, double etol, int pop_opt, boost::python::list bands) = &populate_bands;


  //----------------- Density_Matrix.cpp ------------------------
  MATRIX (*expt_compute_density_matrix_v1)(boost::python::list occ, MATRIX C) = &compute_density_matrix;
  boost::python::list (*expt_Fock_to_P_v1)
  (MATRIX Fao, MATRIX Sao, double Nel, double degen, double kT, double etol, int pop_opt) = &Fock_to_P;


  //----------------- Excitations.cpp ---------------------------
  boost::python::list (*expt_excite_v1)(int I, int J, boost::python::list occ_ini) = &excite;


  //----------------- Energy_Electronic.cpp ---------------------
  double (*expt_energy_elec_v1)(MATRIX Pao,MATRIX Hao,MATRIX Fao) = &energy_elec;
  double (*expt_energy_elec_v2)(MATRIX P_alp, MATRIX P_bet, 
                                MATRIX Hao_alp, MATRIX Hao_bet,
                                MATRIX Fao_alp, MATRIX Fao_bet,
                                MATRIX dFao_alp_dP_alp, MATRIX dFao_alp_dP_bet,
                                MATRIX dFao_bet_dP_alp, MATRIX dFao_bet_dP_bet,
                                MATRIX temp
                               ) = &energy_elec;

  //----------------- Energy_Nuclear.cpp ---------------------
  double (*expt_energy_nucl_v1)(vector<VECTOR>& R, vector<double>& Zeff) = &energy_nucl;
  double (*expt_energy_nucl_v2)(vector<VECTOR>& R, vector<double>& Zeff, vector<VECTOR>& G) = &energy_nucl;



 
  def("fermi_population", expt_fermi_population_v1);
  def("fermi_integral", expt_fermi_integral_v1);
  def("fermi_integral", expt_fermi_integral_v2);
  def("fermi_energy", expt_fermi_energy_v1);
  def("fermi_energy", expt_fermi_energy_v2);

  def("order_bands", expt_order_bands_v1);
  def("populate_bands", expt_populate_bands_v1);

  def("compute_density_matrix",expt_compute_density_matrix_v1);
  def("Fock_to_P",expt_Fock_to_P_v1);

  def("excite",expt_excite_v1);

  def("energy_elec",expt_energy_elec_v1);
  def("energy_elec",expt_energy_elec_v2);

  def("energy_nucl",expt_energy_nucl_v1);
  def("energy_nucl",expt_energy_nucl_v2);



}// export_calculators_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygcalculators){
#else
BOOST_PYTHON_MODULE(libcalculators){
#endif

  export_calculators_objects();

}


}// namespace libcalculators

}// liblibra


