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

#define BOOST_PYTHON_MAX_ARITY 30
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libbasis.h"

using namespace boost::python;
using namespace libmmath;


namespace libqchem{
namespace libbasis{




void export_basis_objects(){


  // Basis.cpp
  void (*expt_basis_params_s_v1)
  (int, vector<double>&, vector<double>&) = &basis_params_s;

  void (*expt_basis_params_p_v1)
  (int, vector<double>&, vector<double>&) = &basis_params_p;

  void (*expt_basis_params_d_v1)
  (int, vector<double>&, vector<double>&) = &basis_params_d;


  void (*expt_add_basis_ao_v1)
  (std::string Atom_name, VECTOR& R, std::string Atom_shell, int Nzeta, int Nquant,
  double  IP, double exp1, double exp2, double coeff1, double coeff2, vector<AO>& basis_ao) = &add_basis_ao;

  void (*expt_add_basis_ao_v2)
  (std::string Atom_name, VECTOR& R, std::string Atom_shell, int Nzeta, int Nquant,
  double  IP, double exp1, double exp2, double coeff1, double coeff2, boost::python::list basis_ao) = &add_basis_ao;



  int (*expt_num_valence_elec_v1)(int) = &num_valence_elec;



  // Basis_ovlp.cpp
  void (*expt_update_overlap_matrix_v1)(int,int,int,const VECTOR&,const VECTOR&,const VECTOR&,
  vector<AO>&,MATRIX&) = &update_overlap_matrix;

  // Basis_map.cpp
  void (*expt_show_mapping_v1)(const vector<vector<int> >&) = &show_mapping;



  def("basis_params_s", expt_basis_params_s_v1);
  def("basis_params_p", expt_basis_params_p_v1);
  def("basis_params_d", expt_basis_params_d_v1);

  def("add_basis_ao", expt_add_basis_ao_v1);
  def("add_basis_ao", expt_add_basis_ao_v2);
  def("num_valence_elec", expt_num_valence_elec_v1);

  def("update_overlap_matrix", expt_update_overlap_matrix_v1);

  def("show_mapping", expt_show_mapping_v1);




}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygbasis){
#else
BOOST_PYTHON_MODULE(libbasis){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_basis_objects();

}


}// namespace libbasis
}// namespace libqchem

