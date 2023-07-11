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
  \file libbasis.cpp
  \brief The file implements Python export function
    
*/

#define BOOST_PYTHON_MAX_ARITY 30

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

#include "libbasis.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;
using namespace liblinalg;


/// libbasis namespace
namespace libbasis{

using namespace libergoscf;
//using libergoscf::basisset_shell_struct;


void export_basis_objects(){
/** 
  \brief Exporter of libbasis classes and functions

*/



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

  void (*expt_MO_overlap_v1)(MATRIX& Smo, vector<AO>& ao_i, vector<AO>& ao_j, MATRIX& Ci, MATRIX& Cj,
  vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2) = &MO_overlap;

  void (*expt_MO_overlap_v2)(CMATRIX& Smo, vector<AO>& ao_i, vector<AO>& ao_j, CMATRIX& Ci, CMATRIX& Cj,
  vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2) = &MO_overlap;

  void (*expt_MO_overlap_v3)(MATRIX& Smo, MATRIX& Ci, MATRIX& Cj, 
  vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2) = &MO_overlap;

  void (*expt_MO_overlap_v4)(CMATRIX& Smo, CMATRIX& Ci, CMATRIX& Cj,
  vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2) = &MO_overlap;


  complex<double> (*expt_SD_overlap_v1)(SD& sd_i, SD& sd_j) = &SD_overlap;

  CMATRIX (*expt_SD_overlap_v2)(vector<SD>& sd_i, vector<SD>& sd_j) = &SD_overlap;

  void (*expt_SD_overlap_v3)(CMATRIX& SD_ovlp, vector<SD>& sd_i, vector<SD>& sd_j) = &SD_overlap;




  // Basis_map.cpp
  void (*expt_show_mapping_v1)(const vector<vector<int> >&) = &show_mapping;

  // Basis_nac.cpp
  void (*expt_update_derivative_coupling_matrix_v1)
  (int x_period,int y_period,int z_period,const VECTOR& t1, const VECTOR& t2, const VECTOR& t3,
   vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
   vector<AO>& basis_ao, int c, MATRIX& Dao_x, MATRIX& Dao_y, MATRIX& Dao_z
  ) = &update_derivative_coupling_matrix;





  def("basis_params_s", expt_basis_params_s_v1);
  def("basis_params_p", expt_basis_params_p_v1);
  def("basis_params_d", expt_basis_params_d_v1);

  def("add_basis_ao", expt_add_basis_ao_v1);
  def("add_basis_ao", expt_add_basis_ao_v2);
  def("num_valence_elec", expt_num_valence_elec_v1);

  def("update_overlap_matrix", expt_update_overlap_matrix_v1);
  def("MO_overlap", expt_MO_overlap_v1);
  def("MO_overlap", expt_MO_overlap_v2);
  def("MO_overlap", expt_MO_overlap_v3);
  def("MO_overlap", expt_MO_overlap_v4);

  def("SD_overlap", expt_SD_overlap_v1);
  def("SD_overlap", expt_SD_overlap_v2);
  def("SD_overlap", expt_SD_overlap_v3);


  def("show_mapping", expt_show_mapping_v1);

  def("update_derivative_coupling_matrix", expt_update_derivative_coupling_matrix_v1);



  class_<basisset_shell_struct>("basisset_shell_struct",init<int>())

      .def_readwrite("shell_type", &basisset_shell_struct::type)
      .def_readwrite("contrCount", &basisset_shell_struct::contrCount)
      .def_readwrite("shell_ID", &basisset_shell_struct::shell_ID)
      .def_readwrite("exponentList", &basisset_shell_struct::exponentList)
      .def_readwrite("coeffList", &basisset_shell_struct::coeffList)

  ;

  class_< basisset_shell_structList >("basisset_shell_structList")
      .def(vector_indexing_suite< basisset_shell_structList >())
  ;

  class_< basisset_shell_structMap >("basisset_shell_structMap")
      .def(vector_indexing_suite< basisset_shell_structMap >())
  ;



  class_<basisset_atom_struct>("basisset_atom_struct",init<int, int>())

      .def_readwrite("noOfShells", &basisset_atom_struct::noOfShells)
      .def_readwrite("shells", &basisset_atom_struct::shells)

  ;

  class_< basisset_atom_structList >("basisset_atom_structList")
      .def(vector_indexing_suite< basisset_atom_structList >())
  ;

  class_< basisset_atom_structMap >("basisset_atom_structMap")
      .def(vector_indexing_suite< basisset_atom_structMap >())
  ;



  class_<basisset_struct>("basisset_struct",init<int, int, int>())
      .def_readwrite("atoms", &basisset_struct::atoms)
      .def("read_basisset_file", &basisset_struct::read_basisset_file);

  ;

/*
  class_< basisset_structList >("basisset_structList")
      .def(vector_indexing_suite< basisset_structList >())
  ;

  class_< basisset_structMap >("basisset_structMap")
      .def(vector_indexing_suite< basisset_structMap >())
  ;
*/


//  int (*expt_read_basisset_file_v1)(basisset_struct* result, const char* fileName, int dirc, const char *dirv[], int print_raw) = &read_basisset_file;
//  def("read_basisset_file", expt_read_basisset_file_v1);

/*
  int (*expt_read_basisset_file_v1)(basisset_struct* result, string fileName, int print_raw) = &read_basisset_file;

  basisset_struct (*expt_read_basisset_file_v2)(std::string fileName,int print_raw) = &read_basisset_file;
  def("read_basisset_file", expt_read_basisset_file_v2);
*/

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
}// namespace liblibra

