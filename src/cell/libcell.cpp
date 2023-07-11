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
  \file libcell.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

#include "libcell.h"

/// liblibra namespace
namespace liblibra{


using namespace boost::python;

/// libcell namespace
namespace libcell{

void export_Cell_objects(){
/** 
  \brief Exporter of libcell classes and functions

*/



//int (DATA::*ScaleData1)(double)                = &DATA::ScaleData;
//int (DATA::*ScaleData2)(double,double)         = &DATA::ScaleData;


  class_<Cell>("Cell",init<>())
      .def(init<VECTOR&, VECTOR&, VECTOR&, double>())
      .def("init", &Cell::init)
      .def("__copy__", &generic__copy__<Cell>) 
      .def("__deepcopy__", &generic__deepcopy__<Cell>)
//      .def_readwrite("x",&VECTOR::x)
            
  ;

  VECTOR (*expt_max_vector_v1)(VECTOR t1,VECTOR t2,VECTOR t3) = &max_vector;
  boost::python::list (*expt_apply_pbc_v1)(MATRIX3x3 H, boost::python::list in, boost::python::list t) = &apply_pbc;
  boost::python::list (*expt_serial_to_vector_v1)(int c,int Nx,int Ny,int Nz) = &serial_to_vector;
  boost::python::list (*expt_serial_to_vector_symm_v1)(int c,int Nx,int Ny,int Nz) = &serial_to_vector_symm;
  boost::python::list (*expt_form_neibc_v1)(int c,int Nx,int Ny,int Nz,double cellx,double celly,double cellz,double Roff) = &form_neibc;
  MATRIX (*expt_fold_coords_v1)(MATRIX& R, MATRIX3x3& box, std::string pbc_type) = &fold_coords;

  def("max_vector", expt_max_vector_v1);
  def("apply_pbc", expt_apply_pbc_v1);
  def("serial_to_vector",expt_serial_to_vector_v1);
  def("serial_to_vector_symm",expt_serial_to_vector_symm_v1);
  def("form_neibc",expt_form_neibc_v1);

  def("find_min_shell",find_min_shell);

  def("make_nlist",make_nlist);
  def("make_nlist_auto",make_nlist_auto);

  def("bruteforce",bruteforce);  
  def("energy",energy);

  def("fold_coords",expt_fold_coords_v1);  



} // export_Cell_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygcell){
#else
BOOST_PYTHON_MODULE(libcell){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Cell_objects();

}

}// namespace libcell
}// liblibra
