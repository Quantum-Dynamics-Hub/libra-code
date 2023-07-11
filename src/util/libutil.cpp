/*********************************************************************************
* Copyright (C) 2016-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libutil.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libutil.h"

/// liblibra namespace
namespace liblibra{


/// libutil namespace
namespace libutil{


void export_util_objects(){
/** 
  \brief Exporter of libutil classes and functions
*/

  void (*expt_show_vector_v1)(vector<int>& A) = &show_vector;  
  def("show_vector", expt_show_vector_v1);
  
  int (*expt_is_in_vector_v1)(int a, vector<int>& A) = &is_in_vector;
  int (*expt_is_in_vector_v2)(int a, vector<int>& A, vector<int>& indx) = &is_in_vector;
  boost::python::list (*expt_is_in_vector_v3)(int a, vector<int>& A) = &is_in_vector2;
  def("is_in_vector", expt_is_in_vector_v1);
  def("is_in_vector", expt_is_in_vector_v2);
  def("is_in_vector", expt_is_in_vector_v3);


  boost::python::list (*expt_is_repeating_v1)(vector<int>& A) = &is_repeating;
  def("is_repeating", expt_is_repeating_v1);  
  
  int (*expt_delta_v1)(vector<int>& A,vector<int>& B,int& a,int& b) = &delta;
  boost::python::list (*expt_delta_v2)(vector<int>& A,vector<int>& B) = &delta;
  def("delta", expt_delta_v1);  
  def("delta", expt_delta_v2);  

  void (*expt_split_line_v1)(std::string line, vector<std::string>& arr) = &split_line;
  void (*expt_split_line_v2)(std::string line,vector<std::string>& arr,char delim) = &split_line;
  def("split_line", expt_split_line_v1);  
  def("split_line", expt_split_line_v2);  


  int (*expt_is_equal_v1)
  (vector<int>& vec1, vector<int>& vec2) = &is_equal;
  def("is_equal", expt_is_equal_v1);
 
  int (*expt_is_included_v1)
  (vector<int>& vec1, vector<vector<int> >& vec) = &is_included;
  int (*expt_is_included_v2)
  (vector<int>& vec1, vector<vector<int> >& vec, int start, int num_of_elts) = &is_included; 
  def("is_included", expt_is_included_v1);
  def("is_included", expt_is_included_v2);

  int (*expt_is_present_v1)
  (vector< vector<int> >& vec, int i, int start, int end) = &is_present;
  def("is_present", expt_is_present_v1);

  int (*expt_sum_row_v1)
  (int row, vector<vector<int> >& vec) = &sum_row;
  def("sum_row", expt_sum_row_v1);




  vector<int> (*expt_allocate_1D_v1)
  (int sz1) = &allocate_1D;
  def("allocate_1D", expt_allocate_1D_v1);

  vector< vector<int> > (*expt_allocate_2D_v1)
  (int sz1, int sz2) = &allocate_2D;
  def("allocate_2D", expt_allocate_2D_v1);

  vector< vector< vector<int> > > (*expt_allocate_3D_v1)
  (int sz1, int sz2, int sz3) = &allocate_3D;
  def("allocate_3D", expt_allocate_3D_v1);




  void (*expt_check_input_v1)(boost::python::dict params, boost::python::dict default_params, boost::python::list critical_params) = &check_input;
  def("check_input", expt_check_input_v1);  
 
  //std::string int2str(int inp);
  //int find_section(vector<std::string>& A,std::string marker_beg,std::string marker_end,int min_line,int max_line,int& beg,int& end);
  //std::string extract_s(std::string line, std::string marker);    
  
  //void extract_1D(vector<double>& in, vector<double>& out, vector<int>& templ,int shift);
  //void extract_2D(vector< vector<double> >& in, vector< vector<double> >& out, int minx,int maxx, int miny, int maxy );
  //void extract_2D(vector< vector<double> >& in, vector< vector<double> >& out, vector<int>& templ,int shift);



}// export_util_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygutil){
#else
BOOST_PYTHON_MODULE(libutil){
#endif

  export_util_objects();

}



}// namespace libutil
}// liblibra



