/*********************************************************************************
* Copyright (C) 2016 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libutil.cpp
  \brief The file implements Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libutil.h"

/// libutil namespace
namespace libutil{


void export_util_objects(){
/** 
  \brief Exporter of libutil classes and functions
*/

  void (*expt_show_vector_v1)(vector<int>& A) = &show_vector;  
  def("show_vector", expt_show_vector_v1);
  
  int (*expt_is_in_vector_v1)(int a, vector<int>& A) = &is_in_vector;
  int (*expt_is_in_vector_v2)(int a, vector<int>& A, int& pos) = &is_in_vector;
  int (*expt_is_in_vector_v3)(int a, vector<int>& A, vector<int>& indx) = &is_in_vector2;
  def("is_in_vector", expt_is_in_vector_v1);
  def("is_in_vector", expt_is_in_vector_v2);
  def("is_in_vector", expt_is_in_vector_v3);

  int (*expt_is_repeating_v1)(vector<int>& A,int& reap) = &is_repeating;
  def("is_repeating", expt_is_repeating_v1);  
  
  int (*expt_delta_v1)(vector<int>& A,vector<int>& B,int& a,int& b) = &delta;
  def("delta", expt_delta_v1);  

  void (*expt_split_line_v1)(std::string line, vector<std::string>& arr) = &split_line;
  void (*expt_split_line_v2)(std::string line,vector<std::string>& arr,char delim) = &split_line;
  def("split_line", expt_split_line_v1);  
  def("split_line", expt_split_line_v2);  
 
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




