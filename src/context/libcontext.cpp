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

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libcontext.h"

/// liblibra namespace
namespace liblibra{

namespace libcontext{


//------------------ Export -------------------------

void export_Context_objects(){

  void (Context::*expt_set_path)(std::string new_path) = &Context::set_path;
  void (Context::*expt_set_path_separator)(char _path_separator) = &Context::set_path_separator;

  std::string (Context::*expt_get_path)() = &Context::get_path;


/** 
  Note: the add function defined in the class is a template method. However, when we export this method to Python
  we need to define the pointers to the method to be called with different argument types
 
  Read more about the pointer to a template method of a class at:
  http://goodliffe.blogspot.com/2011/07/c-declaring-pointer-to-template-method.html
 
*/

  void (Context::*expt_add_v1)(std::string varname, int varval) = &Context::add;
  void (Context::*expt_add_v2)(std::string varname, vector<int> varval) = &Context::add;

  void (Context::*expt_add_v3)(std::string varname, std::string varval) = &Context::add;
  void (Context::*expt_add_v4)(std::string varname, vector<std::string> varval) = &Context::add;

  void (Context::*expt_add_v5)(std::string varname, double varval) = &Context::add;
  void (Context::*expt_add_v6)(std::string varname, vector<double> varval) = &Context::add;

  void (Context::*expt_add_v7)(std::string varname, VECTOR varval) = &Context::add;
  void (Context::*expt_add_v8)(std::string varname, vector<VECTOR> varval) = &Context::add;

  void (Context::*expt_add_v9)(std::string varname, QUATERNION varval) = &Context::add;
  void (Context::*expt_add_v10)(std::string varname, vector<QUATERNION> varval) = &Context::add;

  void (Context::*expt_add_v11)(std::string varname, MATRIX3x3 varval) = &Context::add;
  void (Context::*expt_add_v12)(std::string varname, vector<MATRIX3x3> varval) = &Context::add;

  void (Context::*expt_add_v13)(std::string varname, MATRIX varval) = &Context::add;
  void (Context::*expt_add_v14)(std::string varname, vector<MATRIX> varval) = &Context::add;

  void (Context::*expt_add_v15)(Context varval) = &Context::add_context;



  Context (Context::*expt_get_child_v1)(std::string _path, std::string varname, Context default_val) = &Context::get_child;
  Context (Context::*expt_get_child_v2)(std::string varname, Context default_val) = &Context::get_child;

  vector<Context> (Context::*expt_get_children_v1)(std::string _path, std::string varname) = &Context::get_children;
  vector<Context> (Context::*expt_get_children_v2)(std::string varname) = &Context::get_children;
  vector<Context> (Context::*expt_get_children_all_v1)(std::string _path) = &Context::get_children_all;
  vector<Context> (Context::*expt_get_children_all_v2)() = &Context::get_children_all;

  void (Context::*expt_show_children_v1)(std::string _path) = &Context::show_children;
  void (Context::*expt_show_children_v2)() = &Context::show_children;




  int (Context::*expt_get_v1)(std::string varname,int default_val) = &Context::get1;
  vector<int> (Context::*expt_get_v2)(std::string varname,vector<int> default_val) = &Context::get1;

  std::string (Context::*expt_get_v3)(std::string varname,std::string default_val) = &Context::get1;
  vector<std::string> (Context::*expt_get_v4)(std::string varname,vector<std::string> default_val) = &Context::get1;

  double (Context::*expt_get_v5)(std::string varname,double default_val) = &Context::get1;
  vector<double> (Context::*expt_get_v6)(std::string varname,vector<double> default_val) = &Context::get1;


  VECTOR (Context::*expt_get_v7)(std::string varname,VECTOR& default_val) = &Context::get2;
  vector<VECTOR> (Context::*expt_get_v8)(std::string varname,vector<VECTOR>& default_val) = &Context::get2;

  QUATERNION (Context::*expt_get_v9)(std::string varname,QUATERNION& default_val) = &Context::get2;
  vector<QUATERNION> (Context::*expt_get_v10)(std::string varname,vector<QUATERNION>& default_val) = &Context::get2;

  MATRIX3x3 (Context::*expt_get_v11)(std::string varname,MATRIX3x3& default_val) = &Context::get2;
  vector<MATRIX3x3> (Context::*expt_get_v12)(std::string varname,vector<MATRIX3x3>& default_val) = &Context::get2;

  MATRIX (Context::*expt_get_v13)(std::string varname,MATRIX& default_val) = &Context::get2;
  vector<MATRIX> (Context::*expt_get_v14)(std::string varname,vector<MATRIX>& default_val) = &Context::get2;

//  Context (Context::*expt_get_v15)(std::string varname, Context default_val) = &Context::get_child;




//  int (Context::*expt_get_value_v1)() = &Context::get_value;
//  vector<int> (Context::*expt_get_value_v2)() = &Context::get_value;

  std::string (Context::*expt_get_value_v3)() = &Context::get_value;
//  vector<std::string> (Context::*expt_get_value_v4)() = &Context::get_value;

//  double (Context::*expt_get_value_v5)() = &Context::get_value;
//  vector<double> (Context::*expt_get_value_v6)() = &Context::get_value;

/*
  VECTOR (Context::*expt_get_value_v6)() = &Context::get_value;
  vector<VECTOR> (Context::*expt_get_value_v7)() = &Context::get_value;

  QUATERNION (Context::*expt_get_value_v8)() = &Context::get_value;
  vector<QUATERNION> (Context::*expt_get_value_v9)() = &Context::get_value;

  MATRIX (Context::*expt_get_value_v10)() = &Context::get_value;
  vector<MATRIX> (Context::*expt_get_value_v11)() = &Context::get_value;
*/




  class_<Context>("Context",init<>())
      .def(init<std::string>())
      .def(init<const Context&>())
//      .def("__copy__", &generic__copy__<Context>)
//      .def("__deepcopy__", &generic__deepcopy__<Context>)

      .def("set_path",expt_set_path)
      .def("set_path_separator",expt_set_path_separator)
      .def("get_path",expt_get_path)

       /**
        Note the difference: add_v1 - add_v14 are the specializations of the template function
        "add" in C++, but add_v15 is the "add_context" in C++. All types are still exported to 
        Python uder generic "add" name
       */

      .def("add",expt_add_v1)
      .def("add",expt_add_v1)
      .def("add",expt_add_v2)
      .def("add",expt_add_v3)
      .def("add",expt_add_v4)
      .def("add",expt_add_v5)
      .def("add",expt_add_v6)
      .def("add",expt_add_v7)
      .def("add",expt_add_v8)
      .def("add",expt_add_v9)
      .def("add",expt_add_v10)
      .def("add",expt_add_v11)
      .def("add",expt_add_v12)
      .def("add",expt_add_v13)
      .def("add",expt_add_v14)
      .def("add",expt_add_v15)

      .def("get_child",expt_get_child_v1)
      .def("get_child",expt_get_child_v2)
      .def("get_children",expt_get_children_v1)
      .def("get_children",expt_get_children_v2)
      .def("get_children_all",expt_get_children_all_v1)
      .def("get_children_all",expt_get_children_all_v2)
      .def("show_children",expt_show_children_v1)
      .def("show_children",expt_show_children_v2)

      .def("size",&Context::size)


       /**
        Note the difference: get_v1 - get_v14 are the specializations of the template function
        "add" in C++, but get_v15 is the "add_context" in C++. All types are still exported to 
        Python uder generic "add" name
       */

      .def("get",expt_get_v1)
      .def("get",expt_get_v2)
      .def("get",expt_get_v3)
      .def("get",expt_get_v4)
      .def("get",expt_get_v5)
      .def("get",expt_get_v6)
      .def("get",expt_get_v7)
      .def("get",expt_get_v8)
      .def("get",expt_get_v9)
      .def("get",expt_get_v10)
      .def("get",expt_get_v11)
      .def("get",expt_get_v12)
      .def("get",expt_get_v13)
      .def("get",expt_get_v14)
//      .def("get",expt_get_v15)


//      .def("get_value",expt_get_value_v1)
//      .def("get_value",expt_get_value_v2)
      .def("get_value",expt_get_value_v3)
//      .def("get_value",expt_get_value_v4)
//      .def("get_value",expt_get_value_v5)
//      .def("get_value",expt_get_value_v6)
//      .def("get_value",expt_get_value_v7)
//      .def("get_value",expt_get_value_v8)
//      .def("get_value",expt_get_value_v9)
//      .def("get_value",expt_get_value_v10)
//      .def("get_value",expt_get_value_v11)



       // Members
      .def("save_xml",&Context::save_xml)
      .def("load_xml",&Context::load_xml)
  ;


  class_< ContextList >("ContextList")
      .def(vector_indexing_suite< ContextList >())
  ;

  class_< ContextMap >("ContextMap")
      .def(vector_indexing_suite< ContextMap >())
  ;


}



void export_context_objects(){

  export_Context_objects();
  export_ctx_Control_Parameters_objects();

}// export_context_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygcontext){
#else
BOOST_PYTHON_MODULE(libcontext){
#endif

  export_context_objects();

}



}// namespace libcontext
}// liblibra



