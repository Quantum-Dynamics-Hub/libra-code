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

#include "Context.h"


namespace libcontext{


//-------------- Class methods implementation ------------------------

void Context::set_path(std::string new_path){
  path = new_path;
  int i= 0;
  BOOST_FOREACH(ptree::value_type &v, ctx_pt){ 
    if(i==0){ v.first = new_path; } i++;  
  } 

}
std::string Context::get_path(){  return path; }


//------------------ Add functions ----------------------
//-------------------------------------------------------

// I know, using templates would be more logical, but ::save function would need to be defined as template, but it is specialized
void Context::add(std::string varname, int varval){   libio::save(ctx_pt, path+"."+varname, varval);  }
void Context::add(std::string varname, vector<int> varval){   libio::save(ctx_pt, path+"."+varname, varval);  }

void Context::add(std::string varname, std::string varval){   libio::save(ctx_pt, path+"."+varname, varval);  }
void Context::add(std::string varname, vector<std::string> varval){   libio::save(ctx_pt, path+"."+varname, varval);  }

void Context::add(std::string varname, double varval){   libio::save(ctx_pt, path+"."+varname, varval);  }
void Context::add(std::string varname, vector<double> varval){   libio::save(ctx_pt, path+"."+varname, varval);  }

void Context::add(std::string varname, VECTOR varval){   libmmath::liblinalg::save(ctx_pt, path+"."+varname, varval);  }
void Context::add(std::string varname, vector<VECTOR> varval){   libmmath::liblinalg::save(ctx_pt, path+"."+varname, varval);  }

void Context::add(std::string varname, QUATERNION varval){   libmmath::liblinalg::save(ctx_pt, path+"."+varname, varval);  }
void Context::add(std::string varname, vector<QUATERNION> varval){   libmmath::liblinalg::save(ctx_pt, path+"."+varname, varval);  }

void Context::add(std::string varname, MATRIX3x3 varval){   libmmath::liblinalg::save(ctx_pt, path+"."+varname, varval);  }
void Context::add(std::string varname, vector<MATRIX3x3> varval){   libmmath::liblinalg::save(ctx_pt, path+"."+varname, varval);  }

void Context::add(std::string varname, MATRIX varval){   libmmath::liblinalg::save(ctx_pt, path+"."+varname, varval);  }
void Context::add(std::string varname, vector<MATRIX> varval){   libmmath::liblinalg::save(ctx_pt, path+"."+varname, varval);  }

void Context::add(Context ctxt){

  int i= 0;
  BOOST_FOREACH(ptree::value_type &v, ctx_pt){ 
    if(i==0){ 
      int j = 0;
      BOOST_FOREACH(ptree::value_type &v1, ctxt.ctx_pt){
        if(j==0){ v.second.put_child(ctxt.path, v1.second); } j++;
      }
    } i++;  
  } // foreach
}// add


//------------------ Get functions ----------------------
//-------------------------------------------------------

int Context::get(std::string varname,int default_val){ 
  int st;
  int varval; 

  libio::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

vector<int> Context::get(std::string varname,vector<int> default_val){ 
  int st;
  vector<int> varval; 

  libio::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

std::string Context::get(std::string varname,std::string default_val){ 
  int st;
  std::string varval; 

  libio::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

vector<std::string> Context::get(std::string varname,vector<std::string> default_val){ 
  int st;
  vector<std::string> varval; 

  libio::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

double Context::get(std::string varname,double default_val){ 
  int st;
  double varval; 

  libio::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

vector<double> Context::get(std::string varname,vector<double> default_val){ 
  int st;
  vector<double> varval; 

  libio::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

VECTOR Context::get(std::string varname,VECTOR default_val){
  int st;
  VECTOR varval; 

  libmmath::liblinalg::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

vector<VECTOR> Context::get(std::string varname,vector<VECTOR> default_val){
  int st;
  vector<VECTOR> varval; 

  libmmath::liblinalg::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

QUATERNION Context::get(std::string varname,QUATERNION default_val){
  int st;
  QUATERNION varval; 

  libmmath::liblinalg::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

vector<QUATERNION> Context::get(std::string varname,vector<QUATERNION> default_val){
  int st;
  vector<QUATERNION> varval; 

  libmmath::liblinalg::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}


MATRIX3x3 Context::get(std::string varname,MATRIX3x3 default_val){
  int st;
  MATRIX3x3 varval; 

  libmmath::liblinalg::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

vector<MATRIX3x3> Context::get(std::string varname,vector<MATRIX3x3> default_val){
  int st;
  vector<MATRIX3x3> varval; 

  libmmath::liblinalg::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}


MATRIX Context::get(std::string varname,MATRIX default_val){
  int st;
  MATRIX varval; 

  libmmath::liblinalg::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

vector<MATRIX> Context::get(std::string varname,vector<MATRIX> default_val){
  int st;
  vector<MATRIX> varval; 

  libmmath::liblinalg::load(ctx_pt, path+"."+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}




Context Context::get(std::string varname, Context default_val){  // here, we actually ask for path name
// This is gonna be a recursive function

  Context res;
  boost::property_tree::ptree x;

  int status = 1;
  try{ 
    x = ctx_pt.get_child(path+"."+varname);
  } catch(std::exception& e){ status = 0; }

   
  if(status){
    res.path = varname;
    BOOST_FOREACH(ptree::value_type &v, x){ res.ctx_pt.put_child(varname+"."+v.first, v.second);  } 
  }
  else{
    res = Context(default_val);
  }
 

  return res;
}


//------------------ Export -------------------------

void export_Context_objects(){

  void (Context::*expt_set_path)(std::string new_path) = &Context::set_path;
  std::string (Context::*expt_get_path)() = &Context::get_path;

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

  void (Context::*expt_add_v15)(Context varval) = &Context::add;




  int (Context::*expt_get_v1)(std::string varname,int default_val) = &Context::get;
  vector<int> (Context::*expt_get_v2)(std::string varname,vector<int> default_val) = &Context::get;

  std::string (Context::*expt_get_v3)(std::string varname,std::string default_val) = &Context::get;
  vector<std::string> (Context::*expt_get_v4)(std::string varname,vector<std::string> default_val) = &Context::get;

  double (Context::*expt_get_v5)(std::string varname,double default_val) = &Context::get;
  vector<double> (Context::*expt_get_v6)(std::string varname,vector<double> default_val) = &Context::get;

  VECTOR (Context::*expt_get_v7)(std::string varname,VECTOR default_val) = &Context::get;
  vector<VECTOR> (Context::*expt_get_v8)(std::string varname,vector<VECTOR> default_val) = &Context::get;

  QUATERNION (Context::*expt_get_v9)(std::string varname,QUATERNION default_val) = &Context::get;
  vector<QUATERNION> (Context::*expt_get_v10)(std::string varname,vector<QUATERNION> default_val) = &Context::get;

  MATRIX3x3 (Context::*expt_get_v11)(std::string varname,MATRIX3x3 default_val) = &Context::get;
  vector<MATRIX3x3> (Context::*expt_get_v12)(std::string varname,vector<MATRIX3x3> default_val) = &Context::get;

  MATRIX (Context::*expt_get_v13)(std::string varname,MATRIX default_val) = &Context::get;
  vector<MATRIX> (Context::*expt_get_v14)(std::string varname,vector<MATRIX> default_val) = &Context::get;

  Context (Context::*expt_get_v15)(std::string varname, Context default_val) = &Context::get;





  class_<Context>("Context",init<>())
      .def(init<std::string>())
      .def(init<const Context&>())
//      .def("__copy__", &generic__copy__<Context>)
//      .def("__deepcopy__", &generic__deepcopy__<Context>)

      .def("set_path",expt_set_path)
      .def("get_path",expt_get_path)
      
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
      .def("get",expt_get_v15)



       // Members
      .def("save_xml",&Context::save_xml)
      .def("load_xml",&Context::load_xml)
  ;


}


}// namespace libcontext

