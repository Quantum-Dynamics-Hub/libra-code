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
void Context::add(std::string varname, int varval){   libio::save(ctx_pt, path+varname, varval);  }
void Context::add(std::string varname, vector<int> varval){   libio::save(ctx_pt, path+varname, varval);  }
void Context::add(std::string varname, std::string varval){   libio::save(ctx_pt, path+varname, varval);  }
void Context::add(std::string varname, vector<std::string> varval){   libio::save(ctx_pt, path+varname, varval);  }
void Context::add(std::string varname, double varval){   libio::save(ctx_pt, path+varname, varval);  }
void Context::add(std::string varname, vector<double> varval){   libio::save(ctx_pt, path+varname, varval);  }

void Context::add(Context ctxt){

  int i= 0;
  BOOST_FOREACH(ptree::value_type &v, ctx_pt){ 
    if(i==0){ v.second.add_child(ctxt.path, ctxt.ctx_pt); } i++;  
  } 
   
//  ctx_pt.add_child(ctxt.path, ctxt.ctx_pt);
}


//------------------ Get functions ----------------------
//-------------------------------------------------------

int Context::get(std::string varname,int default_val){ 
  int st;
  int varval; 

  libio::load(ctx_pt, path+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

vector<int> Context::get(std::string varname,vector<int> default_val){ 
  int st;
  vector<int> varval; 

  libio::load(ctx_pt, path+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

std::string Context::get(std::string varname,std::string default_val){ 
  int st;
  std::string varval; 

  libio::load(ctx_pt, path+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

vector<std::string> Context::get(std::string varname,vector<std::string> default_val){ 
  int st;
  vector<std::string> varval; 

  libio::load(ctx_pt, path+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

double Context::get(std::string varname,double default_val){ 
  int st;
  double varval; 

  libio::load(ctx_pt, path+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
}

vector<double> Context::get(std::string varname,vector<double> default_val){ 
  int st;
  vector<double> varval; 

  libio::load(ctx_pt, path+varname, varval, st); 
  if(st){ return varval; }else{ return default_val; }
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
  void (Context::*expt_add_v7)(Context varval) = &Context::add;

  int (Context::*expt_get_v1)(std::string varname,int default_val) = &Context::get;
  vector<int> (Context::*expt_get_v2)(std::string varname,vector<int> default_val) = &Context::get;
  std::string (Context::*expt_get_v3)(std::string varname,std::string default_val) = &Context::get;
  vector<std::string> (Context::*expt_get_v4)(std::string varname,vector<std::string> default_val) = &Context::get;
  double (Context::*expt_get_v5)(std::string varname,double default_val) = &Context::get;
  vector<double> (Context::*expt_get_v6)(std::string varname,vector<double> default_val) = &Context::get;





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

      .def("get",expt_get_v1)
      .def("get",expt_get_v2)
      .def("get",expt_get_v3)
      .def("get",expt_get_v4)
      .def("get",expt_get_v5)
      .def("get",expt_get_v6)


       // Members
      .def("save_xml",&Context::save_xml)
      .def("load_xml",&Context::load_xml)
  ;


}


}// namespace libcontext

