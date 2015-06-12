#include "Context.h"

namespace libcontext{


//-------------- Class methods implementation ------------------------

//void Context::save(boost::property_tree::ptree& pt,std::string path){ ;; }
//void Context::load(boost::property_tree::ptree& pt,std::string path,int& status){ ;; }


//------------------ Export -------------------------

void export_Context_objects(){


  void (Context::*expt_add_v1)(std::string varname, int varval) = &Context::add;
  void (Context::*expt_add_v2)(std::string varname, vector<int> varval) = &Context::add;
  void (Context::*expt_add_v3)(std::string varname, std::string varval) = &Context::add;
  void (Context::*expt_add_v4)(std::string varname, vector<std::string> varval) = &Context::add;
  void (Context::*expt_add_v5)(std::string varname, double varval) = &Context::add;
  void (Context::*expt_add_v6)(std::string varname, vector<double> varval) = &Context::add;


  class_<Context>("Context",init<>())
      .def(init<const Context&>())
//      .def("__copy__", &generic__copy__<Context>)
//      .def("__deepcopy__", &generic__deepcopy__<Context>)

      .def("add",expt_add_v1)
      .def("add",expt_add_v2)
      .def("add",expt_add_v3)
      .def("add",expt_add_v4)
      .def("add",expt_add_v5)
      .def("add",expt_add_v6)

       // Members
      .def("save_xml",&Context::save_xml)
      .def("load_xml",&Context::load_xml)
  ;


}


}// namespace libcontext

