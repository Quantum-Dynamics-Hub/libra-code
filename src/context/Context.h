#ifndef CONTEXT_H
#define CONTEXT_H

#include "../io/libio.h"
using namespace libio;

#include "../mmath/libmmath.h"
using namespace libmmath;


namespace libcontext{

class Context{
  
  std::string path;                   // root pathname for the present variable (instance of the Context class)
  boost::property_tree::ptree ctx_pt; // This is the internal representation of the data

  public:

 
  //------------------------------------------------
  Context() {}
  Context(const Context& c){  ctx_pt = c.ctx_pt; } 
  virtual ~Context(){}


  // Slightly more convenient functions
  void add(std::string varname, int varval){   libio::save(ctx_pt, path+varname, varval);  }
  void add(std::string varname, vector<int> varval){   libio::save(ctx_pt, path+varname, varval);  }
  void add(std::string varname, std::string varval){   libio::save(ctx_pt, path+varname, varval);  }
  void add(std::string varname, vector<std::string> varval){   libio::save(ctx_pt, path+varname, varval);  }
  void add(std::string varname, double varval){   libio::save(ctx_pt, path+varname, varval);  }
  void add(std::string varname, vector<double> varval){   libio::save(ctx_pt, path+varname, varval);  }

/*
  void add(std::string varname, VECTOR& varval){   libio::save(ctx_pt, path+varname, varval);  }
  void add(std::string varname, QUATERNION& varval){   libio::save(ctx_pt, path+varname, varval);  }
  void add(std::string varname, MATRIX3x3& varval){   libio::save(ctx_pt, path+varname, varval);  }
  void add(std::string varname, MATRIX& varval){   libio::save(ctx_pt, path+varname, varval);  }
*/
  

  void save_xml(std::string filename){ libio::save_xml(filename, ctx_pt); }
  void load_xml(std::string filename){ libio::load_xml(filename, ctx_pt); }


};

void export_Context_objects();

}// namespace libcontext

#endif // CONTEXT_H
