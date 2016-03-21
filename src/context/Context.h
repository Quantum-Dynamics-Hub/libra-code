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

#ifndef CONTEXT_H
#define CONTEXT_H

#include "../io/libio.h"
using namespace libio;

#include "../mmath/libmmath.h"
using namespace libmmath;


namespace libcontext{

//class Context;

class Context{

  std::string path;  // the top-most level of the property tree = the name of the variable of the "context" type
  char path_separator; //
  boost::property_tree::ptree ctx_pt; // This is the internal representation of the data

  public:

 
  //------------------------------------------------
  Context() { path = "glob_context"; path_separator = '.'; } 
  Context(std::string filename){ 
    path_separator = '.';
    libio::load_xml(filename, ctx_pt);
    int i= 0; BOOST_FOREACH(ptree::value_type &v, ctx_pt){ if(i==0){ path = v.first; } i++;  }
  }
  Context(const Context& c){  ctx_pt = c.ctx_pt; path = c.path; path_separator = c.path_separator; } 
  virtual ~Context(){}

  // Manupulation of the "path": These functions are essentially for getting and setting the name of the context variable (path)
  void set_path(std::string new_path);
  void set_path_separator(char _path_separator){ path_separator = _path_separator; }
  std::string get_path();


  // Add new variables to data-structure
  template <typename X>
  void add(std::string varname, X varval){   libio::save(ctx_pt, path+path_separator+varname, path_separator, varval);  }

  void add_context(Context ctxt);



  void show_children(std::string _path);


  // Get value for given variable name, if exist in datastructure. Or return default value
  template <typename X>
  X get1(std::string varname, X default_val){   
    int st;
    X varval; 

    libio::load(ctx_pt, path+path_separator+varname, path_separator, varval, st); 
    if(st){ return varval; }else{ return default_val; }

  }

  template <typename X>
  X get2(std::string varname, X& default_val){   
    int st;
    X varval; 

    libmmath::liblinalg::load(ctx_pt, path+path_separator+varname, path_separator, varval, st); 
    if(st){ return varval; }else{ return default_val; }

  }



/*
  // Get value for given variable name, if exist in datastructure. Or return default value
  int get(std::string varname,int default_val);
  vector<int> get(std::string varname,vector<int> default_val);

  std::string Context::get(std::string varname,std::string default_val);
  vector<std::string> Context::get(std::string varname,vector<std::string> default_val);

  double get(std::string varname,double default_val);
  vector<double> get(std::string varname,vector<double> default_val);

  VECTOR get(std::string varname,VECTOR default_val);
  vector<VECTOR> get(std::string varname,vector<VECTOR> default_val);

  QUATERNION get(std::string varname,QUATERNION default_val);
  vector<QUATERNION> get(std::string varname,vector<QUATERNION> default_val);

  MATRIX3x3 get(std::string varname,MATRIX3x3 default_val);
  vector<MATRIX3x3> get(std::string varname,vector<MATRIX3x3> default_val);

  MATRIX get(std::string varname,MATRIX default_val);
  vector<MATRIX> get(std::string varname,vector<MATRIX> default_val);
*/

  Context get_context(std::string varname, Context default_val);


  

  void save_xml(std::string filename){ libio::save_xml(filename, ctx_pt); }
  void load_xml(std::string filename){ libio::load_xml(filename, ctx_pt); }


};

void export_Context_objects();

}// namespace libcontext

#endif // CONTEXT_H
