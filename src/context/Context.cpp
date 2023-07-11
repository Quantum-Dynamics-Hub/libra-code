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

#include "Context.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

namespace libcontext{


//-------------- Class methods implementation ------------------------

void Context::set_path(std::string new_path){
  path = new_path;
  int i= 0;
  BOOST_FOREACH(ptree::value_type& v, ctx_pt){ 
  //BOOST_FOREACH(auto& v, ctx_pt){ 
  //for (auto& v : ctx_pt){  // C++11
   // if(i==0){ v.first = std::move(new_path); } i++;    AVA: Temporary comment it to be able to compile with C++11
  } 

}
std::string Context::get_path(){ 
  return  path;
}

int Context::size(){ 
/**
  The number of direct children on this node
*/

  return ctx_pt.size();
}



//------------------ Add functions ----------------------
//-------------------------------------------------------

void Context::add_context(Context ctxt){
/**
  Copies one context object into another
*/

  int i= 0;
  BOOST_FOREACH(ptree::value_type &v, ctx_pt){ 
    if(i==0){ 
      int j = 0;
      BOOST_FOREACH(ptree::value_type &v1, ctxt.ctx_pt){
        if(j==0){ v.second.put_child(boost::property_tree::ptree::path_type(ctxt.path, ctxt.path_separator), v1.second); } j++;
      }
    } i++;  
  } // foreach
}// add


void Context::show_children(std::string _path){
/**
  Shows all the children of the property tree at the required level 

  BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){ 
  pt.get_child(path) - array of property trees at given level, path
  v - each element (property tree, value type) in this array 
  v.first - the path to this element
  v.secod - the value of the element

*/


  int j = 0;
  cout<<"current path = "<<path<<endl;

  BOOST_FOREACH(ptree::value_type &v1, ctx_pt.get_child(boost::property_tree::ptree::path_type(_path, path_separator))){
    cout<<"key = "<<v1.first<<endl; //
  }

}// add


void Context::show_children(){
/**
  Shows all the children of the property tree at the present level

*/
  show_children(path);

}// add






//------------------ Get functions ----------------------
//-------------------------------------------------------


Context Context::get_child(std::string _path, std::string varname, Context default_val){ 


  boost::property_tree::ptree x;

  x = ctx_pt.get_child(boost::property_tree::ptree::path_type(_path, path_separator));

  BOOST_FOREACH(ptree::value_type &v, x){ 

    if(v.first == varname){
      Context res;
      res.path_separator = path_separator;
      res.path = v.first;
      res.ctx_pt.add_child(boost::property_tree::ptree::path_type(v.first, path_separator), v.second);  
      return res;
    }
  } 

  return Context(default_val);
}


Context Context::get_child(std::string varname, Context default_val){ 

  return get_child(path, varname, default_val);

}


vector<Context> Context::get_children(std::string _path, std::string varname){
/**
   Looks for all childeren named "varname" in the path given by "_path"
*/

  vector<Context> res;

  boost::property_tree::ptree x;

  x = ctx_pt.get_child(boost::property_tree::ptree::path_type(_path, path_separator));

  BOOST_FOREACH(ptree::value_type &v1, x){

    if(v1.first==varname){

      Context y; 
      y.path_separator = path_separator;
      y.path = v1.first;
      y.ctx_pt.add_child(boost::property_tree::ptree::path_type(v1.first, path_separator), v1.second);  
      res.push_back( y );
    }
  }

  return res;

}// add


vector<Context> Context::get_children(std::string varname){
/**
   Looks for all childeren named "varname" in the current path
*/

  return get_children(path, varname);

}// add



vector<Context> Context::get_children_all(std::string _path){
/**
   Looks for all childeren named "varname" in the path given by "_path"
*/

  vector<Context> res;

  boost::property_tree::ptree x;

  x = ctx_pt.get_child(boost::property_tree::ptree::path_type(_path, path_separator));

  BOOST_FOREACH(ptree::value_type &v1, x){

    Context y; 
    y.path_separator = path_separator;
    y.path = v1.first; 
    y.ctx_pt.add_child(boost::property_tree::ptree::path_type(v1.first, path_separator), v1.second);  
    res.push_back( y );

  }

  return res;

}// add


vector<Context> Context::get_children_all(){
/**
   Looks for all childeren named in the current path
*/

  return get_children_all(path);

}// add






}// namespace libcontext
}// liblibra

