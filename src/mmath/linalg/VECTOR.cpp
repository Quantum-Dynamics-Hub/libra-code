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

#include "VECTOR.h"

using namespace libio;

namespace libmmath{
namespace liblinalg{


void set_value(int& is_defined, VECTOR& value,boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<VECTOR>(obj.attr(attrName.c_str()));
      is_defined = 1;
  }
}

// ----------- Save --------------
void save(boost::property_tree::ptree& pt,std::string path,VECTOR& vt){
  pt.put(path+".x",vt.x);  pt.put(path+".y",vt.y);  pt.put(path+".z",vt.z);
}

void save(boost::property_tree::ptree& pt,std::string path,vector<VECTOR>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+"."+rt,vt[i]);
  }
}


// ----------- Load --------------
void load(boost::property_tree::ptree& pt,std::string path, VECTOR& vt, int& status){
  status = 0;
  int st;
  libio::load(pt,path+".x",vt.x, st); if(st==1) {status=1;}
  libio::load(pt,path+".y",vt.y, st); if(st==1) {status=1;}
  libio::load(pt,path+".z",vt.z, st); if(st==1) {status=1;}
}

void load(boost::property_tree::ptree& pt,std::string path,vector<VECTOR>& vt,int& status){
  VECTOR x; int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      load(pt,path+"."+v.first,x,st);
      if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}



}// namespace liblinalg
}// libmmath

