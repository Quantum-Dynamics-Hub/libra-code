/*********************************************************************************
* Copyright (C) 2015-2020 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "QUATERNION.h"
#include "VECTOR.h"


/// liblibra namespace
namespace liblibra{

using namespace libio;


/// liblinalg namespace
namespace liblinalg{


VECTOR QUATERNION::vect(){VECTOR v; v.x = Lx; v.y = Ly; v.z = Lz; return v;}

double dot_prod(QUATERNION& q1, QUATERNION& q2){
  return (q1.Lt*q2.Lt + q1.Lx*q2.Lx + q1.Ly*q2.Ly + q1.Lz*q2.Lz);
}
double dot_prod(const QUATERNION& q1, const QUATERNION& q2){
  return (q1.Lt*q2.Lt + q1.Lx*q2.Lx + q1.Ly*q2.Ly + q1.Lz*q2.Lz);
}

QUATERNION operator*(double f,const QUATERNION& q){
  QUATERNION tmp;
  tmp.Lt = f*q.Lt;
  tmp.Lx = f*q.Lx;
  tmp.Ly = f*q.Ly;
  tmp.Lz = f*q.Lz;
  return tmp;
}



void set_value(int& is_defined, QUATERNION& value,boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<QUATERNION>(obj.attr(attrName.c_str()));
      is_defined = 1;
  }
}


// ----------- Save --------------
void save(boost::property_tree::ptree& pt,std::string path,QUATERNION& vt){
  pt.put(path+".Lt",vt.Lt);  pt.put(path+".Lx",vt.Lx);
  pt.put(path+".Ly",vt.Ly);  pt.put(path+".Lz",vt.Lz);
}

void save(boost::property_tree::ptree& pt,std::string path, char path_separator, QUATERNION& vt){

  pt.put(boost::property_tree::ptree::path_type(path+std::string(1,path_separator)+"Lt", path_separator),vt.Lt);
  pt.put(boost::property_tree::ptree::path_type(path+std::string(1,path_separator)+"Lx", path_separator),vt.Lx);
  pt.put(boost::property_tree::ptree::path_type(path+std::string(1,path_separator)+"Ly", path_separator),vt.Ly);
  pt.put(boost::property_tree::ptree::path_type(path+std::string(1,path_separator)+"Lz", path_separator),vt.Lz);
}


void save(boost::property_tree::ptree& pt,std::string path,vector<QUATERNION>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+"."+rt,vt[i]);
  }
}

void save(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<QUATERNION>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+std::string(1,path_separator)+rt, path_separator, vt[i]);
  }
}



// ----------- Load --------------
void load(boost::property_tree::ptree& pt,std::string path, QUATERNION& vt, int& status){
  status = 0;
  int st;
  libio::load(pt,path+".Lt",vt.Lt, st); if(st==1) {status=1;}
  libio::load(pt,path+".Lx",vt.Lx, st); if(st==1) {status=1;}
  libio::load(pt,path+".Ly",vt.Ly, st); if(st==1) {status=1;}
  libio::load(pt,path+".Lz",vt.Lz, st); if(st==1) {status=1;}
}

void load(boost::property_tree::ptree& pt,std::string path, char path_separator, QUATERNION& vt, int& status){
  status = 0;
  int st;
  libio::load(pt,path+std::string(1,path_separator)+"Lt", path_separator, vt.Lt, st); if(st==1) {status=1;}
  libio::load(pt,path+std::string(1,path_separator)+"Lx", path_separator, vt.Lx, st); if(st==1) {status=1;}
  libio::load(pt,path+std::string(1,path_separator)+"Ly", path_separator, vt.Ly, st); if(st==1) {status=1;}
  libio::load(pt,path+std::string(1,path_separator)+"Lz", path_separator, vt.Lz, st); if(st==1) {status=1;}
}


void load(boost::property_tree::ptree& pt,std::string path,vector<QUATERNION>& vt,int& status){
  QUATERNION x; int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      load(pt,path+"."+v.first,x,st);
      if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }

}

void load(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<QUATERNION>& vt,int& status){
  QUATERNION x; int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(boost::property_tree::ptree::path_type(path,path_separator))){
      load(pt,path+std::string(1,path_separator)+v.first, path_separator, x, st);
      if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }

}


}// namespace liblinalg
}// namespace liblibra

