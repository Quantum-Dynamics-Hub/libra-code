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
/**
  \file VECTOR.h
  \brief The file describes the VECTOR class for representing 3D point
    
*/

#include "VECTOR.h"

/// liblibra namespace
namespace liblibra{

using namespace std;
using namespace libio;


/// liblinalg namespace
namespace liblinalg{



VECTOR cross(const double k, const VECTOR &v1, const VECTOR &v2){
/**
  This function returns the cross-product of two vectors times the scalar
*/

  VECTOR tmp;
  tmp.x=k*(v1.y*v2.z-v2.y*v1.z);
  tmp.y=k*(v1.z*v2.x-v2.z*v1.x);
  tmp.z=k*(v1.x*v2.y-v2.x*v1.y);
  return tmp;
}



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
/**
  Save the state of the VECTOR object as a property tree

  x, y, and z components are added as nodes to the property tree. The nodes are added to 
  the level of the tree controlled by the path variable.
 
  \param[in,out] pt The property tree to which the properties of the VECTOR are added
  \param[in] path The parameter controlling the level of the tree to which the VECTOR members will be added.
  \param[in] vt The VECTOR object which we want to add to the property tree
*/

  pt.put(path+".x",vt.x);  pt.put(path+".y",vt.y);  pt.put(path+".z",vt.z);
}

void save(boost::property_tree::ptree& pt,std::string path, char path_separator, VECTOR& vt){
/**
  Save the state of the VECTOR object as a property tree

  x, y, and z components are added as nodes to the property tree. The nodes are added to 
  the level of the tree controlled by the path variable.
 
  \param[in,out] pt The property tree to which the properties of the VECTOR are added
  \param[in] path The parameter controlling the level of the tree to which the VECTOR members will be added.
  \param[in] vt The VECTOR object which we want to add to the property tree
*/

  pt.put(boost::property_tree::ptree::path_type(path+std::string(1,path_separator)+"x", path_separator),vt.x);
  pt.put(boost::property_tree::ptree::path_type(path+std::string(1,path_separator)+"y", path_separator),vt.y);
  pt.put(boost::property_tree::ptree::path_type(path+std::string(1,path_separator)+"z", path_separator),vt.z);

}


void save(boost::property_tree::ptree& pt,std::string path,vector<VECTOR>& vt){
/**
  Save the state of the vector of VECTOR objects as a property tree

  Each VECTOR object is added as a separate branch. 
 
  \param[in,out] pt The property tree to which the list of the VECTOR objects will be added
  \param[in] path The parameter controlling the level of the tree to which the list of VECTOR will be added.
  \param[in] vt The list of VECTOR objects to be printed out into property tree
*/

  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+"."+rt,vt[i]);
  }
}

void save(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<VECTOR>& vt){
/**
  Save the state of the vector of VECTOR objects as a property tree

  Each VECTOR object is added as a separate branch. 
 
  \param[in,out] pt The property tree to which the list of the VECTOR objects will be added
  \param[in] path The parameter controlling the level of the tree to which the list of VECTOR will be added.
  \param[in] vt The list of VECTOR objects to be printed out into property tree
*/

  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+std::string(1,path_separator)+rt, path_separator,vt[i]);
  }
}



// ----------- Load --------------
void load(boost::property_tree::ptree& pt,std::string path, VECTOR& vt, int& status){
/**
  Load the state of the VECTOR object from a property tree

  Each data member found in the property tree is extracted as the member of the VECTOR object (only x, y, z members). The
  status of each found data member is set to 1.
 
  \param[in] pt The property tree from which the properties of the VECTOR object will be extracted
  \param[in] path The parameter controlling from which level of the tree we try to extract the VECTOR object
  \param[out] status Is the global status of the success of the operation. It is 1 if at least one VECTOR member is found at
              given level of the property tree.
*/

  status = 0;
  int st;
  libio::load(pt,path+".x",vt.x, st); if(st==1) {status=1;}
  libio::load(pt,path+".y",vt.y, st); if(st==1) {status=1;}
  libio::load(pt,path+".z",vt.z, st); if(st==1) {status=1;}
}

void load(boost::property_tree::ptree& pt,std::string path, char path_separator, VECTOR& vt, int& status){
/**
  Load the state of the VECTOR object from a property tree

  Each data member found in the property tree is extracted as the member of the VECTOR object (only x, y, z members). The
  status of each found data member is set to 1.
 
  \param[in] pt The property tree from which the properties of the VECTOR object will be extracted
  \param[in] path The parameter controlling from which level of the tree we try to extract the VECTOR object
  \param[out] status Is the global status of the success of the operation. It is 1 if at least one VECTOR member is found at
              given level of the property tree.
*/

  status = 0;
  int st;
  libio::load(pt,path+std::string(1,path_separator)+"x", path_separator,vt.x, st); if(st==1) {status=1;}
  libio::load(pt,path+std::string(1,path_separator)+"y", path_separator,vt.y, st); if(st==1) {status=1;}
  libio::load(pt,path+std::string(1,path_separator)+"z", path_separator,vt.z, st); if(st==1) {status=1;}
}


void load(boost::property_tree::ptree& pt,std::string path,vector<VECTOR>& vt,int& status){
/**
  Load the vector of VECTOR objects from a property tree

  Each VECTOR object is extracted from a separate branch. 
 
  \param[in] pt The property tree from which the vector of VECTOR objects will be extracted
  \param[in] path The parameter controlling from which level of the property tree we will try to extract the vector of VECTOR objects
  \param[out] vt The vector of created VECTOR objects
  \param[out] status Is the global status of the success of the operation. It is 1 if at least one VECTOR object is extracted
*/

  VECTOR x; int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      load(pt,path+"."+v.first,x,st);
      if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}

void load(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<VECTOR>& vt,int& status){
/**
  Load the vector of VECTOR objects from a property tree

  Each VECTOR object is extracted from a separate branch. 
 
  \param[in] pt The property tree from which the vector of VECTOR objects will be extracted
  \param[in] path The parameter controlling from which level of the property tree we will try to extract the vector of VECTOR objects
  \param[out] vt The vector of created VECTOR objects
  \param[out] status Is the global status of the success of the operation. It is 1 if at least one VECTOR object is extracted
*/

  VECTOR x; int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(boost::property_tree::ptree::path_type(path, path_separator))){
      load(pt,path+std::string(1,path_separator)+v.first, path_separator, x,st);
      if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}





}// liblinalg

}// liblibra

