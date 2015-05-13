#include "Utility.h"

using namespace boost::python;


namespace libmmath{
namespace liblinalg{


bool hasattr(boost::python::object obj, std::string attrName) {
     return PyObject_HasAttrString(obj.ptr(), (char*)attrName.c_str());
} 


void set_value(int& defined, int& value, boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);        
  if(has_attr){
      value = extract<int>(obj.attr(attrName.c_str()));          
      defined = 1; 
  }

}

void set_value(int& defined, double& value, boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<double>(obj.attr(attrName.c_str()));
      defined = 1;
  }

}

void set_value(int& defined, std::string& value,boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<std::string>(obj.attr(attrName.c_str()));
      defined = 1;
  }
}

void set_value(int& defined, VECTOR& value,boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<VECTOR>(obj.attr(attrName.c_str()));
      defined = 1;
  }
}

void set_value(int& defined, MATRIX& value,boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<MATRIX>(obj.attr(attrName.c_str()));
      defined = 1;
  }
}

void set_value(int& defined, QUATERNION& value,boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<QUATERNION>(obj.attr(attrName.c_str()));
      defined = 1;
  }
}




void set_list(int& defined, vector<int>& value,boost::python::object obj,std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
     boost::python::list lst= extract<boost::python::list>(obj.attr(attrName.c_str()));
     defined = 1;

     for(int i=0;i<len(lst);i++){
     int x = extract<int>(lst[i]);
     value.push_back(x);
     }

  }

}
void set_list(int& defined, vector<double>& value,boost::python::object obj,std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
     boost::python::list lst= extract<boost::python::list>(obj.attr(attrName.c_str()));
     defined = 1;

     for(int i=0;i<len(lst);i++){
     double x = extract<double>(lst[i]);
     value.push_back(x);
     }

  }

}
void set_list(int& defined, vector<std::string>& value,boost::python::object obj,std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
     boost::python::list lst= extract<boost::python::list>(obj.attr(attrName.c_str()));
     defined = 1;

     for(int i=0;i<len(lst);i++){
     std::string x = extract<std::string>(lst[i]);
     value.push_back(x);
     }

  }

}

//----------------- XML via property_tree ---------------------------
// Overload for basic types
// double
void save(boost::property_tree::ptree& pt,std::string path,double& vt ){
  pt.put(path,vt);
}

void save(boost::property_tree::ptree& pt,std::string path,vector<double>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+"."+rt,vt[i]);
  }
}

// int
void save(boost::property_tree::ptree& pt,std::string path,int& vt ){
  pt.put(path,vt);
}

void save(boost::property_tree::ptree& pt,std::string path,vector<int>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+"."+rt,vt[i]);
  }
}

// std::string
void save(boost::property_tree::ptree& pt,std::string path,std::string& vt ){
  pt.put(path,vt);
}

void save(boost::property_tree::ptree& pt,std::string path,vector<std::string>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+"."+rt,vt[i]);
  }
}


// Overloaded for VECTOR
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

// Overloaded for QUATERNION
void save(boost::property_tree::ptree& pt,std::string path,QUATERNION& vt){
  pt.put(path+".Lt",vt.Lt);  pt.put(path+".Lx",vt.Lx);
  pt.put(path+".Ly",vt.Ly);  pt.put(path+".Lz",vt.Lz);
}

void save(boost::property_tree::ptree& pt,std::string path,vector<QUATERNION>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+"."+rt,vt[i]);
  }
}

// Overloaded for MATRIX3x3
void save(boost::property_tree::ptree& pt,std::string path,MATRIX3x3& vt){
  pt.put(path+".xx",vt.xx);  pt.put(path+".xy",vt.xy);  pt.put(path+".xz",vt.xz);
  pt.put(path+".yx",vt.yx);  pt.put(path+".yy",vt.yy);  pt.put(path+".yz",vt.yz);
  pt.put(path+".zx",vt.zx);  pt.put(path+".zy",vt.zy);  pt.put(path+".zz",vt.zz);
}

void save(boost::property_tree::ptree& pt,std::string path,vector<MATRIX3x3>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+"."+rt,vt[i]);
  }
}


// Overloaded for MATRIX
void save(boost::property_tree::ptree& pt,std::string path,MATRIX& vt){
  for(int i=0;i<vt.num_of_elems;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    pt.put(path+"."+rt,vt.M[i]);
  }
}

void save(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+"."+rt,vt[i]);
  }
}


//--------------------------------------------------------
// Overload for basic types
// double
void load(boost::property_tree::ptree& pt,std::string path,double& vt, int& status){
  status = 1;
  try{ vt = pt.get<double>(path); } catch(std::exception& e){ status = 0; }
}

void load(boost::property_tree::ptree& pt,std::string path,vector<double>& vt,int& status){
  double x; int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      load(pt,path+"."+v.first,x,st);
      if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}


// int
void load(boost::property_tree::ptree& pt,std::string path,int& vt, int& status){
  status = 1;
  try{ vt = pt.get<int>(path); } catch(std::exception& e){ status = 0; }
}

void load(boost::property_tree::ptree& pt,std::string path,vector<int>& vt,int& status){
  int x; int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      load(pt,path+"."+v.first,x,st);
      if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}


// std::string
void load(boost::property_tree::ptree& pt,std::string path,std::string& vt, int& status){
  status = 1;
  try{ vt = pt.get<std::string>(path); } catch(std::exception& e){ status = 0; }
}

void load(boost::property_tree::ptree& pt,std::string path,vector<std::string>& vt,int& status){
  std::string x; int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      load(pt,path+"."+v.first,x,st);
      if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}


// VECTOR
void load(boost::property_tree::ptree& pt,std::string path, VECTOR& vt, int& status){
  status = 0;
  int st;
  load(pt,path+".x",vt.x, st); if(st==1) {status=1;}
  load(pt,path+".y",vt.y, st); if(st==1) {status=1;}
  load(pt,path+".z",vt.z, st); if(st==1) {status=1;}
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


// QUATERNION
void load(boost::property_tree::ptree& pt,std::string path, QUATERNION& vt, int& status){
  status = 0;
  int st;
  load(pt,path+".Lt",vt.Lt, st); if(st==1) {status=1;}
  load(pt,path+".Lx",vt.Lx, st); if(st==1) {status=1;}
  load(pt,path+".Ly",vt.Ly, st); if(st==1) {status=1;}
  load(pt,path+".Lz",vt.Lz, st); if(st==1) {status=1;}
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

// MATRIX3x3
void load(boost::property_tree::ptree& pt,std::string path, MATRIX3x3& vt, int& status){
  status = 0;
  int st;
  load(pt,path+".xx",vt.xx, st); if(st==1) {status=1;}
  load(pt,path+".xy",vt.xy, st); if(st==1) {status=1;}
  load(pt,path+".xz",vt.xz, st); if(st==1) {status=1;}

  load(pt,path+".yx",vt.yx, st); if(st==1) {status=1;}
  load(pt,path+".yy",vt.yy, st); if(st==1) {status=1;}
  load(pt,path+".yz",vt.yz, st); if(st==1) {status=1;}

  load(pt,path+".zx",vt.zx, st); if(st==1) {status=1;}
  load(pt,path+".zy",vt.zy, st); if(st==1) {status=1;}
  load(pt,path+".zz",vt.zz, st); if(st==1) {status=1;}

}
void load(boost::property_tree::ptree& pt,std::string path,vector<MATRIX3x3>& vt,int& status){
  MATRIX3x3 x; int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      load(pt,path+"."+v.first,x,st);
      if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }


}

// MATRIX
void load(boost::property_tree::ptree& pt,std::string path, MATRIX& vt, int& status){ }
void load(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt,int& status){ }


}// namespace liblinalg
}// libmmath
