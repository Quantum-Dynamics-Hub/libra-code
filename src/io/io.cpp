#include "io.h"


namespace libio{

///=========== Check if Python object has a particularly-named datamember ==================

bool hasattr(boost::python::object obj, std::string attrName) {
     return PyObject_HasAttrString(obj.ptr(), (char*)attrName.c_str());
} 


///=========== For extracting values from Python object to C++ representation ==================

void set_value(int& is_defined, int& value, boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);        
  if(has_attr){
      value = extract<int>(obj.attr(attrName.c_str()));          
      is_defined = 1; 
  }

}

void set_value(int& is_defined, double& value, boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<double>(obj.attr(attrName.c_str()));
      is_defined = 1;
  }

}

void set_value(int& is_defined, std::string& value,boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<std::string>(obj.attr(attrName.c_str()));
      is_defined = 1;
  }
}


void set_list(int& is_defined, vector<int>& value,boost::python::object obj,std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
     boost::python::list lst= extract<boost::python::list>(obj.attr(attrName.c_str()));
     is_defined = 1;

     for(int i=0;i<len(lst);i++){
     int x = extract<int>(lst[i]);
     value.push_back(x);
     }

  }

}


void set_list(int& is_defined, vector<double>& value,boost::python::object obj,std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
     boost::python::list lst= extract<boost::python::list>(obj.attr(attrName.c_str()));
     is_defined = 1;

     for(int i=0;i<len(lst);i++){
     double x = extract<double>(lst[i]);
     value.push_back(x);
     }

  }

}
void set_list(int& is_defined, vector<std::string>& value,boost::python::object obj,std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
     boost::python::list lst= extract<boost::python::list>(obj.attr(attrName.c_str()));
     is_defined = 1;

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

// vector<double>
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

// vector<int>
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

// vector<std::string>
void save(boost::property_tree::ptree& pt,std::string path,vector<std::string>& vt){
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

// vector<double>
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

// vector<int>
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

// vector<std::string>
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



//----------------------------------------------------

void save_xml(std::string filename, boost::property_tree::ptree& pt){

  boost::property_tree::xml_writer_settings<char> settings(' ', 4);
  write_xml(filename, pt, std::locale(), settings);

}

void load_xml(std::string filename, boost::property_tree::ptree& pt){

  read_xml(filename, pt);

}

/*
std::string load_xml(std::string filename, boost::property_tree::ptree& pt){

  read_xml(filename, pt);

}
*/



}// libio
