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
/**
  \file io.h
  \brief The file describes auxiliary (but basic) functions for Python input/output operations
    
*/

#ifndef IO_H
#define IO_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <string>
#include <vector>
#include <map>
#include <set>
#include <exception>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <boost/python.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#endif 



/// liblibra 
namespace liblibra{


using namespace std;
using namespace boost::python;
using boost::property_tree::ptree;

// The signature of boost::property_tree::xml_parser::write_xml() changed in Boost 1.56
// See https://github.com/PointCloudLibrary/pcl/issues/864
#include <boost/version.hpp>
#if (BOOST_VERSION >= 105600)
  typedef boost::property_tree::xml_writer_settings<std::string> xml_writer_settings;
#else
  typedef boost::property_tree::xml_writer_settings<char> xml_writer_settings;
#endif




/// libio namespace
namespace libio{


bool hasattr(boost::python::object obj, std::string attrName);

// For extracting values from Python object to C++ representation
void set_value(int& is_defined, int& value,         boost::python::object obj, std::string attrName);
void set_value(int& is_defined, double& value,      boost::python::object obj, std::string attrName);
void set_value(int& is_defined, std::string& value, boost::python::object obj, std::string attrName);
void set_list(int& is_defined, vector<int>& value,        boost::python::object obj,std::string attrName);
void set_list(int& is_defined, vector<double>& value,     boost::python::object obj,std::string attrName);
void set_list(int& is_defined, vector<std::string>& value,boost::python::object obj,std::string attrName);




//------------------------ XML via property_tree ---------------------------
/*
template<typename X> void save(boost::property_tree::ptree& pt,std::string path, char path_separator, X& vt );
template<typename X> void save(boost::property_tree::ptree& pt,std::string path, X& vt );

template<typename X> void save(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<X>& vt );
template<typename X> void save(boost::property_tree::ptree& pt,std::string path, vector<X>& vt );



template<typename X> void load(boost::property_tree::ptree& pt,std::string path, char path_separator, X& vt, int& status);
template<typename X> void load(boost::property_tree::ptree& pt,std::string path, X& vt, int& status);

template<typename X> void load(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<X>& vt,int& status);
template<typename X> void load(boost::property_tree::ptree& pt,std::string path,vector<X>& vt,int& status);

*/


template<typename X>
void save(boost::property_tree::ptree& pt,std::string path, char path_separator, X& vt ){
/** 
  \brief Saves double in the property tree

  \param[in, out] pt Is the property tree to which we add the variable
  \param[in] path Is the path of the added variable in the property tree. This is essentially an extended name of the added variable.
  \param[in] vt The double value of the variable we are adding to the property tree.
*/
  pt.put(boost::property_tree::ptree::path_type(path, path_separator),vt);
}

template<typename X>
void save(boost::property_tree::ptree& pt,std::string path, X& vt ){
/** 
  \brief Saves double in the property tree

  \param[in, out] pt Is the property tree to which we add the variable
  \param[in] path Is the path of the added variable in the property tree. This is essentially an extended name of the added variable.
  \param[in] vt The double value of the variable we are adding to the property tree.
*/
  pt.put(path, vt);
}


template<typename X>
void save(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<X>& vt ){
/** 
  \brief Saves double in the property tree

  \param[in, out] pt Is the property tree to which we add the variable
  \param[in] path Is the path of the added variable in the property tree. This is essentially an extended name of the added variable.
  \param[in] vt The double value of the variable we are adding to the property tree.
*/

  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save<X>(pt,path+std::string(1,path_separator)+rt, path_separator,vt[i]);
  }

}

template<typename X>
void save(boost::property_tree::ptree& pt,std::string path, vector<X>& vt ){
/** 
  \brief Saves double in the property tree

  \param[in, out] pt Is the property tree to which we add the variable
  \param[in] path Is the path of the added variable in the property tree. This is essentially an extended name of the added variable.
  \param[in] vt The double value of the variable we are adding to the property tree.
*/

  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save<X>(pt,path+"."+rt, vt[i]);
  }

}




template<typename X>
void load(boost::property_tree::ptree& pt,std::string path, X& vt, int& status){
/** 
  \brief Extracts the double value from the property tree

  \param[in] pt Is the property tree from which we want to extract the value
  \param[in] path Is the path of the variable to extract from the property tree. 
                  This is essentially an extended name of the variable to extract.
  \param[out] vt The double value of the variable we are extrating from the property tree.
  \param[out] status Is the flag showing the sucess of the extraction: 0 - is no such variable name was found, 
           1 - if we have completed the extraction sucessfully.
*/

  status = 1;
  try{ vt = pt.get<X>(path); } catch(std::exception& e){ status = 0; }
}

template<typename X>
void load(boost::property_tree::ptree& pt,std::string path, char path_separator, X& vt, int& status){
/** 
  \brief Extracts the double value from the property tree

  \param[in] pt Is the property tree from which we want to extract the value
  \param[in] path Is the path of the variable to extract from the property tree. 
                  This is essentially an extended name of the variable to extract.
  \param[in] path_separator The character used to denote different levels of organization, 
  "." is default in versions without explicit path separator
  \param[out] vt The double value of the variable we are extrating from the property tree.
  \param[out] status Is the flag showing the sucess of the extraction: 0 - is no such variable name was found, 
           1 - if we have completed the extraction sucessfully.
*/

  status = 1;
  try{ vt = pt.get<X>(boost::property_tree::ptree::path_type(path, path_separator)); }
  catch(std::exception& e){ status = 0; }
}

template<typename X>
void load(boost::property_tree::ptree& pt,std::string path,vector<X>& vt,int& status){
/** 
  \brief Extracts the vector of double values from the property tree

  \param[in] pt Is the property tree from which we want to extract the value
  \param[in] path Is the path of the variable to extract from the property tree. 
                  This is essentially an extended name of the variable to extract.
  \param[out] vt The vector of double values of the variable we are extrating from the property tree.
  \param[out] status Is the flag showing the sucess of the extraction: 0 - is no such variable name was found, 
           1 - if we have completed the extraction sucessfully.
*/

  X x; int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      load<X>(pt,path+"."+v.first,x,st);
      if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}

template<typename X>
void load(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<X>& vt,int& status){
/** 
  \brief Extracts the vector of double values from the property tree

  \param[in] pt Is the property tree from which we want to extract the value
  \param[in] path Is the path of the variable to extract from the property tree. 
                  This is essentially an extended name of the variable to extract.
  \param[in] path_separator The character used to denote different levels of organization, 
  "." is default in versions without explicit path separator
  \param[out] vt The vector of double values of the variable we are extrating from the property tree.
  \param[out] status Is the flag showing the sucess of the extraction: 0 - is no such variable name was found, 
           1 - if we have completed the extraction sucessfully.
*/

  X x; int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      load<X>(pt,path+std::string(1,path_separator)+v.first, path_separator,x,st);
      if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}






/*
//-------- MATRIX3x3 --------
void set_value(int& defined, MATRIX3x3& value, boost::python::object obj, std::string attrName);
void save(boost::property_tree::ptree& pt,std::string path,MATRIX3x3& vt);
void save(boost::property_tree::ptree& pt,std::string path, char path_separator, MATRIX3x3& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<MATRIX3x3>& vt);
void save(boost::property_tree::ptree& pt,std::string path, char path_separation, vector<MATRIX3x3>& vt);

void load(boost::property_tree::ptree& pt,std::string path,MATRIX3x3& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path, char path_separator, MATRIX3x3& vt, int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<MATRIX3x3>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path, char path_separator, vector<MATRIX3x3>& vt,int& status);
*/


/*
// For convertion between XML datafile and internal C++ representation
void save(boost::property_tree::ptree& pt,std::string path,double& vt);
void save(boost::property_tree::ptree& pt,std::string path,char path_separator, double& vt );

void save(boost::property_tree::ptree& pt,std::string path,vector<double>& vt);
void save(boost::property_tree::ptree& pt,std::string path,char path_separator, vector<double>& vt);

void save(boost::property_tree::ptree& pt,std::string path,int& vt);
void save(boost::property_tree::ptree& pt,std::string path,char path_separator, int& vt );

void save(boost::property_tree::ptree& pt,std::string path,vector<int>& vt);
void save(boost::property_tree::ptree& pt,std::string path,char path_separator, vector<int>& vt );

void save(boost::property_tree::ptree& pt,std::string path,std::string& vt);
void save(boost::property_tree::ptree& pt,std::string path,char path_separator, string& vt );

void save(boost::property_tree::ptree& pt,std::string path,vector<std::string>& vt);
void save(boost::property_tree::ptree& pt,std::string path,char path_separator, vector<std::string>& vt );


void load(boost::property_tree::ptree& pt,std::string path,double& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,char path_separator, double& vt,int& status);

void load(boost::property_tree::ptree& pt,std::string path,vector<double>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,char path_separator,vector<double>& vt,int& status);

void load(boost::property_tree::ptree& pt,std::string path,int& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<int>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,std::string& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<std::string>& vt,int& status);


*/

void save_xml(std::string filename, boost::property_tree::ptree& pt);
void load_xml(std::string filename, boost::property_tree::ptree& pt);



int read_file(std::string filename,int verbose,vector<std::string>& A);

void file2matrix(std::string filename,vector< vector<double> >& m);
void file2matrix(std::string filename,vector< vector<double> >& m,double scl);
void file2matrix(std::string filename,vector< vector<int> >& m);

void show_2D(vector< vector<double> >& in);




}// libio
}// liblibra

#endif // IO_H
