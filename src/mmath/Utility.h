#ifndef UTILITY_H
#define UTILITY_H

#include <boost/python.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <exception>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <set>
#include <exception>
#include <iostream>
#include <stdlib.h>
#include <sstream>


#include "libmmath.h"
#include "Units.h"

using namespace boost::python;



namespace libmmath{

bool hasattr(boost::python::object obj, std::string attrName);

void set_value(int& defined, int& value,         boost::python::object obj, std::string attrName);
void set_value(int& defined, double& value,      boost::python::object obj, std::string attrName);
void set_value(int& defined, std::string& value, boost::python::object obj, std::string attrName);
void set_value(int& defined, VECTOR& value,      boost::python::object obj, std::string attrName);
void set_value(int& defined, MATRIX3x3& value,   boost::python::object obj, std::string attrName);
void set_value(int& defined, MATRIX& value,      boost::python::object obj, std::string attrName);
void set_value(int& defined, QUATERNION& value,  boost::python::object obj, std::string attrName);

void set_list(int& defined, vector<int>& value,boost::python::object obj,std::string attrName);
void set_list(int& defined, vector<double>& value,boost::python::object obj,std::string attrName);
void set_list(int& defined, vector<std::string>& value,boost::python::object obj,std::string attrName);

typedef std::vector<int> IntList;


//------------------------ XML via property_tree ---------------------------
void save(boost::property_tree::ptree& pt,std::string path,double& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<double>& vt);
void save(boost::property_tree::ptree& pt,std::string path,int& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<int>& vt);
void save(boost::property_tree::ptree& pt,std::string path,std::string& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<std::string>& vt);
void save(boost::property_tree::ptree& pt,std::string path,VECTOR& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<VECTOR>& vt);
void save(boost::property_tree::ptree& pt,std::string path,QUATERNION& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<QUATERNION>& vt);
void save(boost::property_tree::ptree& pt,std::string path,MATRIX3x3& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<MATRIX3x3>& vt);
void save(boost::property_tree::ptree& pt,std::string path,MATRIX& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt);

void load(boost::property_tree::ptree& pt,std::string path,double& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<double>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,int& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<int>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,std::string& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<std::string>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,VECTOR& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<VECTOR>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,QUATERNION& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<QUATERNION>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,MATRIX3x3& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<MATRIX3x3>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,MATRIX& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt,int& status);

}// libmmath

#endif // UTILITY_H
