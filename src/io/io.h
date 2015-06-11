#ifndef IO_H
#define IO_H


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

using namespace std;
using namespace boost::python;
using boost::property_tree::ptree;


namespace libio{


bool hasattr(boost::python::object obj, std::string attrName);

/// For extracting values from Python object to C++ representation
void set_value(int& is_defined, int& value,         boost::python::object obj, std::string attrName);
void set_value(int& is_defined, double& value,      boost::python::object obj, std::string attrName);
void set_value(int& is_defined, std::string& value, boost::python::object obj, std::string attrName);
void set_list(int& is_defined, vector<int>& value,        boost::python::object obj,std::string attrName);
void set_list(int& is_defined, vector<double>& value,     boost::python::object obj,std::string attrName);
void set_list(int& is_defined, vector<std::string>& value,boost::python::object obj,std::string attrName);




//------------------------ XML via property_tree ---------------------------
/// For convertion between XML datafile and internal C++ representation
void save(boost::property_tree::ptree& pt,std::string path,double& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<double>& vt);
void save(boost::property_tree::ptree& pt,std::string path,int& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<int>& vt);
void save(boost::property_tree::ptree& pt,std::string path,std::string& vt);
void save(boost::property_tree::ptree& pt,std::string path,vector<std::string>& vt);

void load(boost::property_tree::ptree& pt,std::string path,double& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<double>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,int& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<int>& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,std::string& vt,int& status);
void load(boost::property_tree::ptree& pt,std::string path,vector<std::string>& vt,int& status);


void save_xml(std::string filename, boost::property_tree::ptree& pt);
void load_xml(std::string filename, boost::property_tree::ptree& pt);


}// libio

#endif // IO_H
