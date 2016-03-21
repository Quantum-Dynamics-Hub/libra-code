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
/**
  \file io.cpp
  \brief The file implements auxiliary (but basic) functions for Python input/output operations
    
*/

#include "io.h"

/// libio namespace
namespace libio{


bool hasattr(boost::python::object obj, std::string attrName) {
/** 
  \brief Check if the Python object has a particularly-named datamember 

  \param[in] obj An object of Python class (note, this is not just any Python object). We will try to extract date from
                 that object
  \param[in] attrName The string containing the name of the attribute to be looked for in the object of Python class  
*/
     return PyObject_HasAttrString(obj.ptr(), (char*)attrName.c_str());
} 



void set_value(int& is_defined, int& value, boost::python::object obj, std::string attrName){
/** 
  \brief For extracting integer values from Python object to C++ representation 

  \param[out] is_defined In the end, it is set to 0, if no data member with requested name has been found in the Python object
              It is set to 1, if such data member is found. In this case, we also extract the value of the found member.
  \param[out] value This is the value to be extracted
  \param[in] obj The Python object (instance of a Python class), from which we extract the value
  \param[in] attrName The string containing the name of the attribute to be looked for in the object of Python class  
*/

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);        
  if(has_attr){
      value = extract<int>(obj.attr(attrName.c_str()));          
      is_defined = 1; 
  }

}

void set_value(int& is_defined, double& value, boost::python::object obj, std::string attrName){
/** 
  \brief For extracting double values from Python object to C++ representation 

  \param[out] is_defined In the end, it is set to 0, if no data member with requested name has been found in the Python object
              It is set to 1, if such data member is found. In this case, we also extract the value of the found member.
  \param[out] value This is the value to be extracted
  \param[in] obj The Python object (instance of a Python class), from which we extract the value
  \param[in] attrName The string containing the name of the attribute to be looked for in the object of Python class  
*/

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<double>(obj.attr(attrName.c_str()));
      is_defined = 1;
  }

}

void set_value(int& is_defined, std::string& value,boost::python::object obj, std::string attrName){
/** 
  \brief For extracting string values from Python object to C++ representation 

  \param[out] is_defined In the end, it is set to 0, if no data member with requested name has been found in the Python object
              It is set to 1, if such data member is found. In this case, we also extract the value of the found member.
  \param[out] value This is the value to be extracted
  \param[in] obj The Python object (instance of a Python class), from which we extract the value
  \param[in] attrName The string containing the name of the attribute to be looked for in the object of Python class  
*/

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<std::string>(obj.attr(attrName.c_str()));
      is_defined = 1;
  }
}


void set_list(int& is_defined, vector<int>& value,boost::python::object obj,std::string attrName){
/** 
  \brief For extracting vector of integers values from Python object to C++ representation 

  \param[out] is_defined In the end, it is set to 0, if no data member with requested name has been found in the Python object
              It is set to 1, if such data member is found. In this case, we also extract the value of the found member.
  \param[out] value This is the value to be extracted
  \param[in] obj The Python object (instance of a Python class), from which we extract the value
  \param[in] attrName The string containing the name of the attribute to be looked for in the object of Python class  

  The vector of integers is represented in the Python object by a list of integers. This list is a class member itself.
*/


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
/** 
  \brief For extracting vector of double values from Python object to C++ representation 

  \param[out] is_defined In the end, it is set to 0, if no data member with requested name has been found in the Python object
              It is set to 1, if such data member is found. In this case, we also extract the value of the found member.
  \param[out] value This is the value to be extracted
  \param[in] obj The Python object (instance of a Python class), from which we extract the value
  \param[in] attrName The string containing the name of the attribute to be looked for in the object of Python class  

  The vector of integers is represented in the Python object by a list of doubles. This list is a class member itself.
*/


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
/** 
  \brief For extracting vector of string values from Python object to C++ representation 

  \param[out] is_defined In the end, it is set to 0, if no data member with requested name has been found in the Python object
              It is set to 1, if such data member is found. In this case, we also extract the value of the found member.
  \param[out] value This is the value to be extracted
  \param[in] obj The Python object (instance of a Python class), from which we extract the value
  \param[in] attrName The string containing the name of the attribute to be looked for in the object of Python class  

  The vector of integers is represented in the Python object by a list of strings. This list is a class member itself.
*/


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












//----------------------------------------------------

void save_xml(std::string filename, boost::property_tree::ptree& pt){
/** 
  \brief Save the property tree object as an XML file

  \param[in] filename The name of the file to which the property tree object will be printed out
  \param[in] pt The property tree containing the internal representation of the XML tree (data)

*/


  boost::property_tree::xml_writer_settings<char> settings(' ', 4);
  write_xml(filename, pt, std::locale(), settings);

}

void load_xml(std::string filename, boost::property_tree::ptree& pt){
/** 
  \brief Load the property tree object from an XML file

  \param[in] filename The name of the file from which the property tree object will be loaded
  \param[in,out] pt The property tree to which the loaded data will be added.

*/


  read_xml(filename, pt);

}

/*
std::string load_xml(std::string filename, boost::property_tree::ptree& pt){

  read_xml(filename, pt);

}
*/



}// libio
