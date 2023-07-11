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
  \file io.cpp
  \brief The file implements auxiliary (but basic) functions for Python input/output operations
    
*/

#include "io.h"


/// liblibra 
namespace liblibra{

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


//  boost::property_tree::xml_writer_settings<char> settings(' ', 4);
//  write_xml(filename, pt, std::locale(), settings);
  write_xml(filename, pt, std::locale(), xml_writer_settings('\t', 1));

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




int read_file(std::string filename,int verbose,vector<std::string>& A){
/**********************************************************************
  This function reads file <filename> and stores it as a vector of strings
  each string is a line of the file
**********************************************************************/

  // Prepare A
  if(A.size()>0){ A.clear(); }

  // Read the file
  if(verbose==1){ cout<<"Reading file"<<filename<<endl; }
  ifstream f;
  f.open(filename.c_str(),ios::in);
  if(f.is_open()){
    // Estimate the number of lines
    int nlines = 0;
    while(!f.eof()){ std::string line; getline(f,line); nlines++; }
    A.reserve(int(nlines*1.25));
    // Reset file position to beginning
    f.clear();
    f.seekg (0, ios::beg);

    // Actually put the file content in memory
    while(!f.eof()){ std::string line; getline(f,line); A.push_back(line);  }
  }else{ cout<<"Error: Can not open file "<<filename<<endl; }
  f.close();

  return A.size();
}


void file2matrix(std::string filename,vector< vector<double> >& M){
/*****************************************************************
  This function reads the content of the tabular (2D) file into matrix M
*****************************************************************/
  if(M.size()>0){ M.clear(); }
  ifstream in;
  in.open(filename.c_str(),ios::in);
  if(in.is_open()){
    std::string s;
    while(!in.eof()){
      vector<double> line; // line of the file
      getline(in,s);
      stringstream ss(s,stringstream::in|stringstream::out);
      while(ss>>s){  line.push_back(atof(s.c_str()));   }
      if(line.size()>0) { M.push_back(line); }
    }
  }else{ cout<<"Error: Can not open file "<<filename<<". Check if this file exists\n"; exit(0);}
  in.close();
}

void file2matrix(std::string filename,vector< vector<double> >& M,double scl){
/*****************************************************************
  This function reads the contend of the tabular (2D) file into matrix M and scales
  the read data by factor scl
*****************************************************************/
  if(M.size()>0){ M.clear(); }
  ifstream in;
  in.open(filename.c_str(),ios::in);
  if(in.is_open()){
    std::string s;
    while(!in.eof()){
      vector<double> line; // line of the file
      getline(in,s);
      stringstream ss(s,stringstream::in|stringstream::out);
      while(ss>>s){  line.push_back(scl*atof(s.c_str()));   }
      if(line.size()>0) { M.push_back(line); }
    }
  }else{ cout<<"Error: Can not open file "<<filename<<". Check if this file exists\n"; exit(0);}
  in.close();
}


void file2matrix(std::string filename,vector< vector<int> >& M){
/*****************************************************************
  This function reads the content of the tabular (2D) file into matrix M
  Version overloaded for int
*****************************************************************/
  if(M.size()>0){ M.clear(); }
  ifstream in;
  in.open(filename.c_str(),ios::in);
  if(in.is_open()){
    std::string s;
    while(!in.eof()){
      vector<int> line; // line of the file
      getline(in,s);
      stringstream ss(s,stringstream::in|stringstream::out);
      while(ss>>s){  line.push_back(atoi(s.c_str()));   }
      if(line.size()>0) { M.push_back(line);     }
    }
  }else{ cout<<"Error: Can not open file "<<filename<<". Check if this file exists\n"; exit(0);}
  in.close();
}

void show_2D(vector< vector<double> >& in){
/******************************************************************
  This function prints out the matrix in a tabular form
******************************************************************/
  for(int i=0;i<in.size();i++){
    for(int j=0;j<in[i].size();j++){
      cout<<"in["<<i<<"]["<<j<<"]="<<in[i][j]<<" ";
    }
    cout<<endl;
  }
}





}// libio
}// liblibra

