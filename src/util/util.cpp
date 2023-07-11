/***********************************************************
 * Copyright (C) 2013,2016-2022 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/
/**
  \file util.cpp
  \brief The file implements the auxiliary functions for various purposes
    
*/

#include "util.h"


/// liblibra namespace
namespace liblibra{

/// libutil namespace
namespace libutil{


//---------------- Operations on vectors of integers ---------------

void show_vector(vector<int>& A){
/**
  This function simply prints out the vector A
  \param[in] A vector to print out  
*/
  int sz = A.size();  for(int i=0;i<sz;i++){  cout<<A[i]<<"  "; }
}

int is_in_vector(int a, vector<int>& A){
/**
  This function checks whether the integer a is present in the vector of integers A
  The function returns 1 if the element is found, 0 otherwise
  \param[in] a the element we search for in vector of integers
  \param[in] A the vector of integers in which we search for a specific integer
*/

  int res = 0;
  int sz = A.size();
  for(int i=0;i<sz;i++){ if(a==A[i]){ res=1; } }
  return res;
}

int is_in_vector(int a, vector<int>& A, int& pos){
/**
  This function checks whether the integer a is present in the vector of integers A
  The function returns 1 if the element is found, 0 otherwise
  \param[in] a the element we search for in vector of integers
  \param[in] A the vector of integers in which we search for a specific integer
  \param[out] pos  the position of the element a in vector A, if present. It is set to -1 otherwise
  If there are more occurences of a in A, the position of the last one will be returned into pos
*/
  pos = -1;
  int res = 0;
  int sz = A.size();
  for(int i=0;i<sz;i++){ if(a==A[i]){ res=1; pos = i; } }
  return res;
}

boost::python::list is_in_vector2(int a, vector<int>& A){
/**
  This function checks whether the integer a is present in the vector of integers A
  The function returns a list of results, X:
  X[0] -  1 if the element is found, 0 otherwise
  X[1] - the position of the element a in vector A, if present. It is set to -1 otherwise
  If there are more occurences of a in A, the position of the last one will be returned into pos
  \param[in] a the element we search for in vector of integers
  \param[in] A the vector of integers in which we search for a specific integer
*/

  int pos = -1;
  int res = 0;
  int sz = A.size();
  for(int i=0;i<sz;i++){ if(a==A[i]){ res=1; pos = i; } }

  boost::python::list X;
  X.append(res); X.append(pos);

  return X;
}


int is_in_vector(int a, vector<int>& A, vector<int>& indx){
/** 
  This function returns how many times integer-valued element a is found in the vector of integers A
  the argument indx will contain the positions at which a is found in A
  \param[in] a the element we search for in vector of integers
  \param[in] A the vector of integers in which we search for a specific integer
  \param[out] indx the vector containing the positions of vector A in which the element a is found, the size of this vector
  should be equal to the number returned by this function
*/
  int res = 0;
  int sz = A.size();
  if(indx.size()>0){ indx.clear(); }
  for(int i=0;i<sz;i++){ if(a==A[i]){ res++; indx.push_back(i); } }
  return res;
}


int is_repeating(vector<int>& A,int& reap){
/**
  This function finds out if there are repeating elements in vector A
  \param[in] A the original input vector of integers
  \param[out] reap the first integer found repeating
*/
  int res = 0;
  int sz = A.size();

  for(int i=0;i<sz-1;i++){
    vector<int> tmp = std::vector<int>(A.begin()+i+1,A.end());
    if(is_in_vector(A[i],tmp)){ res = 1; reap = A[i]; break; }
  }

  return res;
}

boost::python::list is_repeating(vector<int>& A){
/**
  This function finds out if there are repeating elements in vector A
  Returns the results in the list X:
  X[0] - 1-if there are repeating elements, 0 - if there are no
  X[1] - the first value that is present in the input vetor A more than 1 time
  \param[in] A the original input vector of integers
*/
  int res = 0;
  int sz = A.size();
  int reap;

  for(int i=0;i<sz-1;i++){
    vector<int> tmp = std::vector<int>(A.begin()+i+1,A.end());
    if(is_in_vector(A[i],tmp)){ res = 1; reap = A[i]; break; }
  }

  boost::python::list X;
  X.append(res); X.append(reap);

  return X;
}



int delta(vector<int>& A,vector<int>& B,int& a,int& b){
/** 
  This function searches for pair of indices in A and B
  which are different, they are then assigned to a and b correspondingly
  All other indexes in A and B should be equivalent (not necessarily in the same order)

  The function returns 1 if only 1 pair of indices is different between A and B,
  return 0 otherwise

  That is res = 1 means A and B are coupled, while res = 0 - are not

  Algorithm: first construct the list of overlapping orbitals
  second - go through this list and determine the populations of
  each orbital in this list in each of the configurations

  THIS VERSION REGARDS SIGN OF THE ORBITAL INDEX

  \param[in] A first vector to compare
  \param[in] B second vector to compare
  \param[out] a the integer from A which is not present in B
  \param[out] b the integer from B which is not present in A
*/
  int i;
  int res = 1;
  int sz = A.size();
  int nA = 0; // number of elements in A which do not exist in B
  int nB = 0; // number of elements in B which do not exist in A
  std::vector<int> aA, aB, aC; // modules of A and B and overlap C

  for(i=0;i<sz;i++){ // the size of A and B is assumed to be the same
    int mA = A[i];
    int mB = B[i];
    aA.push_back(mA);
    aB.push_back(mB);
    if(!is_in_vector(mA,aC)){ aC.push_back(mA);}
    if(!is_in_vector(mB,aC)){ aC.push_back(mB);}
  }// for i

  sz = aC.size();

  int nexc = 0; // Number of excitations between 2 states
  for(i=0;i<sz;i++){
    int n_in_a,n_in_b;
    vector<int> tmpa,tmpb;
    n_in_a = is_in_vector(aC[i],aA,tmpa);
    n_in_b = is_in_vector(aC[i],aB,tmpb);
    int d = n_in_a - n_in_b;
    if (d>0){ nexc += d;   }
    if(d==1){ a = aC[i]; }
    if(d==-1){ b = aC[i]; }
  }// for i

  if(nexc==1){  res = 1;
    // The orbitals should already be known, because only one pair is different
  }
  else{ res = 0; }


  return res;
}

boost::python::list delta(vector<int>& A,vector<int>& B){

  int a, b, res;
  res = delta(A,B,a,b);

  boost::python::list X;
  X.append(res); X.append(a); X.append(b);

  return X;
}


//-------------------- Operations on strings ----------------------------

void split_line(std::string line, vector<std::string>& arr){
/**
  This function splits a long string into a vector of smaller sub-strings (delimited by the newline)
  \param[in] line The original long string
  \param[out] arr The generated vector of sub-strings
*/

  stringstream ss(line,stringstream::in|stringstream::out);
  std::string s;
  while(ss>>s){ arr.push_back(s); }

}

void split_line(std::string line,vector<std::string>& arr,char delim){
/**
  This function splits a long string into a vector of smaller sub-strings using a
  specified delimiter
  \param[in] line The original long string
  \param[out] arr The generated vector of sub-strings
  \param[in] delim The delimeter to use
*/


  std::istringstream f(line);
  std::string s;    
  while(std::getline(f, s, delim)){ arr.push_back(s); }

}

std::string int2str(int inp){
/**
  This is a function that converts integer into a corresponding string
  \param[in] inp the input value
*/

  stringstream ss(stringstream::in | stringstream::out);
  std::string out;
  (ss << inp);  ss >> out;
  return out;
}


int find_section(vector<std::string>& A,std::string marker_beg,std::string marker_end,int min_line,int max_line,int& beg,int& end){ 

  beg = end = -1;
  size_t found;
  int status = 0;

  for(int i=min_line;i<max_line;i++){
    if(beg==-1){
      found = A[i].find(marker_beg);
      if(found!=string::npos){  beg = i;}
    }
    if(end==-1){
      found = A[i].find(marker_end);
      if(found!=string::npos){  end = i;}
    }
    if((beg!=-1) && (end!=-1)){ 
      if(beg<end){status = 1;}
      break;
    }
  }// for i

  return status;
}

std::string extract_s(std::string line, std::string marker){
// line - is the line from which we want to extract some data
// marker - is the keyword before that data, format is:  marker="data"
  size_t pos,pos1,pos2;
  string res;

  pos = line.find(marker);
  if(pos!=string::npos){  // Marker is found
    pos1 = line.find("\"",pos+1);
    pos2 = line.find("\"",pos1+1);
    if(pos1!=string::npos && pos2!=string::npos){
      res = line.substr(pos1+1,pos2-pos1-1); // extract the value
    }
  }
  return res;
}



//-------------------- Operations on arrays ----------------------------


void extract_2D(vector< vector<double> >& in, vector< vector<double> >& out, int minx,int maxx, int miny, int maxy ){
/*****************************************************************
  This function extracts 2D sub-array out from 2D array in
  The region is degined by indices: [minx,maxx]x[miny,maxy]
  Axes are: in[x][y], out[x][y]
******************************************************************/
  if(out.size()>0){ out.clear(); }
  for(int x=minx;x<=maxx;x++){
    vector<double> line = std::vector<double>(in[x].begin()+miny,in[x].begin()+maxy+1);
    out.push_back(line);
  }
}

void extract_2D(vector< vector<double> >& in, vector< vector<double> >& out, vector<int>& templ,int shift){
/*****************************************************************
  This function extracts 2D sub-array out from 2D array in
  The region is degined by indices: templ x templ
  temp - is a template for extraction
******************************************************************/
  if(out.size()>0){ out.clear(); }
  int sz = templ.size();
  for(int i=0;i<sz;i++){
    vector<double> line = std::vector<double>(sz,0.0);
    for(int j=0;j<sz;j++){  line[j] = in[templ[i]+shift][templ[j]+shift];   }
    out.push_back(line);
  }
}

void extract_1D(vector<double>& in, vector<double>& out, vector<int>& templ,int shift){
/*****************************************************************
  This function extracts 1D sub-array out from 1D array in
  The region is degined by indices: templ x templ
  temp - is a template for extraction
******************************************************************/
  if(out.size()>0){ out.clear(); }
  int sz = templ.size();
  out = std::vector<double>(sz,0.0);
  for(int i=0;i<sz;i++){  out[i] = in[templ[i]+shift];   }
}
 



int is_equal(vector<int>& vec1, vector<int>& vec2){
/*****************************************************************
  This function checks if two vectors of ints are equal

  Returns: 
  0 - if not equal = either theis sizes are different or at least 1 element is different
  1 - if equal - the sizes are the same, and all elements are equal to each other
******************************************************************/

  int res = 1; 
  int sz1 = vec1.size();
  int sz2 = vec2.size();

  if(sz1!=sz2){ res = 0; }
  else{
    for(int i=0; i<sz1; i++){
      if(vec1[i] != vec2[i]){ res = 0; break; }
    }
  }
  return res;
}

int is_included(vector<int>& vec1, vector<vector<int> >& vec){
/*****************************************************************
  This function checks if a vector of ints `vec1` is already present in a list
  of such vectors `vec`

  Returns: 
  0 - if not present
  1 - if present
******************************************************************/

  int res = 0;
  for(int i=0; i<vec.size(); i++){
    if(is_equal(vec1, vec[i])){ res = 1; break; }
  }
  return res;
}


int is_included(vector<int>& vec1, vector<vector<int> >& vec, int start, int num_of_elts){
/*****************************************************************
  This function checks if a vector of ints `vec1` is already present 
  in the interval [start, start+num_of_elts) in the list of such vectors `vec`

  Returns: 
  0 - if not present
  1 - if present
******************************************************************/

  int res = 0;
  for(int i=start; i<(start + num_of_elts); i++){
    if(is_equal(vec1, vec[i])){ res = 1; break; }
  }
  return res;
}



int is_present(vector< vector<int> >& vec, int i, int start, int end){
/**
  Checks whether the vector of the row i is already present among the rows [start, end)

  return 1 if the row i is already present among those rows
*/

  int d = vec[0].size();
 
  int res = 0;
  for(int j=start; j<end; j++){

    int are_same = 1;
    for(int k=0;k<d;k++){
      if( vec[i][k] != vec[j][k] ){ are_same = 0; break; }
    }

    if(are_same == 1){ res = 1; break; }
  }

  return res;

}


int sum_row(int row, vector<vector<int> >& vec){
/**
  Computes the sum of the elements in a given row
*/

  int res = 0;
  int ncols = vec[0].size();
 
  for(int i=0;i<ncols;i++){ res += vec[row][i]; }

  return res;

}



vector<int> allocate_1D(int sz1){

  vector<int> res(sz1, 0);

  return res;

}

vector< vector<int> > allocate_2D(int sz1, int sz2){

  vector< vector<int> > res(sz1, vector<int>(sz2, 0) );

  return res;
}

vector< vector< vector<int> > > allocate_3D(int sz1, int sz2, int sz3){

  vector< vector< vector<int> > > res(sz1, vector< vector<int> >(sz2,  vector<int>(sz3, 0) ) );

  return res;
}






void check_input(boost::python::dict params, boost::python::dict default_params, boost::python::list critical_params){
/*****************************************************************

    This function checks the input Python dictionary and verifies if it contains
    all the required keys. Exits if there are no. It also checks other variables 
    and sets the values of these keys to the default values, if they weren't defined

    Args:
        params ( Python dictionary ): the dictionary (usually containing parameters) to be checked
            This dictionary is passed by reference and can be changed inside this function.

        default_params ( Python dictionary ): the dictionary with the default key-value pairs to 
            be inserted into the `params` variable, if those pairs are not defined in it.

        critical_params ( list of strings ): the names of the keys that must be present in the
            `params` variable. The function exits and prints out an error, if any of the keys in 
            the `critical_params` variable are not found in `params`.

    Returns:
        None: 

            This function checks whether the required keys are defined in the input.
            If not - print out a worning and set the parameters to the default values provided.

            But the critical parameters should be present - otherwise will exit

            Revises the parameters dictionary - but only according to the subset of 
            default parameters - other ones are not affects

    Example:

        Assume there is a function called `my_function`, which prints a greeting message.
        So, it needs to know the name. There is too much of an ambiguity to set up a default name.
        Therefore, the "Name" key is required (critical) for the function not to crash. 
        Also, the function prints out the age of the person. If that info is provided, no problem.
        But if the user forgets about it, one can guess something more-or-less reasonable. So,
        in this example, the "Age" key is set up to a default value before we call the `my_function`.
          
        >>> def my_function(data, params):
        >>>     print "Hello ", params["Name"]
        >>>     print "You age is ", params["Age"]


        The following snippet 

        >>> params = {"Name":"Alexey", "Age":33 }
        >>> default_params = { "Department": "Chemistry", "Institution":"UB", "Age":1 }
        >>> critical_parames = ["Name"]
        >>> check_input(params, default_params, critical_params)
        >>> my_function(data, params)

        will produce:
        
        >>> Hello Alexey
        >>> Your age is 33


        The following snippet 

        >>> params = {"Name":"Alexey" }
        >>> default_params = { "Department": "Chemistry", "Institution":"UB", "Age":1 }
        >>> critical_parames = ["Name"]
        >>> check_input(params, default_params, critical_params)
        >>> my_function(data, params)

        will produce:
        
        >>> Hello Alexey
        >>> Your age is 1


        The following snippet 

        >>> params = { }
        >>> default_params = { "Department": "Chemistry", "Institution":"UB", "Age":1 }
        >>> critical_parames = ["Name"]
        >>> check_input(params, default_params, critical_params)
        >>> my_function(data, params)

        will crash:

        >>> ERROR: The critical parameter Name must be defined!
        >>> Exiting now...


******************************************************************/

  int i,j, is_found;
  std::string key;


  for(i=0;i<len(critical_params);i++){
      key = boost::python::extract<std::string>(critical_params[i]);

      if(!params.has_key(key)){
          std::cout<<"ERROR: The critical parameter "<< key <<" must be defined!\nExiting now...\n";
          exit(0);
      }
  }


  for(i=0;i<len(default_params.keys());i++){
      key = boost::python::extract<std::string>(default_params.keys()[i]);

      if(!params.has_key(key)){
          std::cout<<"WARNING: Parameter "<< key<< " is not defined! in the input parameters"
                   <<"Use the default value \n"; //<< default_params.values()[i]<<"\n";

          params.setdefault(key, default_params.values()[i]);
      }
  }
}



}// libutil

}// liblibra



