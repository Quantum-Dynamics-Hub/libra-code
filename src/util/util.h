/***********************************************************
 * Copyright (C) 2013, 2016-2017 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/
/**
  \file util.h
  \brief The file describes some auxiliary functions for various general purposes
    
*/

#ifndef UTIL_H
#define UTIL_H

#if defined(USING_PCH)
#include "../pch.h"
#else

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <boost/python.hpp>

#endif 

/// liblibra namespace
namespace liblibra{

using namespace std;
using namespace boost::python;


/// liutil namespace
namespace libutil{


// Operations on the the vectors of integers
void show_vector(vector<int>& A);

int is_in_vector(int a, vector<int>& A);
int is_in_vector(int a, vector<int>& A, int& pos);
int is_in_vector(int a, vector<int>& A, vector<int>& indx);
boost::python::list is_in_vector2(int a, vector<int>& A);


int is_repeating(vector<int>& A,int& reap);
boost::python::list is_repeating(vector<int>& A);

int delta(vector<int>& A,vector<int>& B,int& a,int& b);
boost::python::list delta(vector<int>& A,vector<int>& B);

// Operations on strings ("lines")
void split_line(std::string line, vector<std::string>& arr);
void split_line(std::string line,vector<std::string>& arr,char delim);
std::string int2str(int inp);
int find_section(vector<std::string>& A,std::string marker_beg,std::string marker_end,int min_line,int max_line,int& beg,int& end);
std::string extract_s(std::string line, std::string marker);    

// Operations on arrays - reformings, etc.
void extract_1D(vector<double>& in, vector<double>& out, vector<int>& templ,int shift);
void extract_2D(vector< vector<double> >& in, vector< vector<double> >& out, int minx,int maxx, int miny, int maxy );
void extract_2D(vector< vector<double> >& in, vector< vector<double> >& out, vector<int>& templ,int shift);

// Checking the presence of vectors in lists of vectors
int is_equal(vector<int>& vec1, vector<int>& vec2);
int is_included(vector<int>& vec1, vector<vector<int> >& vec);
int is_included(vector<int>& vec1, vector<vector<int> >& vec, int start, int num_of_elts);
int is_present(vector< vector<int> >& vec, int i, int start, int end);

// Properties of vectors
int sum_row(int row, vector<vector<int> >& vec);

// Allocating storage
vector<int> allocate_1D(int sz1);
vector< vector<int> > allocate_2D(int sz1, int sz2);
vector< vector< vector<int> > > allocate_3D(int sz1, int sz2, int sz3);



void check_input(boost::python::dict params, boost::python::dict default_params, boost::python::list critical_params);


}// libutil

}// liblibra

#endif // UTIL_H

