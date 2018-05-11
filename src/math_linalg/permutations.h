/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file permutations.h
  \brief The file describes the functions and data types for dealing with permutations
*/


#ifndef permutations_H
#define permutations_H

#include <complex>
#include <vector>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>


#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>



/// liblibra 
namespace liblibra{

using namespace boost::python;
using namespace std;


/// liblinalg namespace
namespace liblinalg{


vector<int> composite_permutation(vector<int>& perm_t, vector<int>& perm_cum);
void update_permutation(vector<int>& perm_t, vector<int>& perm_cum);
void update_permutation(vector<int>& perm_t, vector<int>* perm_cum);
void check_permutation(vector<int>& perm, int n);



}//namespace liblinalg
}// liblibra

#endif // permutations_H

