/*********************************************************************************
* Copyright (C) 2018-2020 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

/**
  \file permutations.cpp
  \brief The file describes the functions and data types for dealing with permutations
*/


#include "permutations.h"

/// liblibra namespace
namespace liblibra{

/// liblinalg namespace
namespace liblinalg{


vector<int> id_permutation(int sz){
  vector<int> res(sz,0);
  for(int i=0;i<sz;i++){ res[i] = i; }
  return res;
}

vector<int> inverse_permutation(vector<int>& perm){
/**
  Compute the inverse of the permutation <perm>

  Forward permutation:  i -> perm[i]
  
  this  inv_perm[ perm[i] ] = i
*/

  int sz = perm.size();
  vector<int> res(sz,0);

  for(int i=0;i<sz;i++){    res[ perm[i] ] = i;   }

  return res;
}

vector<int> composite_permutation(vector<int>& perm_t, vector<int>& perm_cum){
    /**
    \param[in] perm_cum - Cumulative permutation
    \param[in] perm_t - permutation at a given point

    This function computes a composition of 
    the two permutations: perm_t (x) perm_cum

    E.g. if perm_cum = [1,0], meaning perm_cum(0) = 1, perm_cum(1) = 0
    and perm_t = [0,1], meaning perm_t(0) = 0, perm_t(1) = 1
    Then:
    perm_t (x) perm_cum = [1, 0]
    (perm_t * perm_cum)(0) = perm_t( perm_cum(0) ) = perm_t(1) = 1
    (perm_t * perm_cum)(1) = perm_t( perm_cum(1) ) = perm_t(0) = 0
    */

    int n = perm_cum.size();
    vector<int> tmp(n, 0);

    for(int i=0;i<n;i++){  tmp[i] = perm_t[ perm_cum[i] ];    }

    return tmp;

}

void update_permutation(vector<int>& perm_t, vector<int>* perm_cum){

    /**
    \param[in, out] perm_cum - Cumulative permutation
    \param[in] perm_t - permutation at a given point

    This function computes a composition of 
    the two permutations: perm_t (x) perm_cum

    E.g. if perm_cum = [1,0], meaning perm_cum(0) = 1, perm_cum(1) = 0
    and perm_t = [0,1], meaning perm_t(0) = 0, perm_t(1) = 1
    Then:
    perm_t (x) perm_cum = [1, 0]
    (perm_t * perm_cum)(0) = perm_t( perm_cum(0) ) = perm_t(1) = 1
    (perm_t * perm_cum)(1) = perm_t( perm_cum(1) ) = perm_t(0) = 0
    */

    int i;
    int n = perm_cum[0].size();
    vector<int> tmp(n, 0);

    for(i=0;i<n;i++){  tmp[i] = perm_t[ perm_cum[0][i] ];    }
 
    // Copy the temporary vector into the input one:
    for(i=0;i<n;i++){  perm_cum[0][i] = tmp[i];   }


}

void update_permutation(vector<int>& perm_t, vector<int>& perm_cum){

  update_permutation(perm_t, &perm_cum);

}

void check_permutation(vector<int>& perm, int n){

    int sz = perm.size();

    if(sz!=n){
      cout<<"ERROR in void check_permutation(vector<int>& perm, int n)\n";
      cout<<"The number of items in the input permutation ("<<sz<<") is not equal to the expected size of the permutation ("<<n<<")\n";
      cout<<"Exiting...\n";
      exit(0);
    }

    int min_indx = perm[0];
    int max_indx = perm[0];
    for(int i=0;i<sz; i++){ 
      if(perm[i] < min_indx){  min_indx = perm[i]; }
      if(perm[i] > max_indx){  max_indx = perm[i]; }
    }

    if(min_indx < 0){
      cout<<"ERROR in void check_permutation(vector<int>& perm, int n)\n";
      cout<<"The minimal mapping index of the permutation ("<<min_indx<<") is smaller than 0\n";
      cout<<"Exiting...\n";
      exit(0);
    }

    if(max_indx > n-1){
      cout<<"ERROR in void check_permutation(vector<int>& perm, int n)\n";
      cout<<"The maximal mapping index of the permutation ("<<max_indx<<") is larger than the \
             maximal allowed index of permutation ("<<n-1<<")\n";
      cout<<"Exiting...\n";
      exit(0);
    }

}





}//namespace liblinalg
}// liblibra
