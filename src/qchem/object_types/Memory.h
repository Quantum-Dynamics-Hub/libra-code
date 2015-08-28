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

#ifndef MEMORY_H
#define MEMORY_H

#include "Mathematics.h"


class Memory{

  public:

  // Initialize general-purpose auxiliary arrays
  vector<double*> aux;  int n_aux;
  vector<VECTOR*> auxv; int n_auxv;

  // For INDO - for n fragments
  vector<double> eri; 
  vector<double> V_AB;


  Memory(int n_,int szn_,int nv_,int sznv_){
    n_aux = n_;
    for(int i=0;i<szn_;i++){    // szn_ auxiliary arrays with double
      double* x; x = new double[n_aux];
      aux.push_back(x);
    }

    n_auxv = nv_;
    for(int i=0;i<sznv_;i++){    // sznv_ auxiliary arrays with VECTOR
      VECTOR* x; x = new VECTOR[n_auxv];
      auxv.push_back(x);
    }
  }
  Memory(){
    Memory(20,40, 20,40);
  }

  ~Memory(){
    for(int i=0;i<n_aux;i++) {  delete aux[i];  }  aux.clear();
    for(int i=0;i<n_auxv;i++){  delete auxv[i]; }  auxv.clear();
  }

};

#endif // MEMORY_H