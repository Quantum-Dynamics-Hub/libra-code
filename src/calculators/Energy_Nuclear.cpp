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
  \file Energy_Nuclear.cpp
  \brief The file implement functions for nuclear interaction energy/force calculations
    
*/

#include "Energy_Nuclear.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libcalculators namespace
namespace libcalculators{

double energy_nucl(vector<VECTOR>& R, vector<double>& Zeff){
/**
  \brief Compute nuclear energy of a sub-system

  Simple nuclear-nuclear Coulombic interaction

  \param[in]  R  Coordinates of all atoms in the system
  \param[in] Zeff effective charges of all atoms in the system
*/



  int I,J;
  int sz = R.size();
  double en = 0.0;

  for(int i=0;i<sz;i++){
    for(int j=i+1;j<sz;j++){

      en += Zeff[i] * Zeff[j] / (R[i] - R[j]).length();
      
    }// for j
  }// for i

  return en;

}// energy_nucl


double energy_nucl(vector<VECTOR>& R, vector<double>& Zeff, vector<VECTOR>& G){
/**
  \brief Compute nuclear energy and energy gradients of a sub-system 

  Simple nuclear-nuclear Coulombic interaction. Derivatives w.r.t. all coordinates

  \param[in]  R  Coordinates of all atoms in the system
  \param[in] Zeff Effective charges of all atoms in the system
  \param[out] G  Will contain gradients w.r.t. to each nuclear coordinates
*/


  int I,J;
  int sz = R.size();
  double en = 0.0;
  //if(G.size()!=sz){  G = vector<VECTOR>(sz,VECTOR(0.0,0.0,0.0)); }
  //else{ for(int i=0;i<sz;i++){  G[i] = 0.0; }  }
  if(G.size()>0) { G.clear(); }
  for(int i=0;i<sz;i++){  G.push_back(VECTOR(0.0, 0.0, 0.0)); } 

  for(int i=0;i<sz;i++){
    for(int j=i+1;j<sz;j++){

      double rij = (R[i] - R[j]).length();
      double tmp = Zeff[i] * Zeff[j] / rij;
      en += tmp;

      double modg = -tmp/(rij * rij);
      VECTOR Gij; Gij = (R[i] - R[j]) * modg;

      G[i] += Gij;
      G[j] -= Gij;
      
      
    }// for j
  }// for i

  return en;

}// energy_nucl


}//namespace libcalculators
}// liblibra


