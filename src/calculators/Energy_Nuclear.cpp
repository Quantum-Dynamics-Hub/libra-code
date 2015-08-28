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

#include "Energy_Nuclear.h"

namespace libcalculators{

double energy_nucl(vector<VECTOR>& R, vector<double>& Zeff){
/// Compute nuclear energy of a sub-system

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

}//namespace libcalculators
