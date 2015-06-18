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
