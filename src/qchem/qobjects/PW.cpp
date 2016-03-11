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
  \file PW.cpp
  \brief The file implements basic operations on/with plane-wave objects    
*/

#include "PW.h"

/// libqchem namespace
namespace libqchem{

/// libqobjects namespace
namespace libqobjects{


//------------------ Members of the PW class -------------------------
PW PW::operator-(){
/**
  Negation operator: Returns a new wavefunction with all coefficients c_new[i] = -c_old[i]
*/

  PW res(*this); 
  for(int i=0;i<npw;i++){ res.coeff[i] = -coeff[i]; }
  return res;
}

PW PW::operator+(const PW& m){
/**
  Addition operator: Returns a new wavefunction with all coefficients c_new[i] = c_old[i] + m[i],
  but only if the k-points are the same
*/

  if(!(m.kx==this->kx && m.ky==this->ky && m.kz==this->kz)){  
    cout<<"Error in PW::operator+: can add orbitals only if the k-points are the same!\n";
    exit(0);
  }

  PW res(*this);
  for(int i=0;i<npw;i++){ res.coeff[i] = coeff[i] + m.coeff[i]; }
  return res;

}

PW PW::operator-(const PW& m){
/**
  Subtraction operator: Returns a new wavefunction with all coefficients c_new[i] = c_old[i] - m[i],
  but only if the k-points are the same
*/

  if(!(m.kx==this->kx && m.ky==this->ky && m.kz==this->kz)){  
    cout<<"Error in PW::operator-: can subtract orbitals only if the k-points are the same!\n";
    exit(0);
  }

  PW res(*this);
  for(int i=0;i<npw;i++){ res.coeff[i] = coeff[i] - m.coeff[i]; }
  return res;

}

void PW::operator+=(const PW& m){
/**
  Increment operator: Modifies existing coefficient by: c[i] = c[i] + m[i],
  but only if the k-points are the same
*/

  if(!(m.kx==this->kx && m.ky==this->ky && m.kz==this->kz)){  
    cout<<"Error in PW::operator+=: can increment orbitals only if the k-points are the same!\n";
    exit(0);
  }

  for(int i=0;i<npw;i++){ coeff[i] += m.coeff[i]; }

}

void PW::operator-=(const PW& m){
/**
  Decrement operator: Modifies existing coefficient by: c[i] = c[i] - m[i],
  but only if the k-points are the same
*/

  if(!(m.kx==this->kx && m.ky==this->ky && m.kz==this->kz)){  
    cout<<"Error in PW::operator-=: can increment orbitals only if the k-points are the same!\n";
    exit(0);
  }

  for(int i=0;i<npw;i++){ coeff[i] -= m.coeff[i]; }

}


PW PW::operator/(double num){
/**
  Returns a new PW object with the coefficients divided by a real number: c_new[i] = c_old[i] / num
*/

  PW res(*this);
  for(int i=0;i<npw;i++){ res.coeff[i] /= num;  }
  return res;
}

PW PW::operator/(complex<double> num){
/**
  Returns a new PW object with the coefficients divided by a complex number: c_new[i] = c_old[i] / num
*/

  PW res(*this);
  for(int i=0;i<npw;i++){ res.coeff[i] /= num; }
  return res;
}


PW PW::conj(){
/**
  Returns a PW object which is conjugate to original one
*/

  PW res(*this);
  for(int i=0;i<npw;i++){ res.coeff[i] = std::conj(coeff[i]); }
  return res;
}

void PW::normalize(){
/**
  Normalizes the PW object. Function modifies the coefficients of the original object
*/

  double norm = 0.0;
  for(int i=0;i<npw;i++){ norm += (std::conj(coeff[i]) * coeff[i] ).real();  }
  norm = sqrt(1.0/norm);
  for(i=0;i<npw;i++){ coeff[i] *= norm; }
}

void PW::complete(){
/** 
  Complete the PW object by adding the complex conjugate part
  This is the storage saving trick, due to c(-G) = conj(G), which is usually applied only to
  gamma-point (k_point==0). We also assume that the very first coefficient coeff[0] corresponds
  to G-point = 0,0,0 This is true for QE wavefunctions
  
  The function changes the original object
*/
  if(kx==0 && ky==0 && kz==0){

    // Increase storage
    coeff.resize(2*npw-1);

    // Add new coefficients and compute norm of the 'completed' wavefunction
    double norm = (std::conj(coeff[0]) * coeff[0] ).real();
    for(int i=1;i<npw;i++){
      coeff[npw+i] = std::conj(coeff[i]); 
      norm += 2.0*(std::conj(coeff[i]) * coeff[i] ).real();
    }

    // Update the number of the planewaves in new (completed) wfc
    npw =  2*npw - 1;

    // Finally, normalize the completed wfc
    norm = sqrt(1.0/norm);
    for(i=0;i<npw;i++){ coeff[i] *= norm; }

  }
  else{
    cout<<"Error in PW::complete() - can't complete non-Gamma-point orbitals\n";
    exit(0);
  }

    
}


template<class T>
PW multiply(T& f, const PW& m1){
  PW res(m1);
  for(int i=0;i<res.npw;i++){ res.coeff[i] *= f; }
  return res;
}

PW operator*(const double& f, const PW& m1){ return multiply<const double>(f,m1); }
PW operator*(const PW& m1, const double& f){ return multiply<const double>(f,m1); }
PW operator*(const float& f, const PW& m1){ return multiply<const float>(f,m1); }
PW operator*(const PW& m1, const float& f){ return multiply<const float>(f,m1); }
PW operator*(const complex<double>& f, const PW& m1){ return multiply<const complex<double> >(f,m1); }
PW operator*(const PW& m1, const complex<double>& f){ return multiply<const complex<double> >(f,m1); }
PW operator*(const complex<float>& f, const PW& m1){ return multiply<const complex<float> >(f,m1); }
PW operator*(const PW& m1, const complex<float>& f){ return multiply<const complex<float> >(f,m1); }


complex<double> operator*(const PW& m1, const PW& m2){
  complex<double> res(0.0,0.0);
  if(m1.npw!=m2.npw){ cout<<"Error: Can not multiply MOs with different basis sizes\n"; exit(0); }
  
  else{   for(int i=0;i<m1.npw;i++){ res += m1.coeff[i] * m2.coeff[i];}  }
  return res;
}


}// namespace libqobjects
}// namespace libqchem

