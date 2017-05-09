/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
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

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;


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
  int i;
  double norm = 0.0;
  for(i=0;i<npw;i++){ norm += (std::conj(coeff[i]) * coeff[i] ).real();  }
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
  int i;
  if(kx==0 && ky==0 && kz==0){

    // Increase storage
    coeff.resize(2*npw-1);

    // Add new coefficients and compute norm of the 'completed' wavefunction
    double norm = (std::conj(coeff[0]) * coeff[0] ).real();
    for(i=1;i<npw;i++){
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





complex<double> I_1D(double kx, double kxp, double gx, double gxp){
/*
    # All arguments are float
    # kx,  gx  - refer to k-point 1 and its corresponding grid points
    # kxp, gxp - refer to k-point 2 and its corresponding grid points
*/

    double zero = 1e-8;
    complex<double> res(0.0, 0.0);
    complex<double> one(0.0, 1.0);

    double delt = kx + gx - kxp - gxp;

    if(fabs(delt) <= zero){   res = complex<double>(1.0, 0.0); }
    else{  res = -one * ( ( complex<double>( cos(2.0*M_PI*delt) - 1.0 ), sin(2.0*M_PI*delt)) / (2.0 * M_PI * delt) );  }

    return res;
}


complex<double> I_3D(VECTOR& k, VECTOR& kp, VECTOR& g, VECTOR& gp){
/*
    # All arguments are VECTOR (float)
    # k,  g  - refer to k-point 1 and its corresponding grid points
    # kp, gp - refer to k-point 2 and its corresponding grid points
*/

    complex<double> res = I_1D(k.x, kp.x, g.x, gp.x) * I_1D(k.y, kp.y, g.y, gp.y) * I_1D(k.z, kp.z, g.z, gp.z);

    return res;

}


CMATRIX pw_overlap(VECTOR& k1, VECTOR& k2, CMATRIX& coeff1, CMATRIX& coeff2, vector<VECTOR>& grid1, vector<VECTOR>& grid2){
/*
    # all k- and g-points are in units of 2*pi/a
    # k1, k2 - are k-point vectors (VECTOR of float) 
    # coeff1 - is a matrix (complex) of coefficeints for all states for given k-point (1), dimensions: npw1 x nbands1
    # coeff2 - is a matrix (complex) of coefficeints for all states for given k-point (2), dimensions: npw2 x nbands2
    # grid1 - a list of vectors for all G-points for given k-point (1): dimension npw1
    # grid2 - a list of vectors for all G-points for given k-point (2): dimension npw2
*/

    int npw1 = coeff1.n_rows;
    int nbands1 = coeff1.n_cols;

    int npw2 = coeff2.n_rows;
    int nbands2 = coeff2.n_cols;


    CMATRIX* S;
    S = new CMATRIX(nbands1, nbands2);  // all orbitals for given pair of k-points (a block of entire matrix)


    // A double sum over the grid points (may be different for the two k-points)
    for(int g1=0; g1<npw1; g1++){

        for(int g2=0; g2<npw2; g2++){

            complex<double> s = I_3D(k1, k2, grid1[g1], grid2[g2]);

            for(int i1=0; i1<nbands1; i1++){
                for(int i2=0; i2<nbands2; i2++){

                    complex<double> tmp = std::conj(coeff1.get(g1,i1)) * s * coeff2.get(g2,i2);

                    S->set(i1,i2, S->get(i1,i2) + tmp);

                }// for i2
            }// for i1
        }// for g2
    }// for g1

    return *S;

}

    


}// namespace libqobjects
}// namespace liblibra

