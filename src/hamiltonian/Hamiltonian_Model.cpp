/*********************************************************************************
* Copyright (C) 2014 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Hamiltonian_Model.h"
#include <complex>
#include <cmath>

namespace libhamiltonian{

using std::complex;
using std::sin;
using std::cos;
using std::exp;
using std::sqrt;


Hamiltonian_Model::Hamiltonian_Model(int ham_indx_){

  int i;
  ham_indx = ham_indx_;

  // 2-level models
  if(ham_indx==0 || ham_indx==1 || ham_indx==2 || ham_indx==3){  // SAC, DAC, ECWR, Marcus
    n_elec = 2;
    n_nucl = 1;
  }
  else if(ham_indx==4){    // superexchange
    n_elec = 3;
    n_nucl = 1;
  }

  ham_dia = new MATRIX(n_elec,n_elec); *ham_dia = 0.0;
  ham_adi = new MATRIX(n_elec,n_elec); *ham_adi = 0.0;

  d1ham_dia = vector<MATRIX*>(n_nucl);
  d2ham_dia = vector<MATRIX*>(n_nucl);
  for(i=0;i<n_nucl;i++){  d1ham_dia[i] = new MATRIX(n_elec,n_elec); *d1ham_dia[i] = 0.0;  }
  for(i=0;i<n_nucl*n_nucl;i++){  d2ham_dia[i] = new MATRIX(n_elec,n_elec); *d2ham_dia[i] = 0.0;  }
  

}


void Hamiltonian_Model::compute_diabatic(vector<double>& q,vector<double>& v){

  double e;

  //=============== 1D models ========================
  double x = q[0];

  // SAC potetnial
  if(ham_indx==0){
  }// ham_indx == 0  

  // DAC potential
  else if(ham_indx==1){
  }

  // ECWR potential
  else if(ham_indx==2){
  }


/*
  // Marcus spin-boson model
  else if(pot_indx==3){
  // conversion to Subotnik terms:
  // A = omega
  // B = E_r
  // C = V
  // D = eps_0
    double mo2 = 0.5*mol->mass[0] * A*A;
    double M = sqrt(mo2*B);

    H00 = mo2*x*x + M*x;     dH00 = 2.0*mo2*x + M;  d2H00 = 2.0*mo2;
    H11 = mo2*x*x - M*x - D; dH11 = 2.0*mo2*x - M;  d2H11 = 2.0*mo2;
    H01 = C;                 dH01 = 0.0;            d2H01 = 0.0;

  }

  // Superexchange model 
  else if(pot_indx==4){

    H00 = D1;  H01 = A1*exp(-B1*(x-C1)*(x-C1));  H02 = A3*exp(-B3*(x-C3)*(x-C3)); 
               H11 = D2;                         H12 = A2*exp(-B2*(x-C2)*(x-C2));  
                                                 H22 = D3;

    dH00 = 0.0;  dH01 = -2.0*B1*(x-C1)*H01;      dH02 = -2.0*B3*(x-C3)*H02; 
                 dH11 = 0.0;                     dH12 = -2.0*B2*(x-C2)*H12;
                                                 dH22 = 0.0;

    d2H00 = 0.0;  d2H01 = -2.0*B1*H01 + 4.0*B1*B1*(x-C1)*(x-C1)*H01;      d2H02 = -2.0*B3*H02 + 4.0*B3*B3*(x-C3)*(x-C3)*H02; 
                  d2H11 = 0.0;                                            d2H12 = -2.0*B2*H12 + 4.0*B2*B2*(x-C2)*(x-C2)*H12;
                                                                          d2H22 = 0.0;


  }
*/

}


/*
// Constructor
Hamiltonian_Tully::Hamiltonian_Tully(int rep_, int pot_indx_, Nuclear* mol){ 

  rep = rep_;
  set_model(pot_indx_);

  x = mol->R[0].x;
  v = mol->P[0].x/mol->mass[0];

  set_status(0);

}


// Choose Tully model and set corresponding parameters
void Hamiltonian_Tully::set_model(int pot_indx_){

  pot_indx = pot_indx_;
       if(pot_indx==0){   A = 0.010;  B = 1.600; C = 0.005; D = 1.00;           }
  else if(pot_indx==1){   A = 0.100;  B = 0.028; C = 0.015; D = 0.06; E = 0.05; } // alternative set
  else if(pot_indx==2){   A = 0.0006; B = 0.100; C = 0.900; D = 1.00;           }
  else if(pot_indx==3){   A = 3.5e-4; B = 2.39e-2; C = 1.49e-5; D = 1.5e-2;     }
  else if(pot_indx==4){
    D1=0.000;  D2=0.010;  D3=0.005;
    A1=0.001;  A2=0.010;  A3=0.000;
    B1=0.500;  B2=0.500;  B3=0.500;
    C1=0.000;  C2=0.000;  C3=0.000;
  }


}

void Hamiltonian_Tully::set_param(std::string var,double val){

  if(var=="A"){  A = val; }
  else if(var=="B"){  B = val; }
  else if(var=="C"){  C = val; }
  else if(var=="D"){  D = val; }
  else if(var=="E"){  E = val; }

}


// This function computed Hamiltonian and its derivatives in diabatic representation
// D_ij = <psi_i | d^2/dR^2 |psi_j> in addition to normal non-adiabatic coupling
// dd_ij = d/dR  d_ij 
// d_ij = <psi_i | d/dR | psi_j >
void Hamiltonian_Tully::compute_diabatic(Nuclear* mol){

    x = mol->R[0].x;

    double e;

    // SAC potetnial
    if(pot_indx==0){
  
      if(x>=0){ e = exp(-B*x); H00 = A*(1.0 - e); dH00 = A*B*e; d2H00 = -B*dH00; }
      else{  e = exp(B*x); H00 = -A*(1.0 - e); dH00 = A*B*e; d2H00 = -B*dH00; }
      H11 = -H00;           dH11 = -dH00;         d2H11 = -d2H00;
      H01 = C*exp(-D*x*x);  dH01 = -2.0*D*x*H01;  d2H01 = (-2.0*D*H01 - 2.0*D*x*dH01);

    }

    // DAC potential
    else if(pot_indx==1){

      H00  = 0.0; dH00 = 0.0;  d2H00 = 0.0;
      e = A*exp(-B*x*x); H11 = E - e; dH11 = 2.0*B*x*e; d2H11 = (2.0*B*e - 2.0*B*x*dH11);
      H01  = C*exp(-D*x*x); dH01 = -2.0*D*x*H01;  d2H01 = (-2.0*D*H01 - 2.0*D*x*dH01);

    }

    // ECWR potential
    else if(pot_indx==2){

      H00 = A;      dH00 = 0.0;  d2H00 = 0.0;
      H11 = -A;     dH11 = 0.0;  d2H11 = 0.0;
      if(x<=0){ H01 = B*exp(C*x); dH01 = C*H01; d2H01 = C*dH01; }
      else{ e = exp(-C*x); H01 = B*(2.0 - e); dH01 = B*C*e; d2H01 = -C*dH01; }
    }
 
    // Marcus spin-boson model
    else if(pot_indx==3){
    // conversion to Subotnik terms:
    // A = omega
    // B = E_r
    // C = V
    // D = eps_0
      double mo2 = 0.5*mol->mass[0] * A*A;
      double M = sqrt(mo2*B);

      H00 = mo2*x*x + M*x;     dH00 = 2.0*mo2*x + M;  d2H00 = 2.0*mo2;
      H11 = mo2*x*x - M*x - D; dH11 = 2.0*mo2*x - M;  d2H11 = 2.0*mo2;
      H01 = C;                 dH01 = 0.0;            d2H01 = 0.0;

    }

    // Superexchange model 
    else if(pot_indx==4){

      H00 = D1;  H01 = A1*exp(-B1*(x-C1)*(x-C1));  H02 = A3*exp(-B3*(x-C3)*(x-C3)); 
                 H11 = D2;                         H12 = A2*exp(-B2*(x-C2)*(x-C2));  
                                                   H22 = D3;

      dH00 = 0.0;  dH01 = -2.0*B1*(x-C1)*H01;      dH02 = -2.0*B3*(x-C3)*H02; 
                   dH11 = 0.0;                     dH12 = -2.0*B2*(x-C2)*H12;
                                                   dH22 = 0.0;

      d2H00 = 0.0;  d2H01 = -2.0*B1*H01 + 4.0*B1*B1*(x-C1)*(x-C1)*H01;      d2H02 = -2.0*B3*H02 + 4.0*B3*B3*(x-C3)*(x-C3)*H02; 
                    d2H11 = 0.0;                                            d2H12 = -2.0*B2*H12 + 4.0*B2*B2*(x-C2)*(x-C2)*H12;
                                                                            d2H22 = 0.0;


    }


}


void Hamiltonian_Tully::compute_adiabatic(Nuclear* mol){

    // First, do diabatic calculations, if they are not done yet
    compute_diabatic(mol);

    // Common operation - transformation from diabatic to adiabatic representation
    double sum,dsum,dif,ddif,dif2,ddif2,mix2,dmix2,sroot,dsroot,d2dif,dE;

    sum   = (H00 + H11);         dsum  = (dH00 + dH11);
    dif   = (H00 - H11);         ddif  = (dH00 - dH11);
                                 d2dif = (d2H00 - d2H11);

    dif2  = dif*dif;             ddif2 = 2.0*dif*ddif;
    mix2  = 4.0*H01*H01;         dmix2 = 8.0*H01*dH01;
    sroot = sqrt(dif2 + mix2);   dsroot      = 0.5*(ddif2 + dmix2)/sroot;

    E0 = 0.5*(sum - sroot); dE0 = 0.5*(dsum - dsroot);
    E1 = 0.5*(sum + sroot); dE1 = 0.5*(dsum + dsroot);
    dE = E0 - E1;


    // This is more correct (singularity-free) version
    d01 = ( ddif*H01 - dH01*dif )/(dif2 + mix2);
    d10 = -d01;

    // Symmetry: om2_01 = -om2_10
    double om2;
    om2 = (H01*(d2dif - 8.0*d01*dH01) - (d2H01 + 2.0*d01*ddif)*dif )/(dif2 + mix2);

    dd01 = om2;
    dd10 = -om2;

    D01 = (dif/dE)*om2 + 2.0*H01*d01*d01;
    D10 = D01;

  
}// compute_adiabatic



void Hamiltonian_Tully::compute(Nuclear* mol){

  if(get_status()==0){  
    if(rep==0){    compute_diabatic(mol); } //   set_status(status_diabatic); }
    else if(rep==1){ compute_adiabatic(mol); } // set_status(status_adiabatic); }
    
    set_status(1); 
  }

}



std::complex<double> Hamiltonian_Tully::H(Nuclear* mol,int i, int j){
// This function only returns computed values, one must call compute() first!!!

  if(get_status()==0){  compute_adiabatic(mol); set_status(1); }

  std::complex<double> res(0.0,0.0);

  if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian

    if(i==0 && j==0){ res = std::complex<double>(H00,0.0); }
    if(i==0 && j==1){ res = std::complex<double>(H01,0.0); }
    if(i==0 && j==2){ res = std::complex<double>(H02,0.0); }
    if(i==1 && j==0){ res = std::complex<double>(H01,0.0); }
    if(i==1 && j==1){ res = std::complex<double>(H11,0.0); }
    if(i==1 && j==2){ res = std::complex<double>(H12,0.0); }
    if(i==2 && j==0){ res = std::complex<double>(H02,0.0); }
    if(i==2 && j==1){ res = std::complex<double>(H12,0.0); }
    if(i==2 && j==2){ res = std::complex<double>(H22,0.0); }

  }// diabatic

  else if(rep==1){ // Adiabatic Hamiltonian - complex, imaginary part is antisymmetric => Hermitian

    // Diagonal terms - energies
    // off-diagonal terms -   (-i*hbar*NAC_ij * p/m), so this is already a vibronic Hamiltonian 
    v = mol->P[0].x/mol->mass[0];

    if(i==0 && j==0){ res = std::complex<double>(E0, 0.0); }
    if(i==0 && j==1){ res = std::complex<double>(0.0,-v*d01); }
    if(i==1 && j==0){ res = std::complex<double>(0.0, v*d01); }
    if(i==1 && j==1){ res = std::complex<double>(E1,0.0); }

    if(i>=2 || j>=2){  cout<<"Adiabatic representation available only for 2D models\n"; exit(0); }

  }// adiabatic

  return res;

}

std::complex<double> Hamiltonian_Tully::H(Nuclear* mol,int i, int j, int rep_){
// This function only returns computed values, one must call compute() first!!!

  if(get_status()==0){  compute_adiabatic(mol); set_status(1); }

  std::complex<double> res(0.0,0.0);

  if(rep_==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian

    if(i==0 && j==0){ res = std::complex<double>(H00,0.0); }
    if(i==0 && j==1){ res = std::complex<double>(H01,0.0); }
    if(i==0 && j==2){ res = std::complex<double>(H02,0.0); }
    if(i==1 && j==0){ res = std::complex<double>(H01,0.0); }
    if(i==1 && j==1){ res = std::complex<double>(H11,0.0); }
    if(i==1 && j==2){ res = std::complex<double>(H12,0.0); }
    if(i==2 && j==0){ res = std::complex<double>(H02,0.0); }
    if(i==2 && j==1){ res = std::complex<double>(H12,0.0); }
    if(i==2 && j==2){ res = std::complex<double>(H22,0.0); }

  }// diabatic

  else if(rep_==1){ // Adiabatic Hamiltonian - complex, imaginary part is antisymmetric => Hermitian

    // Diagonal terms - energies
    // off-diagonal terms -   (-i*hbar*NAC_ij * p/m), so this is already a vibronic Hamiltonian 
    v = mol->P[0].x/mol->mass[0];

    if(i==0 && j==0){ res = std::complex<double>(E0, 0.0); }
    if(i==0 && j==1){ res = std::complex<double>(0.0,-v*d01); }
    if(i==1 && j==0){ res = std::complex<double>(0.0, v*d01); }
    if(i==1 && j==1){ res = std::complex<double>(E1,0.0); }

    if(i>=2 || j>=2){  cout<<"Adiabatic representation available only for 2D models\n"; exit(0); }


  }// adiabatic

  return res;

}


std::complex<double> Hamiltonian_Tully::Dx(Nuclear* mol,int i, int j, int k){
// This function only returns computed values, one must call compute() first!!!
// In principle, nothing prevents the derivative coupling to be complex (for each component)

  if(get_status()==0){  compute_adiabatic(mol); set_status(1); }
  std::complex<double> res(0.0,0.0);

  if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian

    // In diabatic basis derivative couplings are zero, by definition

  }// diabatic

  else if(rep==1){ // Adiabatic Hamiltonian - complex, imaginary part is antisymmetric => Hermitian

    if(i==0 && j==0){ res = std::complex<double>( 0.0, 0.0); }
    if(i==0 && j==1){ res = std::complex<double>( d01, 0.0); }
    if(i==1 && j==0){ res = std::complex<double>(-d01, 0.0); }
    if(i==1 && j==1){ res = std::complex<double>( 0.0, 0.0); }

    if(i>=2 || j>=2){  cout<<"Adiabatic representation available only for 2D models\n"; exit(0); }


  }// adiabatic

  return res;

}

std::complex<double> Hamiltonian_Tully::Dy(Nuclear* mol,int i, int j, int k){
// This function only returns computed values, one must call compute() first!!!

  std::complex<double> res(0.0,0.0);

  return res;
}

std::complex<double> Hamiltonian_Tully::Dz(Nuclear* mol,int i, int j, int k){
// This function only returns computed values, one must call compute() first!!!

  std::complex<double> res(0.0,0.0);

  return res;
}




std::complex<double> Hamiltonian_Tully::dHdRx(Nuclear* mol,int i, int j,int k){
// This function only returns computed values, one must call compute() first!!!

  if(get_status()==0){  compute_adiabatic(mol); set_status(1); }
  std::complex<double> res(0.0,0.0);

  if(k==0){

    if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian
    
      if(i==0 && j==0){ res = std::complex<double>(dH00,0.0); }
      if(i==0 && j==1){ res = std::complex<double>(dH01,0.0); }
      if(i==0 && j==2){ res = std::complex<double>(dH02,0.0); }
      if(i==1 && j==0){ res = std::complex<double>(dH01,0.0); }
      if(i==1 && j==1){ res = std::complex<double>(dH11,0.0); }
      if(i==1 && j==2){ res = std::complex<double>(dH12,0.0); }
      if(i==2 && j==0){ res = std::complex<double>(dH02,0.0); }
      if(i==2 && j==1){ res = std::complex<double>(dH12,0.0); }
      if(i==2 && j==2){ res = std::complex<double>(dH22,0.0); }
    
    }// diabatic
    
    else if(rep==1){ // Adiabatic Hamiltonian - complex, imaginary part is antisymmetric => Hermitian
      v = mol->P[0].x/mol->mass[0];
    
      if(i==0 && j==0){ res = std::complex<double>(dE0, 0.0); }
      if(i==0 && j==1){ res = std::complex<double>(0.0,-v*dd01); }
      if(i==1 && j==0){ res = std::complex<double>(0.0, v*dd01); }
      if(i==1 && j==1){ res = std::complex<double>(dE1,0.0); }

      if(i>=2 || j>=2){  cout<<"Adiabatic representation available only for 2D models\n"; exit(0); }

    
    }// adiabatic

  }

  return res;

}

std::complex<double> Hamiltonian_Tully::dHdRy(Nuclear* mol,int i, int j,int k){
// This function only returns computed values, one must call compute() first!!!

  std::complex<double> res(0.0,0.0);
  return res;
}

std::complex<double> Hamiltonian_Tully::dHdRz(Nuclear* mol,int i, int j,int k){
// This function only returns computed values, one must call compute() first!!!

  std::complex<double> res(0.0,0.0);
  return res;
}



*/

}// namespace libhamiltonian




