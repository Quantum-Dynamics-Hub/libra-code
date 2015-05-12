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
#include "../mmath_eigen/libmmath_eigen.h"
using namespace libmmath;

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
  d1ham_adi = vector<MATRIX*>(n_nucl);
  for(i=0;i<n_nucl;i++){  
    d1ham_dia[i] = new MATRIX(n_elec,n_elec); *d1ham_dia[i] = 0.0; 
    d1ham_adi[i] = new MATRIX(n_elec,n_elec); *d1ham_adi[i] = 0.0; 
  }
  for(i=0;i<n_nucl*n_nucl;i++){  d2ham_dia[i] = new MATRIX(n_elec,n_elec); *d2ham_dia[i] = 0.0;  }


  rep = 0; // default representation is diabatic  

}

Hamiltonian_Model::~Hamiltonian_Model(){
  int i;
 
  delete ham_dia;
  delete ham_adi;

  for(i=0;i<n_nucl;i++){  
    delete d1ham_dia[i];
    delete d1ham_adi[i];
  }
  for(i=0;i<n_nucl*n_nucl;i++){  delete d2ham_dia[i]; }

  d1ham_dia.clear();
  d2ham_dia.clear();
  d1ham_adi.clear();



}

void Hamiltonian_Model::set_rep(int rep_){
  rep = rep_;
}

void Hamiltonian_Model::set_params(vector<double>& params_){

  int num_params = 0;

  if(ham_indx==0 || ham_indx==2){ // SAC, ECWR
    num_params = 4;
  }
  else if(ham_indx==1 || ham_indx==3){  // DAC, Marcus
    num_params = 5;
  }
  else if(ham_indx==4){  // SEXCH
    num_params = 12;
  }

  params = vector<double>(num_params, 0.0);

  // Now copy input params:
  for(int i=0;i<params_.size() && i<num_params; i++){
    params[i] = params_[i];
  }

}

void Hamiltonian_Model::set_params(boost::python::list params_){

  int num_params = 0;

  if(ham_indx==0 || ham_indx==2){ // SAC, ECWR
    num_params = 4;
  }
  else if(ham_indx==1 || ham_indx==3){  // DAC, Marcus
    num_params = 5;
  }
  else if(ham_indx==4){  // SEXCH
    num_params = 12;
  }

  params = vector<double>(num_params, 0.0);

  // Now copy input params:
  for(int i=0;i<boost::python::len(params_) && i<num_params; i++){
    params[i] = boost::python::extract<double>(params_[i]);
  }

}


void Hamiltonian_Model::compute_diabatic(vector<double>& q,vector<double>& v){

  double e;

  //=============== 1D models ========================
  double x = q[0];

  if(ham_indx==0){  SAC_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);     }    // SAC potetnial
  else if(ham_indx==1){ DAC_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params); }    // DAC potential
  else if(ham_indx==2){ ECWR_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);   } // ECWR potential
  else if(ham_indx==3){ Marcus_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params); } // Marcus potential
  else if(ham_indx==4){ SEXCH_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);  } // SEXCH potential

}

void Hamiltonian_Model::compute_diabatic(boost::python::list q_, boost::python::list v_){

  int i;
  int sz = boost::python::len(q_);

  vector<double> q(sz,0.0);
  for(i=0;i<sz; i++){
    q[i] = boost::python::extract<double>(q_[i]);
  }
 
  sz = boost::python::len(v_);
  vector<double> v(sz,0.0);
  for(i=0;i<sz; i++){
    v[i] = boost::python::extract<double>(v_[i]);
  }

  compute_diabatic(q,v);

//  cout<<"in compute_diabatic: \n";
//  cout<<"ham_dia = \n";
//  cout<<*ham_dia<<endl;

}


void Hamiltonian_Model::compute_adiabatic(vector<double>& q,vector<double>& v){
// This function computes adiabatic PESs and derivative couplings

  compute_diabatic(q,v);

  MATRIX* S; S = new MATRIX(n_elec, n_elec);  S->Init_Unit_Matrix(1.0);
  MATRIX* C; C = new MATRIX(n_elec, n_elec);  *C = 0.0;

  // Transformation to adiabatic basis
  solve_eigen(n_elec, ham_dia, S, ham_adi, C);  // H_dia * C = S * C * H_adi


  // Now compute the derivative couplings (off-diagonal, multiplied by energy difference) and adiabatic gradients (diagonal)
  for(int n=0;n<n_nucl;n++){

    *d1ham_adi[n] = (*C).T() * (*d1ham_dia[n]) * (*C);

  }// for n


  delete S;
  delete C;
  
}

void Hamiltonian_Model::compute_adiabatic(boost::python::list q_, boost::python::list v_){

  int i;
  int sz = boost::python::len(q_);

  vector<double> q(sz,0.0);
  for(i=0;i<sz; i++){
    q[i] = boost::python::extract<double>(q_[i]);
  }
 
  sz = boost::python::len(v_);
  vector<double> v(sz,0.0);
  for(i=0;i<sz; i++){
    v[i] = boost::python::extract<double>(v_[i]);
  }

  compute_adiabatic(q,v);

}


std::complex<double> Hamiltonian_Model::H(int i,int j){

  std::complex<double> res(0.0,0.0);

  if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( ham_dia->get(i,j), 0.0 );
  }
  else if(rep==1){    // Adiabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( ham_adi->get(i,j), 0.0 );
  }

  return res;
}

/*


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




