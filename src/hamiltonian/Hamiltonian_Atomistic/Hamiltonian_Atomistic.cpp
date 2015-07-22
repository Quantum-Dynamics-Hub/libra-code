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


#include <complex>
#include <cmath>
#include "Hamiltonian_Atomistic.h"


namespace libhamiltonian{
namespace libhamiltonian_atomistic{


using namespace libmmath;
using namespace libmmath::libmeigen;
using std::complex;
using std::sin;
using std::cos;
using std::exp;
using std::sqrt;


Hamiltonian_Atomistic::Hamiltonian_Atomistic(int _nelec, int _nnucl){

  int i;

  nelec = _nelec;
  nnucl = _nnucl;

  ham_dia = new MATRIX(nelec,nelec); *ham_dia = 0.0;
  ham_adi = new MATRIX(nelec,nelec); *ham_adi = 0.0;

  d1ham_dia = vector<MATRIX*>(nnucl);
  d2ham_dia = vector<MATRIX*>(nnucl);
  d1ham_adi = vector<MATRIX*>(nnucl);
  for(i=0;i<nnucl;i++){  
    d1ham_dia[i] = new MATRIX(nelec,nelec); *d1ham_dia[i] = 0.0; 
    d1ham_adi[i] = new MATRIX(nelec,nelec); *d1ham_adi[i] = 0.0; 
  }
  for(i=0;i<nnucl*nnucl;i++){  d2ham_dia[i] = new MATRIX(nelec,nelec); *d2ham_dia[i] = 0.0;  }


  rep = 0; // default representation is diabatic  
  status_dia = 0;
  status_adi = 0;

}

Hamiltonian_Atomistic::~Hamiltonian_Atomistic(){
  int i;

  delete ham_dia;
  delete ham_adi;

  for(i=0;i<nnucl;i++){  
    delete d1ham_dia[i];
    delete d1ham_adi[i];
  }
  for(i=0;i<nnucl*nnucl;i++){  delete d2ham_dia[i]; }

  if(d1ham_dia.size()>0) { d1ham_dia.clear(); }
  if(d2ham_dia.size()>0) { d2ham_dia.clear(); }
  if(d1ham_adi.size()>0) { d1ham_adi.clear(); }


}

void Hamiltonian_Atomistic::set_rep(int rep_){
  rep = rep_;
}

void Hamiltonian_Atomistic::set_params(vector<double>& params_){

  
  params = vector<double>(params_.size(), 0.0);

  // Now copy input params:
  for(int i=0;i<params_.size(); i++){
    params[i] = params_[i];
  }

  // Since the parameters have changed - we need to recompute everything
  status_dia = 0;
  status_adi = 0;

}

void Hamiltonian_Atomistic::set_params(boost::python::list params_){

  int sz = boost::python::len(params_);
  vector<double> tmp_params(sz, 0.0);

  // Now copy input params:
  for(int i=0;i<sz; i++){
    tmp_params[i] = boost::python::extract<double>(params_[i]);
  }

  set_params(tmp_params);

}

void Hamiltonian_Atomistic::set_q(vector<double>& q_){
  q = q_;
  status_dia = 0;
  status_adi = 0;
}

void Hamiltonian_Atomistic::set_q(boost::python::list q_){
 
  int sz = boost::python::len(q_);
  vector<double> tmp_q(sz,0.0);

  for(int i=0;i<sz; i++){
    tmp_q[i] = boost::python::extract<double>(q_[i]);
  }

  set_q(tmp_q);
}

void Hamiltonian_Atomistic::set_v(vector<double>& v_){
  v = v_;
  status_adi = 0;  // only affects adiabatic computations
}

void Hamiltonian_Atomistic::set_v(boost::python::list v_){

  int sz = boost::python::len(v_);
  vector<double> tmp_v(sz,0.0);

  for(int i=0;i<sz; i++){
    tmp_v[i] = boost::python::extract<double>(v_[i]);
  }

  set_v(tmp_v);

}


void Hamiltonian_Atomistic::compute(){

  if(rep==0){  compute_diabatic();  }
  else if(rep==1){  compute_adiabatic(); }

}


void Hamiltonian_Atomistic::compute_diabatic(){

  if(status_dia == 0){ // only compute this is the result is not up to date
  
    //==========================================================

    // class - specific computations: MM, semi-empirical, HF, etc   

    //==========================================================

    // Set status flag 
    status_dia = 1;

  }//   status_dia == 0

}


void Hamiltonian_Atomistic::compute_adiabatic(){
// This function computes adiabatic PESs and derivative couplings

//------------------------------------------------------------------

// actually, here is gonna be SCF solver - to get adiabatic states

//------------------------------------------------------------------

  compute_diabatic();

  if(status_adi == 0){

    MATRIX* S; S = new MATRIX(nelec, nelec);  S->Init_Unit_Matrix(1.0);
    MATRIX* C; C = new MATRIX(nelec, nelec);  *C = 0.0;

    // Transformation to adiabatic basis
    solve_eigen(nelec, ham_dia, S, ham_adi, C);  // H_dia * C = S * C * H_adi


    // Now compute the derivative couplings (off-diagonal, multiplied by energy difference) and adiabatic gradients (diagonal)
    for(int n=0;n<nnucl;n++){

      *d1ham_adi[n] = (*C).T() * (*d1ham_dia[n]) * (*C);

    }// for n


    delete S;
    delete C;

    // Set status flag
    status_adi = 1;

  }// status_adi == 0
  
}


std::complex<double> Hamiltonian_Atomistic::H(int i,int j){

  std::complex<double> res(0.0,0.0);

  if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( ham_dia->get(i,j), 0.0 );
  }
  else if(rep==1){    // Adiabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( ham_adi->get(i,j), 0.0 );
  }

  return res;
}

std::complex<double> Hamiltonian_Atomistic::dHdq(int i,int j,int n){

  std::complex<double> res(0.0,0.0);

  if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( d1ham_dia[n]->get(i,j), 0.0 );
  }
  else if(rep==1){    // Adiabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( d1ham_adi[n]->get(i,j), 0.0 );
  }

  return res;
}


std::complex<double> Hamiltonian_Atomistic::D(int i,int j,int n){

  std::complex<double> res(0.0,0.0);

  if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>(0.0,0.0);
  }
  else if(rep==1){    // Adiabatic Hamiltonian - real, symmetric => Hermitian
    if(i!=j){  

      double dE = (ham_adi->get(j,j) - ham_adi->get(i,i) );
      if(fabs(dE)<1e-10){ dE = 1e-10 * (dE>0.0 ? 1.0 : -1.0); }

      res = std::complex<double>( d1ham_adi[n]->get(i,j)/dE, 0.0 );

    }
  }

  return res;
}

std::complex<double> Hamiltonian_Atomistic::nac(int i,int j){

  std::complex<double> res(0.0,0.0);

  for(int n=0;n<nnucl;n++){
    res += D(i,j,n) * v[n]; 
  }
  return res;
}

std::complex<double> Hamiltonian_Atomistic::Hvib(int i,int j){
 
  const double hbar = 1.0;  // in atomic units

  std::complex<double> ham_ = H(i,j);
  std::complex<double> nac_ = nac(i,j);

  std::complex<double> res(ham_.real(), ham_.imag() - hbar* nac_.real() );

  return res;
}



}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian




