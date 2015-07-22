#include "Hamiltonian.h"
#include <stdlib.h>

namespace libhamiltonian{
namespace libhamiltonian_generic{


Hamiltonian::Hamiltonian(){ 
//cout<<"Base Ham. constructor\n";
}

Hamiltonian::~Hamiltonian(){ 
//cout<<"Base Ham. destructor\n"; 
}

void Hamiltonian::set_rep(int rep_){
  rep = rep_;
}


void Hamiltonian::set_params(boost::python::list params_){

  int sz = boost::python::len(params_);
  vector<double> tmp_params(sz, 0.0);

  // Now copy input params:
  for(int i=0;i<sz; i++){
    tmp_params[i] = boost::python::extract<double>(params_[i]);
  }

  set_params(tmp_params);

}

void Hamiltonian::set_q(vector<double>& q_){
  q = q_;
  status_dia = 0;
  status_adi = 0;
}

void Hamiltonian::set_q(boost::python::list q_){
 
  int sz = boost::python::len(q_);
  vector<double> tmp_q(sz,0.0);

  for(int i=0;i<sz; i++){
    tmp_q[i] = boost::python::extract<double>(q_[i]);
  }

  set_q(tmp_q);
}

void Hamiltonian::set_v(vector<double>& v_){
  v = v_;
  status_adi = 0;  // only affects adiabatic computations
}

void Hamiltonian::set_v(boost::python::list v_){

  int sz = boost::python::len(v_);
  vector<double> tmp_v(sz,0.0);

  for(int i=0;i<sz; i++){
    tmp_v[i] = boost::python::extract<double>(v_[i]);
  }

  set_v(tmp_v);

}


void Hamiltonian::compute(){

  if(rep==0){  compute_diabatic();  }
  else if(rep==1){  compute_adiabatic(); }

}



std::complex<double> Hamiltonian::H(int i,int j){

  std::complex<double> res(0.0,0.0);

  if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( ham_dia->get(i,j), 0.0 );
  }
  else if(rep==1){    // Adiabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( ham_adi->get(i,j), 0.0 );
  }

  return res;
}

std::complex<double> Hamiltonian::dHdq(int i,int j,int n){

  std::complex<double> res(0.0,0.0);

  if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( d1ham_dia[n]->get(i,j), 0.0 );
  }
  else if(rep==1){    // Adiabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( d1ham_adi[n]->get(i,j), 0.0 );
  }

  return res;
}


std::complex<double> Hamiltonian::D(int i,int j,int n){

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

std::complex<double> Hamiltonian::nac(int i,int j){

  std::complex<double> res(0.0,0.0);

  for(int n=0;n<nnucl;n++){
    res += D(i,j,n) * v[n]; 
  }
  return res;
}

std::complex<double> Hamiltonian::Hvib(int i,int j){
 
  const double hbar = 1.0;  // in atomic units

  std::complex<double> ham_ = H(i,j);
  std::complex<double> nac_ = nac(i,j);

  std::complex<double> res(ham_.real(), ham_.imag() - hbar* nac_.real() );

  return res;
}


}// namespace libhamiltonian_generic
}// namespace libhamiltonian
