#include "Hamiltonian.h"
#include <stdlib.h>

namespace libhamiltonian{

Hamiltonian::Hamiltonian(){ /*cout<<"Base Ham. constructor\n";*/ }

Hamiltonian::~Hamiltonian(){ /*cout<<"Base Ham. destructor\n";*/ }


std::complex<double> Hamiltonian::H(int i, int j){
  std::cout<<"Calling the base function H. Level is too abstract\n";
  exit(0);
  return std::complex<double>(0.0,0.0); 
}

std::complex<double> Hamiltonian::dHdq(int i, int j, int n){
  std::cout<<"Calling the base function dHdq. Level is too abstract\n";
  exit(0);
  return std::complex<double>(0.0,0.0); 
}

std::complex<double> Hamiltonian::D(int i, int j, int n){
  std::cout<<"Calling the base function D. Level is too abstract\n";
  exit(0);
  return std::complex<double>(0.0,0.0); 
}

std::complex<double> Hamiltonian::nac(int i, int j){
  std::cout<<"Calling the base function nac. Level is too abstract\n";
  exit(0);
  return std::complex<double>(0.0,0.0); 
}

std::complex<double> Hamiltonian::Hvib(int i, int j){
  std::cout<<"Calling the base function Hvib. Level is too abstract\n";
  exit(0);
  return std::complex<double>(0.0,0.0); 
}



std::complex<double> Hamiltonian::dHdRx(int i, int j, int k){
  std::cout<<"Calling base function dHdRx. Level is too abstract\n"; 
  exit(0);
  return std::complex<double>(0.0,0.0); 
}

std::complex<double> Hamiltonian::dHdRy(int i, int j, int k){
  std::cout<<"Calling base function dHdRy. Level is too abstract\n"; 
  exit(0);
  return std::complex<double>(0.0,0.0);
}

std::complex<double> Hamiltonian::dHdRz(int i, int j, int k){
  std::cout<<"Calling base function dHdRz. Level is too abstract\n"; 
  exit(0);
  return std::complex<double>(0.0,0.0);
}

/*
std::complex<double> Hamiltonian::D1x(int i, int j, int k){
  std::cout<<"Calling base function D1x. Level is too abstract\n"; 
  exit(0);
  return std::complex<double>(0.0,0.0); 
}

std::complex<double> Hamiltonian::D1y(int i, int j, int k){
  std::cout<<"Calling base function D1y. Level is too abstract\n"; 
  exit(0);
  return std::complex<double>(0.0,0.0);
}

std::complex<double> Hamiltonian::D1z(int i, int j, int k){
  std;:cout<<"Calling base function D1z. Level is too abstract\n";
  exit(0);
  return std::complex<double>(0.0,0.0);
}
*/


}// namespace libhamiltonian
