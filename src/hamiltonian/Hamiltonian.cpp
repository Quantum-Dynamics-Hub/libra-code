#include "Hamiltonian.h"
#include <stdlib.h>

namespace libhamiltonian{

Hamiltonian::Hamiltonian(){ /*status = 0; cout<<"Calling base function. Level is too abstract\n"; exit(0);*/}

Hamiltonian::~Hamiltonian(){ /*cout<<"Calling base function. Level is too abstract\n"; exit(0);*/ }


std::complex<double> Hamiltonian::H(int i, int j){
  std::cout<<"Calling base function. Level is too abstract\n";
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
