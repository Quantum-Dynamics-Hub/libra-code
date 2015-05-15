#ifndef ELECTRONIC_H
#define ELECTRONIC_H

#include "../../mmath/libmmath.h"
#include "../../hamiltonian/libhamiltonian.h"

using namespace libmmath;
using namespace libmmath::librandom;
using namespace libhamiltonian;


namespace libdyn{
namespace libelectronic{



class Electronic{
  // This class describes a time-dependent state (multiconfigurational in sense that 
  // it can have time-dependent contributions from different stationary states )
  
  Random* rnd_obj;
  void rnd_phase(double&, double&, double);

  public:

  int nstates;             // number of stationary (basis) states
  int istate;              // index of current basis state (for stochastic)
  std::vector<double> q;   // MMTS variables for all basis states Re(c)
  std::vector<double> p;   // MMTS variables        --- // ---    Im(c)


  //--------------- Class functions ---------------

  void init(int,int);
  void init(int);

  Electronic();
  Electronic(int);
  Electronic(int,int);

  ~Electronic();

  std::complex<double> c(int i);      // return amplitude in the complex format: c_i = q_i + i*p_i
  std::complex<double> rho(int i, int j); // return the density matrix element: rho_ij = c^*_i * c_j


  //------ Methods ------------
  // In Electronic_Dynamics1.cpp
  void propagate_electronic(double dt,Hamiltonian* ham);
  void propagate_electronic(double dt,Hamiltonian& ham);

};


// In Electronic_Dynamics1.cpp
void propagate_electronic(double dt,Electronic* el,Hamiltonian* ham);


}// namespace libelectronic
}// namespace libdyn

#endif // ELECTRONIC_H
