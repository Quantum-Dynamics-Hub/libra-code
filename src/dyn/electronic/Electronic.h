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
  Electronic(const Electronic&); // cctor

  ~Electronic();

  Electronic& operator=(const Electronic& ob);


  std::complex<double> c(int i);      // return amplitude in the complex format: c_i = q_i + i*p_i
  std::complex<double> rho(int i, int j); // return the density matrix element: rho_ij = c^*_i * c_j


  //------ Methods ------------
  // In Electronic_Dynamics1.cpp
  void propagate_electronic(double dt,Hamiltonian* ham);
  void propagate_electronic(double dt,Hamiltonian& ham);


  friend bool operator == (const Electronic& e1, const Electronic& e2){
/*
    bool res = ( (e1.istate == e2.istate) && (e1.nstates == e2.nstates)  );
    for(int i=0;i<e1.nstates;i++){  res = res && (e1.q[i] == e2.q[i]) && (e1.p[i] == e2.p[i]); }
    return res;
*/  return &e1 == &e2;
  }
  friend bool operator != (const Electronic& e1, const Electronic& e2){
    return !(e1==e2);  // only compare addresses
  }


};


typedef std::vector< Electronic > ElectronicList;


// In Electronic_Dynamics1.cpp
void propagate_electronic(double dt,Electronic* el,Hamiltonian* ham);


}// namespace libelectronic
}// namespace libdyn

#endif // ELECTRONIC_H
