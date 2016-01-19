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

#ifndef ELECTRONIC_H
#define ELECTRONIC_H

#include "../../mmath/libmmath.h"
#include "../../hamiltonian/libhamiltonian.h"

using namespace libmmath;
using namespace libmmath::librandom;
using namespace libhamiltonian;


namespace libdyn{
namespace libelectronic{


/*! This class describes a time-dependent state (multiconfigurational in sense that 
    it can have time-dependent contributions from different stationary states )

*/

class Electronic{
  
  Random* rnd_obj;                           
  void rnd_phase(double&, double&, double);

  public:

  int nstates;             ///< number of stationary (basis) states
  int istate;              ///< index of current basis state (for stochastic)
  std::vector<double> q;   ///< MMTS variables for all basis states: q = Re(c)
  std::vector<double> p;   ///< MMTS variables for all basis states: p = Im(c)


  //--------------- Class functions ---------------

  void init(int,int);      ///< initialize object: 1-st parameter - the number of electronic DOF, 2-nd - the initial state index
  void init(int);          ///< initialize object: parameter - the number of electronic DOF

  Electronic();
  Electronic(int);
  Electronic(int,int);
  Electronic(const Electronic&); ///< cctor

  ~Electronic();

  Electronic& operator=(const Electronic& ob);


  std::complex<double> c(int i);          ///< return amplitude in the complex format: c_i = q_i + i*p_i
  std::complex<double> rho(int i, int j); ///< return the density matrix element: rho_ij = c^*_i * c_j


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
void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib);
void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib, MATRIX& S);
void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib, CMATRIX& S);
void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib, CMATRIX& S, CMATRIX& Sdot);


}// namespace libelectronic
}// namespace libdyn

#endif // ELECTRONIC_H
