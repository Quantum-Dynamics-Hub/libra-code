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
/**
  \file Electronic.h
  \brief The file describes Electronic class and functions for propagation of electronic DOF
    
*/

#ifndef ELECTRONIC_H
#define ELECTRONIC_H

#include "../../mmath/libmmath.h"
#include "../../hamiltonian/libhamiltonian.h"

using namespace libmmath;
using namespace libmmath::librandom;
using namespace libhamiltonian;

/// libdyn namespace 
namespace libdyn{

/// libelectronic namespace 
namespace libelectronic{


/*! This class describes a time-dependent state (multiconfigurational in sense that 
    it can have time-dependent contributions from different stationary states )

*/

class Electronic{
/** 
  \brief Electronic class that is designed to represent electronic degrees of freedom

  The electronic DOF in this case are just the (time-dependent) coefficients of the 
  basis functions in the representation of the overall time-dependent wavefunction:
  |PSI> = sum_i { c_i(t) * |i>} , so c_i - is the electronic DOF

  In this case, we use the Meyer-Miller-Thoss-Stock (MMTS) classically-mapped coordinates: q = Re(c), p = Im(c)

  See more: 
  (1) Meyer, H.-D.; Miller, W. H. A Classical Analog for Electronic Degrees of Freedom in Nonadiabatic Collision Processes. J. Chem. Phys. 1979, 70, 3214–3223.
  (2) Meyer, H.-D.; Miller, W. H. Analysis and Extension of Some Recently Proposed Classical Models for Electronic Degrees of Freedom. J. Chem. Phys. 1980, 72, 2272–2281.
  (3) Thoss, M.; Stock, G. Mapping Approach to the Semiclassical Description of Nonadiabatic Quantum Dynamics. Phys. Rev. A 1999, 59, 64–79.
  (4) Stock, G.; Thoss, M. Semiclassical Description of Nonadiabatic Quantum Dynamics. Phys. Rev. Lett. 1997, 78, 578–581.

*/

  
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
  Electronic(const Electronic&); 

  ~Electronic();

  Electronic& operator=(const Electronic& ob);


  std::complex<double> c(int i);          ///< return amplitude in the complex format: c_i = q_i + i*p_i
  std::complex<double> rho(int i, int j); ///< return the density matrix element: rho_ij = c^*_i * c_j


  //------ Methods ------------
  // In Electronic_Dynamics1.cpp
  void propagate_electronic(double dt,Hamiltonian* ham);
  void propagate_electronic(double dt,Hamiltonian& ham);
  void propagate_electronic(double dt,Hamiltonian& ham, CMATRIX& S);


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


typedef std::vector< Electronic > ElectronicList; ///< Type containing the vector of Electronic objects


// In Electronic_Dynamics1.cpp
void propagate_electronic(double dt,Electronic* el,Hamiltonian* ham);
void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib);
void propagate_electronic(double dt,CMATRIX& Coeff, CMATRIX& Hvib);

void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib, MATRIX& S);
void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib, CMATRIX& S);
void propagate_electronic(double dt,CMATRIX& Coeff, CMATRIX& Hvib, CMATRIX& S);

void grid_propagator(double dt, CMATRIX& Hvib, CMATRIX& S, CMATRIX& U);


}// namespace libelectronic

}// namespace libdyn

#endif // ELECTRONIC_H
