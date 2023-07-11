/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
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

#include "../../math_linalg/liblinalg.h"
#include "../../math_random/librandom.h"
#include "../../nhamiltonian/libnhamiltonian.h"
/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace librandom;
using namespace libnhamiltonian;


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

  
  void rnd_phase(double&, double&, double, double);

  public:

  int nstates;             ///< number of stationary (basis) states
  int istate;              ///< index of current basis state (for stochastic)
  std::vector<double> q;   ///< MMTS variables for all basis states: q = Re(c)
  std::vector<double> p;   ///< MMTS variables for all basis states: p = Im(c)


  //--------------- Class functions ---------------

  void init(int,int,double);///< initialize object: 1-st parameter - the number of electronic DOF, 2-nd - the initial state index, 3-rd - the phase
  void init(int,int);       ///< initialize object: 1-st parameter - the number of electronic DOF, 2-nd - the initial state index
  void init(int);           ///< initialize object: parameter - the number of electronic DOF

  Electronic();
  Electronic(int);
  Electronic(int,int);
  Electronic(int,int,double);
  Electronic(const Electronic&); 

  ~Electronic();

  Electronic& operator=(const Electronic& ob);


  std::complex<double> c(int i) const;          ///< return amplitude in the complex format: c_i = q_i + i*p_i
  std::complex<double> rho(int i, int j) const; ///< return the density matrix element: rho_ij = c^*_i * c_j
  CMATRIX C() const;   ///< return the amplitutes in a vector format
  CMATRIX RHO() const; ///< Return the density matrix in a matrix format



  //------ Methods ------------
  // In Electronic_Dynamics1.cpp
  void project_out(int i);
  void project_out(int i, int renorm_flag);
  void collapse(int i);
  void collapse(int i, int phase_flag);

//  void propagate_electronic(double dt,Hamiltonian* ham);
//  void propagate_electronic(double dt,Hamiltonian& ham);
//  void propagate_electronic(double dt,Hamiltonian& ham, CMATRIX& S);


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


//=================== In Electronic_Dynamics1.cpp ======================
// Auxiliary functions for rotations-based integrator
void iL2_action(double dt, CMATRIX& Coeff, CMATRIX& Hvib, int i, int j);
void iL3_action(double dt, CMATRIX& Coeff, CMATRIX& Hvib, int i, int j);

// Elementary integration algorithms
void propagate_electronic_rot(double dt, CMATRIX& Coeff, CMATRIX& Hvib);
void propagate_electronic_eig(double dt, CMATRIX& Coeff, CMATRIX& Hvib);
void propagate_electronic_eig(double dt, CMATRIX& Coeff, CMATRIX& Hvib, CMATRIX& S);
void propagate_electronic_nonHermitian(double dt, CMATRIX& Coeff, CMATRIX& Hvib);
void propagate_electronic_qtag(double dt, CMATRIX& Coeff, CMATRIX& Hvib, CMATRIX& S);
void propagate_electronic_qtag2(double dt, CMATRIX& Coeff, CMATRIX& Hvib, CMATRIX& Hvib_old, CMATRIX& S, CMATRIX& S_old);


// For grid integration
void grid_propagator(double dt, CMATRIX& Hvib, CMATRIX& S, CMATRIX& U);

// For Liouvillian integration
CMATRIX vectorize_density_matrix(CMATRIX* rho);
CMATRIX vectorize_density_matrix(CMATRIX& rho);
CMATRIX unvectorize_density_matrix(CMATRIX& rho_vec);
CMATRIX make_Liouvillian(CMATRIX& ham);


}// namespace libelectronic

}// namespace libdyn
}// liblibra

#endif // ELECTRONIC_H
