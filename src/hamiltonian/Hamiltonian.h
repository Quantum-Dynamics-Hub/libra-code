/*********************************************************************************
* Copyright (C) 2012 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <complex>
#include "../mmath/libmmath.h"
using namespace libmmath;


namespace libhamiltonian{


class Hamiltonian{
//
// This is our functor class: it manages customized Hamiltonian calculations
//                            customization is realized via inheritance
//                            common interface is realized via virtual functions
//

// This class implements an interface with Hamiltonian calculations and other associated tasks
// for example getting it from quantum simulations or reading this information from files
//
// It  impractical to store entire Hamiltonian and, especially all its derivatives at once
// Example: for 100 C atoms with 4 basis functions on each there are 400 basis states, so 
// the 400 x 400 x 2 x 8 = 2560000 bytes = 2.56 Mb are required for Ham, and 
// 2.56 x 100 = 256 Mb for each of the derivative
// These numbers go far beyond practical limits for somewhat larger systems of ~1000 atoms
// so it may be better to get all the elements algorithmically, although it can be slower
//
//
// The commented members below is what we would want ideally, but in practice we have to find better way
//  matrix* Ham;               // Electronic energy levels and couplings
//  vector<matrix* > dHamdRx;  // derivatives of Ham wrt X coordinate of all atoms
//  vector<matrix* > dHamdRy;  // derivatives of Ham wrt Y coordinate of all atoms
//  vector<matrix* > dHamdRz;  // derivatives of Ham wrt Z coordinate of all atoms


//  To solve this problem we use virtual functions for computation of H and dHdR 
//  specific implementation is controlled by the derived classes

//  int status;   //   1 - if Hamiltonian is up to date for position and velocity (velocity enters couplings)

public:

  // Constructors
  Hamiltonian(); //{ ;; }

  // Use dfault copy constructor
  //Hamiltonian(const Hamiltonian&);

  // Destructor
  virtual ~Hamiltonian();    // This destructor does nothing - but in fact this is intended so the el and mol objects are not 
                             // destroyed, and the references are copied verbatim (not their content) - this helps prevent 
                             // additional overhead due to construction of el and mol objects

//  void set_status(int st_){ status = st_; }
//  int get_status(){ return status; }


  // NOTE!!! It is very important to provide a minimal implementation of virtual methods
  // in the source file (Hamiltonian.cpp)
  // Otherwise the symbols remain undefined in the resulting library
  // This is a nasty bug to be aware of. 
  // See for example: http://stackoverflow.com/questions/1458180/vtable-for-referenced-from-compile-error-xcode/1478553#1478553
  // Updates coordinates, velocities, etc. - all nuclear information
//  virtual void update_nuclear(Nuclear* mol);

  // This function performs actual computations
  virtual void compute(){ ;; }

//  virtual void set_param(std::string var, double val){  ;; }

  // Calculation methods
  virtual std::complex<double> H(int, int);          // Hamiltonian
  virtual std::complex<double> dHdRx(int, int, int); // derivative of Hamiltonian w.r.t. Cartesian coordinate x
  virtual std::complex<double> dHdRy(int, int, int); // derivative of Hamiltonian w.r.t. Cartesian coordinate y
  virtual std::complex<double> dHdRz(int, int, int); // derivative of Hamiltonian w.r.t. Cartesian coordinate z


// Nonadiabatic couplings are not the part of the Hamiltonian functor - they should be related to representation !
/*
  virtual std::complex<double> D1x(int, int, int); // first-order derivative coupling w.r.t. Cartesian coordinate x
  virtual std::complex<double> D1y(int, int, int); // first-order derivative coupling w.r.t. Cartesian coordinate y
  virtual std::complex<double> D1z(int, int, int); // first-order derivative coupling w.r.t. Cartesian coordinate z

  virtual std::complex<double> D2x(int, int, int); // second-order derivative coupling w.r.t. Cartesian coordinate x
  virtual std::complex<double> D2y(int, int, int); // second-order derivative coupling w.r.t. Cartesian coordinate y
  virtual std::complex<double> D2z(int, int, int); // second-order derivative coupling w.r.t. Cartesian coordinate z
*/

};


}// namespace libhamiltonian

#endif // HAMILTONIAN_H
