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
  \file Hamiltonian.h
  \brief The file describes the generic Hamiltonian class
    
*/

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <complex>
#endif

#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{


/// libhamiltonian namespace
namespace libatomistic{

/// libhamiltonian_generic namespace
namespace libhamiltonian_generic{


using namespace liblinalg;

class Hamiltonian{
/**
  This is our functor class: it manages customized Hamiltonian calculations
                             customization is realized via inheritance
                             common interface is realized via virtual functions

  This class implements an interface with Hamiltonian calculations and other associated tasks
  for example, getting it from quantum simulations or reading this information from files

  It is impractical to store the entire Hamiltonian and, especially all its derivatives at once

  Example: for 100 C atoms with 4 basis functions on each there are 400 basis states, so 
  the 400 x 400 x 2 x 8 = 2560000 bytes = 2.56 Mb are required for Ham, and 
  2.56 x 100 = 256 Mb for each of the derivative
  These numbers go far beyond practical limits for somewhat larger systems of ~1000 atoms
  so it may be better to get all the elements algorithmically, although it can be slower

  The commented members below is what we would want ideally, but in practice we have to find better way
  matrix* Ham;               // Electronic energy levels and couplings
  vector<matrix* > dHamdRx;  // derivatives of Ham wrt X coordinate of all atoms
  vector<matrix* > dHamdRy;  // derivatives of Ham wrt Y coordinate of all atoms
  vector<matrix* > dHamdRz;  // derivatives of Ham wrt Z coordinate of all atoms

  To solve this problem we use virtual functions for computation of H and dHdR 
  specific implementation is controlled by the derived classes

  Still - these data members are used in derived classes, so we include them here
  in efficient classes, they will simply not be used
*/

//protected:
public:

  int rep;                   ///< representation = 0 - for diabatic, 1 - for adiabatic
  int nelec;                 ///< number of electronic degrees of freedom (energy levels)
  int nnucl;                 ///< number of nuclear degrees of freedom - expected

  // Model-specific parameters
  vector<double> params;     ///< double-valued parameters of the Hamiltonian (the meaning is specific to each Hamiltonian type)
  vector<double> q;          ///< nuclear coordinates - here, they act as parameters
  vector<double> v;          ///< nuclear velocities: v = dq/dt  - here, they act as parameters

  // Diabatic representation
  MATRIX* ham_dia;           ///< Hamiltonian in diabatic representation
  vector<MATRIX*> d1ham_dia; ///< derivatives of the ham_dia w.r.t. all atomic DOFs: q0, q1, .. qN 
  vector<MATRIX*> d2ham_dia; ///< derivatives of the ham_dia w.r.t. all atomic DOFs: q00, q01, ..., q0N, q10, q11, ... qNN

  // Adiabatic representation
  MATRIX* ham_adi;           ///< Hamiltonian in adiabatic representation
  vector<MATRIX*> d1ham_adi; ///< first order derivative couplings: <i|d/dR|j> - is computed from the transformation coefficients

  int status_dia;    ///< Control variable of the diabatic calculations - to keep track of the computational state of the
                     ///< object of this calss - needed to avoid unnecessary computations
                     ///< if 0 - computations are outdated; 1 - computations are up to date
  int status_adi;    ///< Control variable of the adiabatic calculations - to keep track of the computational state of the
                     ///< object of this calss - needed to avoid unnecessary computations                                
                     ///< if 0 - computations are outdated; 1 - computations are up to date                              

//public:

  // Constructors
  Hamiltonian(); 

  // Use default copy constructor
  //Hamiltonian(const Hamiltonian&);

  // Destructor
  virtual ~Hamiltonian();    ///< This destructor does nothing - but in fact this is intended so the el and mol objects are not 
                             ///< destroyed, and the references are copied verbatim (not their content) - this helps prevent 
                             ///< additional overhead due to construction of el and mol objects

//  void set_status(int st_){ status = st_; }
//  int get_status(){ return status; }


  /** NOTE!!! It is very important to provide a minimal implementation of virtual methods
  // in the source file (Hamiltonian.cpp)
  // Otherwise the symbols remain undefined in the resulting library
  // This is a nasty bug to be aware of. 
  // See for example: http://stackoverflow.com/questions/1458180/vtable-for-referenced-from-compile-error-xcode/1478553#1478553
  // Updates coordinates, velocities, etc. - all nuclear information
  */
//  virtual void update_nuclear(Nuclear* mol);

//  virtual void set_q(vector<double>&){ ;; }


  virtual void set_rep(int rep_);

  // Set parameters
  virtual void set_params(vector<double>& params_){ ;; }
  virtual void set_params(boost::python::list params_);
  virtual void set_q(vector<double>& q_);
  virtual void set_q(boost::python::list q_);
  virtual void set_v(vector<double>& v_);
  virtual void set_v(boost::python::list v_);

  // Access to actual data
  virtual MATRIX get_ham_dia(){ return *ham_dia; }
  virtual MATRIX get_ham_adi(){ return *ham_adi; }
  virtual MATRIX get_d1ham_dia(int i){ return *d1ham_dia[i]; }
  virtual MATRIX get_d1ham_adi(int i){ return *d1ham_adi[i]; }


  // This function performs actual computations
  virtual void compute();
  virtual void compute_diabatic(){ ;; }  
  virtual void compute_adiabatic(){ ;; }  


  // Calculation methods
  virtual std::complex<double> H(int, int);             ///< Hamiltonian
  virtual std::complex<double> dHdq(int i,int j,int n); ///< Hamiltonian first-order derivative  
  virtual std::complex<double> D(int i,int j,int n);    ///< derivative coupling                 <i|d/dR_n|j>
  virtual std::complex<double> nac(int i,int j);        ///< non-adiabatic coupling              <i|d/dt|j>
  virtual std::complex<double> Hvib(int i,int j);       ///< vibronic Hamiltonian (for TD-SE)    H - i*hbar*nac


  friend bool operator == (const Hamiltonian& h1, const Hamiltonian& h2){
    return &h1 == &h2;
  }
  friend bool operator != (const Hamiltonian& h1, const Hamiltonian& h2){
    return !(h1 == h2);  // only compare addresses
  }



};

typedef std::vector<Hamiltonian> HamiltonianList;  ///< data type for keeping a list of generic Hamiltonians of their derived classes


}// namespace libhamiltonian_generic
}// namespace libatomistic
}// liblibra

#endif // HAMILTONIAN_H
