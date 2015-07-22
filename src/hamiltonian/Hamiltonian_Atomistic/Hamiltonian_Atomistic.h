/*********************************************************************************
* Copyright (C) 2014 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef HAMILTONIAN_ATOMISTIC_H
#define HAMILTONIAN_ATOMISTIC_H

#include "../Hamiltonian.h"
//#include "Hamiltonian_MM/Hamiltonian_MM.h"


namespace libhamiltonian{
namespace libhamiltonian_atomistic{

using namespace libmmath;

class Hamiltonian_Atomistic : public Hamiltonian{

  // General specification of the model
  int rep;                   // representation = 0 - for diabatic, 1 - for adiabatic
  int nelec;                 // number of electronic degrees of freedom (energy levels)
  int nnucl;                 // number of nuclear degrees of freedom - expected

  // Model-specific parameters
  vector<double> params;
  vector<double> q;          // nuclear coordinates - here, they act as parameters
  vector<double> v;          // nuclear velocities: v = dq/dt  - here, they act as parameters

  // Diabatic representation
  MATRIX* ham_dia;           // Hamiltonian in diabatic representation
  vector<MATRIX*> d1ham_dia; // derivatives of the ham_dia w.r.t. all atomic DOFs: q0, q1, .. qN 
  vector<MATRIX*> d2ham_dia; // derivatives of the ham_dia w.r.t. all atomic DOFs: q00, q01, ..., q0N, q10, q11, ... qNN

  // Adiabatic representation
  MATRIX* ham_adi;           // Hamiltonian in adiabatic representation
  vector<MATRIX*> d1ham_adi; // first order derivative couplings: <i|d/dR|j> - is computed from the transformation coefficients

  // Control variable - to keep track of the computational state of the
  // object of this calss - needed to avoid unnecessary computations
  // if 0 - computations are outdated; 1 - computations are up to date
  int status_dia;  
  int status_adi;
 
  
public:

  // Constructor: only allocates memory and sets up related variables
  Hamiltonian_Atomistic(int, int);

  // Destructor
  ~Hamiltonian_Atomistic();

  // Set properties
  void set_rep(int rep_);

  // Set parameters
  void set_params(vector<double>& params_);
  void set_params(boost::python::list params_);
  void set_q(vector<double>& q_);
  void set_q(boost::python::list q_);
  void set_v(vector<double>& v_);
  void set_v(boost::python::list v_);

  // Perform actual computations - this will construct the internals of the object of this type
  void compute();
  void compute_diabatic();
  void compute_adiabatic();


  // Now call different properties - the call signature is the same, but the result depends in the setters before
  std::complex<double> H(int i,int j);          // Hamiltonian
  std::complex<double> dHdq(int i,int j,int n); // Hamiltonian first-order derivative  
  std::complex<double> D(int i,int j,int n);    // derivative coupling                 <i|d/dR_n|j>
  std::complex<double> nac(int i,int j);        // non-adiabatic coupling              <i|d/dt|j>
  std::complex<double> Hvib(int i,int j);       // vibronic Hamiltonian (for TD-SE)    H - i*hbar*nac


};

typedef std::vector<Hamiltonian_Atomistic> Hamiltonian_AtomisticList;


}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian

#endif // HAMILTONIAN_ATOMISTIC_H
