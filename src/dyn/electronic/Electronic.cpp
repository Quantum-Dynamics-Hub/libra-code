/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Electronic.cpp
  \brief The file implements Electronic class methods and functions for propagation of electronic DOF
    
*/

#include "Electronic.h"

/// liblibra namespace
namespace liblibra{


/// libdyn namespace 
namespace libdyn{

/// libelectronic namespace 
namespace libelectronic{

using namespace liblinalg;


Electronic& Electronic::operator=(const Electronic& ob){  
/**
  \brief An assignment operator

*/

  nstates = ob.nstates;
  istate = ob.istate;
  q = ob.q;
  p = ob.p;

  return *this; // return reference to allow chaining: A = B = C =...
}


void Electronic::rnd_phase(double& x, double& y, double nrm, double phi){ 
/**
  \brief An auxiliary function to generate a random (uniform distribution) phase complex number

  \param[out] x The real component of the complex number generated
  \param[out] y The imaginary component of the complex number generated
  \param[in] nrm The norm of the generated random number
  \param[in] phi The phase of the wavefunction

*/                                         

  x = std::sqrt(nrm) * std::cos(M_PI*phi); 
  y = std::sqrt(nrm) * std::sin(M_PI*phi); 

}

void Electronic::init(int n_,int st, double phi){
/**
  \brief An auxiliary function to initialize (or reinitialize the Electronic object)

  This function allocates memory for time-dependent wfc with n_ stationary states
  Then it initializes the overall multiconfigurational wfc
  to be a 1-configurational, with the weight 1 set to basis state with index st
  and with a random phase

  \param[in] n_ The number of electronic states to include: the dimensionality of the electronic problem
  \param[in] st The index of electronic state to which we initialize the system 
  \param[in] phi The phase of the wavefunction

*/
//  rnd_obj = new Random();

  if(st>=n_){ std::cout<<"Error in Electronic::init - st("<<st<<") must be smaller than n_("<<n_<<")\n"; exit(0); }

  nstates = n_;
  q = std::vector<double>(n_,0.0);
  p = std::vector<double>(n_,0.0); 

  istate = st;
  rnd_phase(q[istate],p[istate],1.0, phi);    // populate only the istate-th state

}

//
// Overloaded version
//
void Electronic::init(int n_, int st){ 
/**
  \brief An auxiliary function to initialize (or reinitialize the Electronic object)

  This function allocates memory for time-dependent wfc with n_ stationary states
  Then it initializes the overall multiconfigurational wfc
  to be a 1-configurational, with the weight 1 set to basis state with index 0
  and with a random phase

  \param[in] n_ The number of electronic states to include: the dimensionality of the electronic problem
  \param[in] st The index of electronic state to which we initialize the system 

*/

  init(n_,st, 0.0);
}  


//
// Overloaded version
//
void Electronic::init(int n_){ 
/**
  \brief An auxiliary function to initialize (or reinitialize the Electronic object)

  This function allocates memory for time-dependent wfc with n_ stationary states
  Then it initializes the overall multiconfigurational wfc
  to be a 1-configurational, with the weight 1 set to basis state with index 0
  and with a random phase

  \param[in] n_ The number of electronic states to include: the dimensionality of the electronic problem

*/

  init(n_,0, 0.0);
}  


Electronic::Electronic(int n_,int st, double phi){ 
/**
  \brief Constructor

  \param[in] n_ The number of electronic states to include: the dimensionality of the electronic problem
  \param[in] st The index of electronic state to which we initialize the system 
  \param[in] phi The phase of the wavefunction
*/

 init(n_,st, phi);
}


//
// Constructors
//
Electronic::Electronic(int n_,int st){ 
/**
  \brief Constructor

  \param[in] n_ The number of electronic states to include: the dimensionality of the electronic problem
  \param[in] st The index of electronic state to which we initialize the system 
*/

 init(n_,st, 0.0);
}

Electronic::Electronic(int n_){ 
/**
  \brief Constructor

  The system is initialized to be in the lowest (index 0) electronic state.
  \param[in] n_ The number of electronic states to include: the dimensionality of the electronic problem
*/

 init(n_,0, 0.0); 
}

Electronic::Electronic(){ 
/**
  \brief Constructor

  The dimensionality of the electronic problem is set to 1 (no excited states)
  The system is initialized to be in the lowest (index 0) electronic state - the only one available
*/

 init(1,0, 0.0);
}


Electronic::Electronic(const Electronic& ob){ /// cctor
/**
  \brief Copy constructor
*/

  nstates = ob.nstates;
  istate = ob.istate;
  q = ob.q;
  p = ob.p;

}


Electronic::~Electronic(){  
/**
  \brief Destructor
*/
  if(q.size()>0){ q.clear(); }
  if(p.size()>0){ p.clear(); }
}



std::complex<double> Electronic::c(int i) const{
/**
  \brief Return the amplitude of a quantum state in the complex format: c_i = q_i + i*p_i

  \param[in] i Index of the quantum state
*/

  return complex<double>(q[i],p[i]);

}

std::complex<double> Electronic::rho(int i, int j) const{
/**
  \brief Returns the density matrix element: rho_ij = c^*_i * c_j

  \param[in] i index of one state
  \param[in] j index of another state
*/

  return complex<double>((q[i]*q[j] + p[i]*p[j]), (q[i]*p[j]-p[i]*q[j]));

}

CMATRIX Electronic::C() const{
/**
  \brief Return the amplitudes of all quantum states in the vector format:
        ( c_0 )
   C =  ( ... )
        ( c_N )
   with  c_i = q_i + i*p_i

*/
  CMATRIX res(nstates,1);
  for(int i=0; i<nstates; i++){  res.set(i, q[i], p[i]);  }
  return res;
}


CMATRIX Electronic::RHO() const{
/**
  \brief Return the density matrix computed for all quantum states:

   P = C * C^+

  The density matrix elements are: P_ij = rho_ij = c^*_i * c_j

*/
  CMATRIX res(nstates,nstates);
  for(int i=0; i<nstates; i++){
    for(int j=0; j<nstates; j++){

      res.set(i, j, (q[i]*q[j] + p[i]*p[j]),  (q[i]*p[j]-p[i]*q[j]) );  
    }
  }
  return res;
}




}//namespace libelectronic
}// namespace libdyn
}// liblibra
