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
  \file Hamiltonian_Model.h
  \brief The file describes the model Hamiltonian class
    
*/

#ifndef HAMILTONIAN_MODEL_H
#define HAMILTONIAN_MODEL_H

#include "../Hamiltonian_Generic/Hamiltonian.h"
#include "Model_SAC.h"
#include "Model_DAC.h"
#include "Model_ECWR.h"
#include "Model_Marcus.h"
#include "Model_SEXCH.h"
#include "Model_Rabi2.h"
#include "Model_sin.h"
#include "Model_sin_2D.h"


/// libhamiltonian namespace
namespace libhamiltonian{

using namespace libhamiltonian_generic;
using namespace libmmath;

/// libhamiltonian_model namespace
namespace libhamiltonian_model{


class Hamiltonian_Model : public Hamiltonian{
/**
  This is a class derived from the generic Hamiltonian class, so it inherits many of its properties
*/

  /// General specification of the model: 0 - SAC, 1 - DAC, 2 - ECWR, 3 - Marcus, 4 - superexchange (SEXCH), 5 - Rabi2
  int ham_indx;    
   
public:

  /// Constructor: only allocates memory and sets up related variables
  Hamiltonian_Model(int ham_indx_);

  /// Destructor
  ~Hamiltonian_Model();

  // Set properties
//  void set_rep(int rep_);

  // Set parameters
  void set_params(vector<double>& params_);
//  void set_params(boost::python::list params_);
//  void set_q(vector<double>& q_);
//  void set_q(boost::python::list q_);
//  void set_v(vector<double>& v_);
//  void set_v(boost::python::list v_);

  // Perform actual computations - this will construct the internals of the object of this type
//  void compute();
  void compute_diabatic();
  void compute_adiabatic();


  // Now call different properties - the call signature is the same, but the result depends in the setters before
//  std::complex<double> H(int i,int j);          // Hamiltonian
//  std::complex<double> dHdq(int i,int j,int n); // Hamiltonian first-order derivative  
//  std::complex<double> D(int i,int j,int n);    // derivative coupling                 <i|d/dR_n|j>
//  std::complex<double> nac(int i,int j);        // non-adiabatic coupling              <i|d/dt|j>
//  std::complex<double> Hvib(int i,int j);       // vibronic Hamiltonian (for TD-SE)    H - i*hbar*nac


};

typedef std::vector<Hamiltonian_Model> Hamiltonian_ModelList;  /// data type for keeping a list of model Hamiltonians of their derived classes


}// namespace libhamiltonian_model
}// namespace libmodel

#endif // HAMILTONIAN_MODEL_H
