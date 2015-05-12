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

#ifndef HAMILTONIAN_MODEL_H
#define HAMILTONIAN_MODEL_H

#include "Hamiltonian.h"
#include "Model_SAC.h"
#include "Model_DAC.h"
#include "Model_ECWR.h"
#include "Model_Marcus.h"
#include "Model_SEXCH.h"
using namespace libmmath;


namespace libhamiltonian{


class Hamiltonian_Model : public Hamiltonian{

  int ham_indx;              // model index: 0 - SAC, 1 - DAC, 2 - ECWR, 3 - Marcus, 4 - superexchange (SEXCH)

  int rep;                   // representation = 0 - for diabatic, 1 - for adiabatic

  int n_elec;                // number of electronic degrees of freedom (energy levels)
  int n_nucl;                // number of nuclear degrees of freedom - expected

  // Parameters
  vector<double> params;

  // Diabatic representation
  MATRIX* ham_dia;           // Hamiltonian in diabatic representation
  vector<MATRIX*> d1ham_dia; // derivatives of the ham_dia w.r.t. all atomic DOFs: q0, q1, .. qN 
  vector<MATRIX*> d2ham_dia; // derivatives of the ham_dia w.r.t. all atomic DOFs: q00, q01, ..., q0N, q10, q11, ... qNN

  // Adiabatic representation
  MATRIX* ham_adi;           // Hamiltonian in adiabatic representation
  vector<MATRIX*> d1ham_adi; // first order derivative couplings: <i|d/dR|j> - is computed from the transformation coefficients
 
  
public:

  // Constructor: only allocates memory and sets up related variables
  Hamiltonian_Model(int ham_indx_);

  // Destructor
  ~Hamiltonian_Model();

  // Set properties
  void set_rep(int rep_);

  void set_params(vector<double>& params_);
  void set_params(boost::python::list params_);

  // Perform actual computations - this will construct the internals of the object of this type
  void compute_diabatic(vector<double>& q_, vector<double>& p_);
  void compute_diabatic(boost::python::list q_, boost::python::list p_);

  void compute_adiabatic(vector<double>& q_, vector<double>& p_);
  void compute_adiabatic(boost::python::list q_, boost::python::list p_);


  // Now call different properties - the call signature is the same, but the result depends in the setters before
  std::complex<double> H(int i,int j);

/*
  void set_model(int pot_indx_);
  void set_param(std::string var,double val);

//  void update_nuclear(Nuclear* mol);
//  void set_position(double x_);
//  void set_velocity(double v_);



  // Computation: first call one of these functions - to compute Hamiltonian for given representation
  void compute(Nuclear* mol);
  void compute_diabatic(Nuclear* mol);
  void compute_adiabatic(Nuclear* mol);


  // Return results: this does not do computations - only return precomputed results
  std::complex<double> H(Nuclear* mol,int i,int j);
  std::complex<double> H(Nuclear* mol,int i,int j, int rep_);
  std::complex<double> dHdRx(Nuclear* mol,int i,int j,int k);
  std::complex<double> dHdRy(Nuclear* mol,int i,int j,int k);
  std::complex<double> dHdRz(Nuclear* mol,int i,int j,int k);

  std::complex<double> Dx(Nuclear* mol,int, int, int); // derivative coupling w.r.t. Cartesian coordinate x
  std::complex<double> Dy(Nuclear* mol,int, int, int); // derivative coupling w.r.t. Cartesian coordinate y
  std::complex<double> Dz(Nuclear* mol,int, int, int); // derivative coupling w.r.t. Cartesian coordinate z
  
*/


};



}// namespace libmodel

#endif // HAMILTONIAN_MODEL_H
