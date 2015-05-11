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
using namespace libmmath;


namespace libhamiltonian{


class Hamiltonian_Model : public Hamiltonian{

  int ham_indx;              // model index: 0 - SAC, 1 - DAC, 2 - ECWR, 3 - Marcus, 4 - superexchange

  int n_elec;                // number of electronic degrees of freedom (energy levels)
  int n_nucl;                // number of nuclear degrees of freedom - expected

  // Diabatic representation
  MATRIX* ham_dia;           // Hamiltonian in diabatic representation
  vector<MATRIX*> d1ham_dia; // derivatives of the ham_dia w.r.t. all atomic DOFs: q0, q1, .. qN 
  vector<MATRIX*> d2ham_dia; // derivatives of the ham_dia w.r.t. all atomic DOFs: q00, q01, ..., q0N, q10, q11, ... qNN

  // Adiabatic representation
  MATRIX* ham_adi;           // Hamiltonian in adiabatic representation

 

/*
  // Status

  // Input parameters
  int rep;             // representation = 0 - for diabatic, 1 - for adiabatic

  double x;            // position
  double v;            // velocity

  // Model parameters
  double A,B,C,D,E;
  double A1,B1,C1,D1, A2,B2,C2,D2, A3,B3,C3,D3;

  // Computed results
  // Diabatic 
  double H00,H11,H22,H01,H02,H12,   dH00,dH11,dH22,dH01,dH02,dH12,  d2H00,d2H11,d2H22,d2H01,d2H02,d2H12;
         

  // Adiabatic
  double E0,E1,dE0,dE1,d01,d10,dd01,dd10,D01,D10;

*/

  
public:

  // Constructor: only allocates memory and sets up related variables
  Hamiltonian_Model(int ham_indx_);

 
  void compute_diabatic(vector<double>&, vector<double>&);
//  Hamiltonian_Model(int rep_, int pot_indx_, Nuclear* mol);

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
