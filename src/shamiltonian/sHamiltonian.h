/*********************************************************************************
* Copyright (C) 2025 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file sHamiltonian.h
  \brief The file describes the generic swarm Hamiltonian class
    
*/

#ifndef sHAMILTONIAN_H
#define sHAMILTONIAN_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <complex>
#endif 


#include <torch/script.h>
#include <boost/python.hpp>
#include <torch/torch.h>
#include <torch/csrc/autograd/python_variable.h>

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{


/// libshamiltonian namespace
namespace libshamiltonian{


using namespace liblinalg;



bp::object tensor_to_python(const torch::Tensor& tensor);
torch::Tensor python_to_tensor(const bp::object& obj);


class sHamiltonian{
/**
  Keep the pointers to the objects (memory) or allocate the memory on demand

  For each variable X (except for the very trivial cases), we'll have a variable "X_mem_status", which 
  which describes the status of the memory for that variable:

  0 - not allocated,
  1 - allocated internally
  2 - allocated externally

*/


public:

  int nbeads;                ///< number of entangled trajectories (e.g. beads in RPMD or generally trajectories
                             /// that should be handled together)
  int nstates;               ///< number of electronic states
  int nnucl;                 ///< number of nuclear degrees of freedom - expected

  /** Diabatic representation

  |Psi> = |psi_dia> C_dia

  |psi_dia> = (psi_dia_0, psi_dia_1, ... psi_dia_{ndia-1})

  ovlp_dia = <psi_dia|psi_dia>

  dc1_dia[n] = <psi_dia| d/dq_n |psi_dia>

  ham_dia = <psi_dia|H|psi_dia>     - non-diagonal matrix

  nac_dia = <psi_dia|d/dt|psi_dia>

  hvib_dia = ham_dia - ihbar * nac_dia

  d1ham_dia[n] = d/dq_n [<psi_dia|H|psi_dia>]

  d2ham_dia[n*nnucl+k] = d^2/(dq_n dq_k)[ <psi_dia|H|psi_dia>]

  density_matrix_operator = |psi_dia> den_mat_dia <psi_dia| 

  So: den_mat_dia - density matrix matrix in the diabatic basis

  */



  // Diabatic properties
  torch::Tensor ovlp_dia;    ///< the overlap of the diabatic basis  [nbeads, nstates, nstates]
  torch::Tensor ham_dia;     ///< Hamiltonian in diabatic representation [nbeads, nstates, nstates]
  torch::Tensor nac_dia;     ///< Nonadiabatic couplings (time-derivative couplings) in diabatic representation [nbeads, nstates, nstates]
  torch::Tensor hvib_dia;    ///< Vibronic Hamiltonian in diabatic representation [nbeads, nstates, nstates]
  torch::Tensor dc1_dia;     ///< first-order derivative coupling matrices in the diabatic basis [nbeads, nnucl, nstates, nstates]
  torch::Tensor d1ham_dia;   ///< first order derivatives of the Hamiltonian matrix in the diabatic basis [nbeads, nnucl, nstates, nstates]
  torch::Tensor d2ham_dia;   ///< second order derivatives of the Hamiltonian in the diabatic basis [nbeads, nnucl, nnucl, nstates, nstates]



//  CMATRIX* den_mat_dia;       ///< Density matrix in the diabatic basis
//  int den_mat_dia_mem_status;


  /** Adiabatic representation: by definition, the states are orthonormal in this representation

  |Psi> = |psi_adi> C_adi

  |psi_adi> = (psi_adi_0, psi_adi_1, ... psi_adi_{nadi-1})

  dc1_adi[n] = <psi_adi| d/dq_n |psi_adi>

  ham_adi = <psi_adi|H|psi_adi>   - diagonal matrix

  nac_adi = <psi_adi|d/dt|psi_adi>

  hvib_adi = ham_adi - ihbar * nac_adi

  d1ham_adi[n] = d/dq_n [<psi_adi|H|psi_adi>]

  d2ham_adi[n*nnucl+k] =d^2/(dq_n dq_k)[ <psia_adi|H|psi_adi> ]

  density_matrix_operator = |psi_adi> den_mat_adi <psi_adi| 
  So: den_mat_adi - density matrix matrix in the adiabatic basis


*/


  torch::Tensor ovlp_adi;    ///< the overlap of the adiabatic basis  [nbeads, nstates, nstates]
  torch::Tensor ham_adi;     ///< Hamiltonian in adiabatic representation [nbeads, nstates, nstates]
  torch::Tensor nac_adi;     ///< Nonadiabatic couplings (time-derivative couplings) in adiabatic representation [nbeads, nstates, nstates]
  torch::Tensor hvib_adi;    ///< Vibronic Hamiltonian in adiabatic representation [nbeads, nstates, nstates]
  torch::Tensor dc1_adi;     ///< first-order derivative coupling matrices in the adiabatic basis [nbeads, nnucl, nstates, nstates]
  torch::Tensor d1ham_adi;   ///< first order derivatives of the Hamiltonian matrix in the adiabatic basis [nbeads, nnucl, nstates, nstates]
  torch::Tensor d2ham_adi;   ///< second order derivatives of the Hamiltonian in the adiabatic basis [nbeads, nnucl, nnucl, nstates, nstates]





  /** Transformations
  
  |psi_adi> = |psi_dia> * U

  Stationary SE:  H |psi_adi> =  |psi_adi> * E

  project on <psi_adi|, keeping in mind that <psi_adi|psi_adi> = I

  <psi_dia| H |psi_adi> =  <psi_dia|psi_adi> E  

  <psi_dia| H |psi_dia> * U =  <psi_dia|psi_dia> * U * E  

  H_dia * U = ovlp_dia * U * H_adi

  Here, U = basis_transform

  Conversion of the density matrices:

  |psi_dia> * U * den_mat_adi * U.H() * <psi_dia| = |psi_dia> * den_mat_dia * <psi_dia|

  so:    U * den_mat_adi * U.H() = den_mat_dia
    
  */


  torch::Tensor basis_transform;  ///< These are the eigenvectors of the generalized eigenvalue problem for diabatic Hamiltonian 
                                  ///< [nbeads, nstates, nstates]
  torch::Tensor time_overlap_adi;  ///< These are <psi_i_adi (t) | psi_j_adi (t+dt)> matrices  [nbeads, nstates, nstates]
  torch::Tensor time_overlap_dia;  ///< These are <psi_i_dia (t) | psi_j_dia (t+dt)> matrices  [nbeads, nstates, nstates]
  

  /**
     All the basic methods: constructor, destructor, getters, setters, etc.
  **/
  ///< In sHamiltonian_basic.cpp


  ///< Constructors
//  sHamiltonian(); 
  sHamiltonian(int nbeads_, int nstates_, int nnucl_);

  ///< Copy Constructor:
  sHamiltonian(const sHamiltonian&); 
 
  ///< Destructor
  ~sHamiltonian();    

  ///< Tree management utilities



  void bind(std::string name, torch::Tensor x);

  void compute(std::string property, bp::object py_funct, torch::Tensor q, bp::object params);
  void compute(bp::object py_funct, torch::Tensor q, bp::object params);

  void dia2adi();
  void compute_nacs_and_grads();

  torch::Tensor forces_adi(); // [nbeads, nstates, nnucl]
  torch::Tensor forces_dia(); // same
  
  

  ///< Setters:
/*
  friend bool operator == (const sHamiltonian& h1, const sHamiltonian& h2){
    return &h1 == &h2;
  }
  friend bool operator != (const sHamiltonian& h1, const sHamiltonian& h2){
    return !(h1 == h2);  // only compare addresses
  }
*/



};





}// namespace libshamiltonian
}// liblibra

#endif // sHAMILTONIAN_H
