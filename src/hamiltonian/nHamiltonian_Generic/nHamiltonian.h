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
  \file nHamiltonian.h
  \brief The file describes the generic Hamiltonian class
    
*/

#ifndef nHAMILTONIAN_H
#define nHAMILTONIAN_H

#include <complex>
#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{


/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_generic namespace
namespace libhamiltonian_generic{


using namespace liblinalg;

void check_mat_dimensions(CMATRIX* x, int nrows, int ncols);
void check_vector_dimensions(vector<CMATRIX*> ptx, vector<CMATRIX>& x_, vector<int>& x_mem_status, int nnucl);


void set_X1_by_ref(CMATRIX* ptx, CMATRIX& x_, int& x_mem_status, int nrows, int ncols);
void set_X1_by_val(CMATRIX* ptx, CMATRIX& x_, int& x_mem_status, int nrows, int ncols);
void init_X1(CMATRIX* ptx, int& x_mem_status, int nrows, int ncols);

void set_X2_by_ref(vector<CMATRIX*> ptx, vector<CMATRIX>& x_, vector<int>& x_mem_status, int nrows, int ncols, int nnucl);
void set_X2_by_val(vector<CMATRIX*> ptx, vector<CMATRIX>& x_, vector<int>& x_mem_status, int nrows, int ncols, int nnucl);
void init_X2(vector<CMATRIX*> ptx, vector<int>& x_mem_status, int nrows, int ncols, int nnucl);






class nHamiltonian{
/**
  Keep the pointers to the objects (memory) or allocate the memory on demand

  For each variable X (except for the very trivial cases), we'll have a variable "X_mem_status", which 
  which describes the status of the memory for that variable:

  0 - not allocated,
  1 - allocated internally
  2 - allocated externally

*/



public:

  int nnucl;                 ///< number of nuclear degrees of freedom - expected

  /** Diabatic representation

  |Psi> = |psi_dia> C_dia

  |psi_dia> = (psi_dia_0, psi_dia_1, ... psi_dia_{ndia-1})

  ovlp_dia = <psi_dia|psi_dia>

  dc1_dia[n] = <psi_dia| d/dq_n |psi_dia>

  ham_dia = <psi_dia|H|psi_dia>     - non-diagonal matrix

  d1ham_dia[n] = d/dq_n [<psi_dia|H|psi_dia>]

  d2ham_dia[n*nnucl+k] = d^2/(dq_n dq_k)[ <psi_dia|H|psi_dia>]

  density_matrix_operator = |psi_dia> den_mat_dia <psi_dia| 

  So: den_mat_dia - density matrix matrix in the diabatic basis

  */
  int ndia;                   ///< the number of electronic DOFs in the diabatic basis

  CMATRIX* ampl_dia;          ///< Amplitudes of the diabatic states in the overal wavefunction
  int ampl_dia_mem_status;

  CMATRIX* ovlp_dia;          ///< the overlap of the diabatic basis 
  int ovlp_dia_mem_status;       

  vector<CMATRIX*> dc1_dia;   ///< first-order derivative coupling matrices in the diabatic basis 
  vector<int> dc1_dia_mem_status;


  CMATRIX* ham_dia;           ///< Hamiltonian in diabatic representation
  int ham_dia_mem_status;


  vector<CMATRIX*> d1ham_dia; ///< first order derivatives of the Hamiltonian matrix in the diabatic basis
  vector<int> d1ham_dia_mem_status;

  vector<CMATRIX*> d2ham_dia; ///< second order derivatives of the Hamiltonian in the diabatic basis
  vector<int> d2ham_dia_mem_status;

  CMATRIX* den_mat_dia;       ///< Density matrix in the diabatic basis
  int den_mat_dia_mem_status;


  /** Adiabatic representation: by definition, the states are orthonormal in this representation

  |Psi> = |psi_adi> C_adi

  |psi_adi> = (psi_adi_0, psi_adi_1, ... psi_adi_{nadi-1})

  dc1_adi[n] = <psi_adi| d/dq_n |psi_adi>

  ham_adi = <psi_adi|H|psi_adi>   - diagonal matrix

  d1ham_adi[n] = d/dq_n [<psi_adi|H|psi_adi>]

  d2ham_adi[n*nnucl+k] =d^2/(dq_n dq_k)[ <psia_adi|H|psi_adi> ]

  density_matrix_operator = |psi_adi> den_mat_adi <psi_adi| 
  So: den_mat_adi - density matrix matrix in the adiabatic basis

  */
  int nadi;                   ///< the number of electronic DOFs in the adiabatic basis

  CMATRIX* ampl_adi;          ///< Amplitudes of the adiabatic states in the overal wavefunction
  int ampl_adi_mem_status;

  vector<CMATRIX*> dc1_adi;   ///< first-order derivative coupling matrices in the adiabatic basis 
  vector<int> dc1_adi_mem_status;

  CMATRIX* ham_adi;           ///< Hamiltonian in adiabatic representation (diagonal)
  int ham_adi_mem_status;

  vector<CMATRIX*> d1ham_adi; ///< first order derivatives of the Hamiltonian matrix in the adiabatic basis (diagonal)
  vector<int> d1ham_adi_mem_status;

  vector<CMATRIX*> d2ham_adi; ///< second order derivatives of the Hamiltonian matrix in the adiabatic basis
  vector<int> d2ham_adi_mem_status;

  CMATRIX* den_mat_adi;       ///< Density matrix in the adiabatic basis
  int den_mat_adi_mem_status;


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

  CMATRIX* basis_transform;  ///< These are the eigenvectors of the generalized eigenvalue problem for diabatic Hamiltonian
  int basis_transform_mem_status;


  /**========= in Hamiltonian.cpp ==================
     
     All the basic methods: constructor, destructor, getters, setters, etc.

  **/
  ///< Constructors
  nHamiltonian(int ndia_, int nadi_, int nnucl_); 

  ///< Destructor
  ~nHamiltonian();    


  ///< Setters:
  // Diabatic
  void set_ampl_dia_by_ref(CMATRIX& ampl_dia_);
  void set_ampl_dia_by_val(CMATRIX& ampl_dia_);

  void set_ovlp_dia_by_ref(CMATRIX& ovlp_dia_);
  void set_ovlp_dia_by_val(CMATRIX& ovlp_dia_);

  void set_dc1_dia_by_ref(vector<CMATRIX>& dc1_dia_);
  void set_dc1_dia_by_val(vector<CMATRIX>& dc1_dia_);

  void set_ham_dia_by_ref(CMATRIX& ham_dia_);
  void set_ham_dia_by_val(CMATRIX& ham_dia_);

  void set_d1ham_dia_by_ref(vector<CMATRIX>& d1ham_dia_);
  void set_d1ham_dia_by_val(vector<CMATRIX>& d1ham_dia_);

  void set_d2ham_dia_by_ref(vector<CMATRIX>& d2ham_dia_);
  void set_d2ham_dia_by_val(vector<CMATRIX>& d2ham_dia_);

  void set_den_mat_dia_by_ref(CMATRIX& den_mat_dia_);
  void set_den_mat_dia_by_val(CMATRIX& den_mat_dia_);

  // Adiabatic
  void set_ampl_adi_by_ref(CMATRIX& ampl_adi_);
  void set_ampl_adi_by_val(CMATRIX& ampl_adi_);

  void set_dc1_adi_by_ref(vector<CMATRIX>& dc1_adi_);
  void set_dc1_adi_by_val(vector<CMATRIX>& dc1_adi_);

  void set_ham_adi_by_ref(CMATRIX& ham_adi_);
  void set_ham_adi_by_val(CMATRIX& ham_adi_);

  void set_d1ham_adi_by_ref(vector<CMATRIX>& d1ham_adi_);
  void set_d1ham_adi_by_val(vector<CMATRIX>& d1ham_adi_);

  void set_d2ham_adi_by_ref(vector<CMATRIX>& d2ham_adi_);
  void set_d2ham_adi_by_val(vector<CMATRIX>& d2ham_adi_);

  void set_den_mat_adi_by_ref(CMATRIX& den_mat_adi_);
  void set_den_mat_adi_by_val(CMATRIX& den_mat_adi_);

  // Transforms
  void set_basis_transform_by_ref(CMATRIX& basis_transform_);
  void set_basis_transform_by_val(CMATRIX& basis_transform_);


  /// Getters
  // Diabatic
  CMATRIX get_ovlp_dia();

  CMATRIX get_ampl_dia();
  CMATRIX get_dc1_dia(int i);
  CMATRIX get_ham_dia();
  CMATRIX get_d1ham_dia(int i);
  CMATRIX get_d2ham_dia(int i);
  CMATRIX get_d2ham_dia(int i,int j);
  CMATRIX get_den_mat_dia();

  // Adiabatic
  CMATRIX get_ampl_adi();
  CMATRIX get_dc1_adi(int i);
  CMATRIX get_ham_adi();
  CMATRIX get_d1ham_adi(int i);
  CMATRIX get_d2ham_adi(int i);
  CMATRIX get_d2ham_adi(int i,int j);
  CMATRIX get_den_mat_adi();

  // Transforms
  CMATRIX get_basis_transform();


  /**========= in Hamiltonian1.cpp ==================

             Computational methods

  **/

  void compute_adiabatic(int lvl);
  void ampl_dia2adi();
  void ampl_adi2dia();


  complex<double> Ehrenfest_energy_adi();
  complex<double> Ehrenfest_energy_dia();


  CMATRIX forces_adi();  // -dE/dR in the adiabatic basis, assuming Cadi = Cadi(t)
  CMATRIX forces_dia();  // -dE/dR in the diabatic basis, assuming Cdia = Cdia(t)

  vector<CMATRIX> forces_tens_adi(); // 
  vector<CMATRIX> forces_tens_dia(); // 


//  CMATRIX forces_adi();  // -dE/dR in the adiabatic basis, assuming Cadi = Cadi(t)
//  CMATRIX forces_dia();  // -dE/dR in the diabatic basis, assuming Cdia = Cdia(t)


  CMATRIX Ehrenfest_forces_adi();  // Ehrenfest forces in adiabatic basis
  CMATRIX Ehrenfest_forces_dia();  // Ehrenfest forces in diabatic basis

  vector<CMATRIX> Ehrenfest_forces_tens_adi();  // Force tensor in adiabatic basis, assuming Cadi = Cadi(t)
  vector<CMATRIX> Ehrenfest_forces_tens_dia();  // Force tensor in diabatic basis, assuming Cdia = Cdia(t)



/*

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
*/

  friend bool operator == (const nHamiltonian& h1, const nHamiltonian& h2){
    return &h1 == &h2;
  }
  friend bool operator != (const nHamiltonian& h1, const nHamiltonian& h2){
    return !(h1 == h2);  // only compare addresses
  }



};

typedef std::vector<nHamiltonian> nHamiltonianList;  ///< data type for keeping a list of generic Hamiltonians of their derived classes


}// namespace libhamiltonian_generic
}// namespace libhamiltonian
}// liblibra

#endif // nHAMILTONIAN_H
