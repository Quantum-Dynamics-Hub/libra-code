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
  \file nHamiltonian.h
  \brief The file describes the generic Hamiltonian class
    
*/

#ifndef nHAMILTONIAN_H
#define nHAMILTONIAN_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <complex>
#endif 

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{


/// libnhamiltonian namespace
namespace libnhamiltonian{


using namespace liblinalg;

void check_mat_dimensions(CMATRIX* x, int nrows, int ncols);
void check_vector_dimensions(vector<CMATRIX*> ptx, vector<CMATRIX>& x_, vector<int>& x_mem_status, int nnucl);


void set_X1_by_ref(CMATRIX* ptx, CMATRIX& x_, int& x_mem_status, int nrows, int ncols);
void set_X1_by_val(CMATRIX* ptx, CMATRIX& x_, int& x_mem_status, int nrows, int ncols);
void init_X1(CMATRIX* ptx, int& x_mem_status, int nrows, int ncols);

void set_X2_by_ref(vector<CMATRIX*> ptx, vector<CMATRIX>& x_, vector<int>& x_mem_status, int nrows, int ncols, int nnucl);
void set_X2_by_val(vector<CMATRIX*> ptx, vector<CMATRIX>& x_, vector<int>& x_mem_status, int nrows, int ncols, int nnucl);
void init_X2(vector<CMATRIX*> ptx, vector<int>& x_mem_status, int nrows, int ncols, int nnucl);



// Forward declaration
//class nHamiltonian; 


class nHamiltonian{
/**
  Keep the pointers to the objects (memory) or allocate the memory on demand

  For each variable X (except for the very trivial cases), we'll have a variable "X_mem_status", which 
  which describes the status of the memory for that variable:

  0 - not allocated,
  1 - allocated internally
  2 - allocated externally

*/


  void check_cmatrix(bp::object obj, std::string matrix_name, int nrows, int ncols);
  void check_cmatrix_list(bp::object obj, std::string matrix_name, int nrows, int ncols, int length);

  void add_branches(int target_level, vector<nHamiltonian*>& res);

public:

  int level;                        ///< level in the tree hierarchy
  int id;                           ///< index of this Hamiltonian in this level of hierarchy
  nHamiltonian* parent;             ///< the Hamiltonian of a higher level
  vector<nHamiltonian*> children;   ///< the Hamiltonians of the lower level


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
  int ndia;                   ///< the number of electronic DOFs in the diabatic basis

//  CMATRIX* ampl_dia;          ///< Amplitudes of the diabatic states in the overal wavefunction
//  int ampl_dia_mem_status;

  CMATRIX* ovlp_dia;          ///< the overlap of the diabatic basis 
  int ovlp_dia_mem_status;       

  vector<CMATRIX*> dc1_dia;   ///< first-order derivative coupling matrices in the diabatic basis 
  vector<int> dc1_dia_mem_status;


  CMATRIX* ham_dia;           ///< Hamiltonian in diabatic representation
  int ham_dia_mem_status;

  CMATRIX* nac_dia;           ///< Nonadiabatic couplings (time-derivative couplings) in diabatic representation
  int nac_dia_mem_status;

  CMATRIX* hvib_dia;          ///< Vibronic Hamiltonian in diabatic representation
  int hvib_dia_mem_status;




  vector<CMATRIX*> d1ham_dia; ///< first order derivatives of the Hamiltonian matrix in the diabatic basis
  vector<int> d1ham_dia_mem_status;

  vector<CMATRIX*> d2ham_dia; ///< second order derivatives of the Hamiltonian in the diabatic basis
  vector<int> d2ham_dia_mem_status;

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
  int nadi;                   ///< the number of electronic DOFs in the adiabatic basis


  vector<CMATRIX*> dc1_adi;   ///< first-order derivative coupling matrices in the adiabatic basis 
  vector<int> dc1_adi_mem_status;

  CMATRIX* ham_adi;           ///< Hamiltonian in adiabatic representation (diagonal)
  int ham_adi_mem_status;


  CMATRIX* nac_adi;           ///< Nonadiabatic couplings (time-derivative couplings) in adiabatic representation
  int nac_adi_mem_status;

  CMATRIX* hvib_adi;          ///< Vibronic Hamiltonian in adiabatic representation
  int hvib_adi_mem_status;


  vector<CMATRIX*> d1ham_adi; ///< first order derivatives of the Hamiltonian matrix in the adiabatic basis (diagonal)
  vector<int> d1ham_adi_mem_status;

  vector<CMATRIX*> d2ham_adi; ///< second order derivatives of the Hamiltonian matrix in the adiabatic basis
  vector<int> d2ham_adi_mem_status;



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

  CMATRIX* time_overlap_adi;  ///< These are <psi_i_adi (t) | psi_j_adi (t+dt)> matrices
  int time_overlap_adi_mem_status;

  CMATRIX* time_overlap_dia;  ///< These are <psi_i_dia (t) | psi_j_dia (t+dt)> matrices
  int time_overlap_dia_mem_status;


  vector<int>* ordering_adi;  ///< the permutation describing how the indices of the internally-stored
                              /// properties like ham_adi or nac_adi relate to the indices of the 
                              /// actual "physical" states.

  CMATRIX* cum_phase_corr;    /// These are the cumulative phase corrections applied to basis_transform (this correction is indexed by
                              /// the inernal indexing of states, not physical)
                              /// Keep in mind that these phases have already been applied to basis_transform variable
                              /// (this is the convention)
  int cum_phase_corr_mem_status;


  /**     
     All the computation control parameters 
  **/
  int eigen_algo;  ///< eigensolver algorithm used to convert diabatic to adiabatic basis 
                   /// 0 - generalized eigensolver (may reorder states) [default]
                   /// 1 - eigensolver without reordering (assumes diabatic basis is orthonormaly)

  double phase_corr_ovlp_tol;  /// phase correction overlap tolerance, we only compute phase corrections 
                               /// if the time overlaps are above (in magnitude) this value, otherwise
                               /// we assume that uncommon state reordering (adiabatic state switching) 
                               /// has happened (which may occur in the stochastic state tracking method)
                               /// then, we don't apply the phase correction since it doesn't make sense



  /**
     
     All the basic methods: constructor, destructor, getters, setters, etc.

  **/
  ///< In nHamiltonian_basic.cpp

  void init_mem_status(int ndia_, int nadi_, int nnucl_);

  ///< Constructors
//  nHamiltonian(); 
  nHamiltonian(int ndia_, int nadi_, int nnucl_); 

  ///< Copy Constructor:
  nHamiltonian(const nHamiltonian&); 
 
  ///< Destructor
  ~nHamiltonian();    

  ///< Tree management utilities
  void set_levels(int lvl_);            ///< Updates node levels downstream
  void add_new_children(int _ndia, int _nadi, int _nnucl, int nchildren);
  void add_child(nHamiltonian& child);  ///< Associate an existing Hamiltonian with the present one
                                        ///< to become its child
  vector<int> get_full_id();            ///< Entire path of the note in the tree
  vector<nHamiltonian*> get_branches(int target_level);


  void copy_content(const nHamiltonian& src);
  void copy_content(nHamiltonian* src);
  void copy_level_content(nHamiltonian* src); ///< To copy content if all the memory is allocated

  void init_all(int der_lvl);
  void init_all(int der_lvl, int lvl);

  void show_memory_status(vector<int>& id_);


  ///< Setters:
  // Diabatic
  void init_ovlp_dia();
  void set_ovlp_dia_by_ref(CMATRIX& ovlp_dia_);
  void set_ovlp_dia_by_val(CMATRIX& ovlp_dia_);

  void init_dc1_dia();
  void set_dc1_dia_by_ref(vector<CMATRIX>& dc1_dia_);
  void set_dc1_dia_by_val(vector<CMATRIX>& dc1_dia_);

  void init_ham_dia();
  void set_ham_dia_by_ref(CMATRIX& ham_dia_);
  void set_ham_dia_by_val(CMATRIX& ham_dia_);

  void init_nac_dia();
  void set_nac_dia_by_ref(CMATRIX& nac_dia_);
  void set_nac_dia_by_val(CMATRIX& nac_dia_);

  void init_hvib_dia();
  void set_hvib_dia_by_ref(CMATRIX& hvib_dia_);
  void set_hvib_dia_by_val(CMATRIX& hvib_dia_);

  void init_d1ham_dia();
  void set_d1ham_dia_by_ref(vector<CMATRIX>& d1ham_dia_);
  void set_d1ham_dia_by_val(vector<CMATRIX>& d1ham_dia_);

  void init_d2ham_dia();
  void set_d2ham_dia_by_ref(vector<CMATRIX>& d2ham_dia_);
  void set_d2ham_dia_by_val(vector<CMATRIX>& d2ham_dia_);




  // Adiabatic
  void init_dc1_adi();
  void set_dc1_adi_by_ref(vector<CMATRIX>& dc1_adi_);
  void set_dc1_adi_by_val(vector<CMATRIX>& dc1_adi_);

  void init_ham_adi();
  void set_ham_adi_by_ref(CMATRIX& ham_adi_);
  void set_ham_adi_by_val(CMATRIX& ham_adi_);

  void init_nac_adi();
  void set_nac_adi_by_ref(CMATRIX& nac_adi_);
  void set_nac_adi_by_val(CMATRIX& nac_adi_);

  void init_hvib_adi();
  void set_hvib_adi_by_ref(CMATRIX& hvib_adi_);
  void set_hvib_adi_by_val(CMATRIX& hvib_adi_);

  void init_d1ham_adi();
  void set_d1ham_adi_by_ref(vector<CMATRIX>& d1ham_adi_);
  void set_d1ham_adi_by_val(vector<CMATRIX>& d1ham_adi_);

  void init_d2ham_adi();
  void set_d2ham_adi_by_ref(vector<CMATRIX>& d2ham_adi_);
  void set_d2ham_adi_by_val(vector<CMATRIX>& d2ham_adi_);

  // Transforms
  void init_basis_transform();
  void set_basis_transform_by_ref(CMATRIX& basis_transform_);
  void set_basis_transform_by_val(CMATRIX& basis_transform_);

  void init_time_overlap_adi();
  void set_time_overlap_adi_by_ref(CMATRIX& time_overlap_adi_);
  void set_time_overlap_adi_by_val(CMATRIX& time_overlap_adi_);

  void init_time_overlap_dia();
  void set_time_overlap_dia_by_ref(CMATRIX& time_overlap_dia_);
  void set_time_overlap_dia_by_val(CMATRIX& time_overlap_dia_);


  void set_ordering_adi_by_ref(vector<int>& ordering_adi_);
  void set_ordering_adi_by_val(vector<int>& ordering_adi_);

  void init_cum_phase_corr();
  void set_cum_phase_corr_by_ref(CMATRIX& cum_phase_corr_);
  void set_cum_phase_corr_by_val(CMATRIX& cum_phase_corr_);


  /// Getters
  // Diabatic
  CMATRIX get_ovlp_dia();
  CMATRIX get_ovlp_dia(vector<int>& id_);
  CMATRIX get_dc1_dia(int i);
  CMATRIX get_dc1_dia(int i, vector<int>& id_);
  CMATRIX get_ham_dia();
  CMATRIX get_ham_dia(vector<int>& id_);
  CMATRIX get_nac_dia();
  CMATRIX get_nac_dia(vector<int>& id_);
  CMATRIX get_hvib_dia();
  CMATRIX get_hvib_dia(vector<int>& id_);
  CMATRIX get_d1ham_dia(int i);
  CMATRIX get_d1ham_dia(int i, vector<int>& id_);
  CMATRIX get_d2ham_dia(int i);
  CMATRIX get_d2ham_dia(int i, vector<int>& id_);
  CMATRIX get_d2ham_dia(int i,int j);
  CMATRIX get_d2ham_dia(int i, int j, vector<int>& id_);

  // Adiabatic
  CMATRIX get_dc1_adi(int i);
  CMATRIX get_dc1_adi(int i, vector<int>& id_);
  CMATRIX get_ham_adi();
  CMATRIX get_ham_adi(vector<int>& id_);
  CMATRIX get_nac_adi();
  CMATRIX get_nac_adi(vector<int>& id_);
  CMATRIX get_hvib_adi();
  CMATRIX get_hvib_adi(vector<int>& id_);
  CMATRIX get_d1ham_adi(int i);
  CMATRIX get_d1ham_adi(int i, vector<int>& id_);
  CMATRIX get_d2ham_adi(int i);
  CMATRIX get_d2ham_adi(int i, vector<int>& id_);
  CMATRIX get_d2ham_adi(int i,int j);
  CMATRIX get_d2ham_adi(int i, int j, vector<int>& id_);

  // Transforms
  CMATRIX get_basis_transform();
  CMATRIX get_basis_transform(vector<int>& id_);

  CMATRIX get_time_overlap_adi();
  CMATRIX get_time_overlap_adi(vector<int>& id_);

  CMATRIX get_time_overlap_dia();
  CMATRIX get_time_overlap_dia(vector<int>& id_);

  vector<int> get_ordering_adi(); 
  vector<int> get_ordering_adi(vector<int>& id_); 

  CMATRIX get_cum_phase_corr();
  CMATRIX get_cum_phase_corr(vector<int>& id_);



  ///< In nHamiltonian_compute_diabatic.cpp

  void compute_diabatic(int model, vector<double>& q, vector<double>& params, int lvl); // for internal model types
  void compute_diabatic(int model, vector<double>& q, vector<double>& params); // for internal model types

//  void compute_diabatic(bp::object py_funct, bp::object q, bp::object params, int lvl); // for models defined in Python
//  void compute_diabatic(bp::object py_funct, bp::object q, bp::object params); // for models defined in Python
  void compute_diabatic(bp::object py_funct, MATRIX& q, bp::object params, int lvl); // for models defined in Python
  void compute_diabatic(bp::object py_funct, MATRIX& q, bp::object params); // for models defined in Python


  ///< In nHamiltonian_compute_ETHD.cpp

  void add_ethd_dia(const MATRIX& q, const MATRIX& invM, int der_lvl);
  void add_ethd_adi(const MATRIX& q, const MATRIX& invM, int der_lvl);

  void add_ethd3_dia(const MATRIX& q, const MATRIX& invM, double alp, int der_lvl);
  void add_ethd3_adi(const MATRIX& q, const MATRIX& invM, double alp, int der_lvl);
  void add_ethd3_adi(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet, int der_lvl);




  ///< In nHamiltonian_compute_adiabatic.cpp
  void update_ordering(vector<int>& perm_t, int lvl);
  void update_ordering(vector<int>& perm_t);

/*
  void apply_phase_corrections(CMATRIX* phase_corr, int lvl);
  void apply_phase_corrections(CMATRIX& phase_corr, int lvl);
  void apply_phase_corrections(CMATRIX* phase_corr);
  void apply_phase_corrections(CMATRIX& phase_corr);

  CMATRIX update_phases(CMATRIX& U_prev, int lvl);
  CMATRIX update_phases(CMATRIX& U_prev);
*/

  void compute_adiabatic(int der_lvl, int lvl);
  void compute_adiabatic(int der_lvl);
//  void compute_adiabatic(bp::object py_funct, bp::object q, bp::object params, int lvl); // for models defined in Python
//  void compute_adiabatic(bp::object py_funct, bp::object q, bp::object params); // for models defined in Python
  void compute_adiabatic(bp::object py_funct, MATRIX& q, bp::object params, int lvl); // for models defined in Python
  void compute_adiabatic(bp::object py_funct, MATRIX& q, bp::object params); // for models defined in Python



  ///< In nHamiltonian_compute_basis_transform.cpp

  void ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi);
  void ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi, vector<int>& id_);
  void ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split);
  void ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi);
  void ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi, vector<int>& id_);
  void ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split);

  void ampl_dia2adi(CMATRIX* ampl_dia, CMATRIX* ampl_adi);
  void ampl_dia2adi(CMATRIX* ampl_dia, CMATRIX* ampl_adi, vector<int>& id_);
  void ampl_dia2adi(CMATRIX* ampl_dia, CMATRIX* ampl_adi, int lvl, int split);
  void ampl_adi2dia(CMATRIX* ampl_dia, CMATRIX* ampl_adi);
  void ampl_adi2dia(CMATRIX* ampl_dia, CMATRIX* ampl_adi, vector<int>& id_);
  void ampl_adi2dia(CMATRIX* ampl_dia, CMATRIX* ampl_adi, int lvl, int split);




  ///< In nHamiltonian_compute_forces.cpp

  CMATRIX forces_adi(vector<int>& act_states);   // -dH_adi/dR in the adiabatic basis for several trajectories
  CMATRIX forces_dia(vector<int>& act_states);   // -dH_dia/dR in the diabatic basis for several trajectories

  CMATRIX all_forces_adi(vector<int>& id_);

  
/*
  These variables are poorly-defined and unnecessary, so we are going to deprecate them, AVA 9/15/2022
  Instead - there will be functions to extract adiabatic and diabatic forces


  vector<CMATRIX> forces_tens_adi(CMATRIX& ampl_adi); // 
  vector<CMATRIX> forces_tens_adi(CMATRIX& ampl_adi, vector<int>& id_); // 
  vector<CMATRIX> forces_tens_dia(CMATRIX& ampl_dia); // 
  vector<CMATRIX> forces_tens_dia(CMATRIX& ampl_dia, vector<int>& id_); // 

  CMATRIX forces_adi(CMATRIX& ampl_adi);  // -dE/dR in the adiabatic basis, assuming Cadi = Cadi(t)
  CMATRIX forces_adi(CMATRIX& ampl_adi, vector<int>& id_);  // -dE/dR in the adiabatic basis, assuming Cadi = Cadi(t)
  CMATRIX forces_adi(vector<int>& act_states);   // -dE/dR in the adiabatic basis for several trajectories
  CMATRIX forces_adi(int act_state);   // -dE/dR in the adiabatic basis for several trajectories, all in the same state

  CMATRIX forces_dia(CMATRIX& ampl_dia);  // -dE/dR in the diabatic basis, assuming Cdia = Cdia(t)
  CMATRIX forces_dia(CMATRIX& ampl_dia, vector<int>& id_);  // -dE/dR in the diabatic basis, assuming Cdia = Cdia(t)
  CMATRIX forces_dia(vector<int>& act_states);   // -dE/dR in the diabatic basis for several trajectories
  CMATRIX forces_dia(int act_state);   // -dE/dR in the diabatic basis for several trajectories, all in the same state
*/


  ///< In nHamiltonian_compute_nac.cpp

  void compute_nac_dia(MATRIX& p, const MATRIX& invM);
  void compute_nac_dia(MATRIX& p, const MATRIX& invM, vector<int>& id_);
  void compute_nac_dia(MATRIX& p, const MATRIX& invM, int lvl, int split);
  void compute_nac_adi(double dt, int method);
  void compute_nac_adi(MATRIX& p, const MATRIX& invM);
  void compute_nac_adi(MATRIX& p, const MATRIX& invM, vector<int>& id_);
  void compute_nac_adi(MATRIX& p, const MATRIX& invM, int lvl, int split);

  void compute_hvib_dia();
  void compute_hvib_dia(vector<int>& id_);
  void compute_hvib_dia(int lvl);
  void compute_hvib_adi();
  void compute_hvib_adi(vector<int>& id_);
  void compute_hvib_adi(int lvl);




  ///< In nHamiltonian_compute_Ehrenfest.cpp and nHamiltonian_compute_Ehrenfest_forces.cpp
  complex<double> Ehrenfest_energy_dia(CMATRIX& ampl_dia);
  complex<double> Ehrenfest_energy_dia(CMATRIX& ampl_dia, vector<int>& id_);
  CMATRIX Ehrenfest_forces_dia_unit(CMATRIX& ampl_dia, int option);      ///< Ehrenfest forces in diabatic basis
  CMATRIX Ehrenfest_forces_dia_unit(CMATRIX& ampl_dia);                  ///< Ehrenfest forces in diabatic basis
  CMATRIX Ehrenfest_forces_dia(CMATRIX& ampl_dia, int lvl, int option);  ///< Ehrenfest forces in diabatic basis
  CMATRIX Ehrenfest_forces_dia(CMATRIX& ampl_dia, int lvl);              ///< Ehrenfest forces in diabatic basis


  ///< In nHamiltonian_compute_Ehrenfest.cpp and nHamiltonian_compute_Ehrenfest_forces.cpp
  complex<double> Ehrenfest_energy_adi(CMATRIX& ampl_adi);
  complex<double> Ehrenfest_energy_adi(CMATRIX& ampl_adi, vector<int>& id_);
  CMATRIX Ehrenfest_forces_adi_unit(CMATRIX& ampl_adi, int option);      ///< Ehrenfest forces in adiabatic basis
  CMATRIX Ehrenfest_forces_adi_unit(CMATRIX& ampl_adi);                  ///< Ehrenfest forces in adiabatic basis
  CMATRIX Ehrenfest_forces_adi(CMATRIX& ampl_adi, int lvl, int option);  ///< Ehrenfest forces in adiabatic basis
  CMATRIX Ehrenfest_forces_adi(CMATRIX& ampl_adi, int lvl);              ///< Ehrenfest forces in adiabatic basis




  vector<CMATRIX> Ehrenfest_forces_tens_adi(CMATRIX& ampl_adi);  // Force tensor in adiabatic basis, assuming Cadi = Cadi(t)
  vector<CMATRIX> Ehrenfest_forces_tens_adi(CMATRIX& ampl_adi, vector<int>& id_);  // Force tensor in adiabatic basis, assuming Cadi = Cadi(t)
  vector<CMATRIX> Ehrenfest_forces_tens_dia(CMATRIX& ampl_dia);  // Force tensor in diabatic basis, assuming Cdia = Cdia(t)
  vector<CMATRIX> Ehrenfest_forces_tens_dia(CMATRIX& ampl_dia, vector<int>& id_);  // Force tensor in diabatic basis, assuming Cdia = Cdia(t)




  friend bool operator == (const nHamiltonian& h1, const nHamiltonian& h2){
    return &h1 == &h2;
  }
  friend bool operator != (const nHamiltonian& h1, const nHamiltonian& h2){
    return !(h1 == h2);  // only compare addresses
  }



};

typedef std::vector<nHamiltonian> nHamiltonianList;  ///< data type for keeping a list of generic Hamiltonians of their derived classes


///< In nHamiltonian_compute adiabatic
/*

*/


///< In nHamiltonian_compute_ETHD.cpp
double ETHD_energy(const MATRIX& q, const MATRIX& invM);
MATRIX ETHD_forces(const MATRIX& q, const MATRIX& invM);


///< In nHamiltonian_compute_ETHD3.cpp
double ETHD3_energy(const MATRIX& q, const MATRIX& invM, double alp);
double ETHD3_energy(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet);

MATRIX ETHD3_forces(const MATRIX& q, const MATRIX& invM, double alp);                                // -dH_{ETHD3}/dQ
MATRIX ETHD3_forces(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet);   // -dH_{ETHD3}/dQ

MATRIX ETHD3_friction(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet); // dH_{ETHD3}/dP




}// namespace libnhamiltonian
}// liblibra

#endif // nHAMILTONIAN_H
