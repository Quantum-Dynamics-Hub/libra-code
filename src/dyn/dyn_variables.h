/*********************************************************************************
* Copyright (C) 2021-2023 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_variables.h
  \brief The file implements a class to store the dynamical variables of need for various methods
*/


#ifndef DYN_VARIABLES_H
#define DYN_VARIABLES_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "../math_linalg/liblinalg.h"
#include "../math_random/librandom.h"
#include "../math_specialfunctions/libspecialfunctions.h"
#include "../nhamiltonian/libnhamiltonian.h"
#include "dyn_control_params.h"



/// liblibra namespace
namespace liblibra{


using namespace liblinalg;
using namespace librandom;
using namespace libnhamiltonian;


/// libdyn namespace
namespace libdyn{

namespace bp = boost::python;



CMATRIX transform_amplitudes(int rep_in, int rep_out, CMATRIX& C, nHamiltonian& ham);


class dyn_variables{

  public:

  ///================= Dimension numbers ===================

  /**
    The number of diabatic states in the data dimension
    
    Options:
     any non-negative integer number
  */
  int ndia;


  /**
    The number of adiabatic states in the data dimension
    
    Options:
     any non-negative integer number
  */
  int nadi;


  /**
    The number of nuclear degrees of freedom in the data dimension
    
    Options:
     any non-negative integer number
  */
  int ndof;


  /**
    The number of trajectories in the data dimension
    
    Options:
     any non-negative integer number
  */
  int ntraj;


  ///================= Electronic variables, for OOP implementation ===================
  /**
    Status of the electronic vars

    0 - not allocated;
    1 - allocated
  */
  int electronic_vars_status; 

  /**
    Electronic amplitudes in diabatic representation
    
    Options:
     CMATRIX(ndia, ntraj)
  */
  CMATRIX* ampl_dia; 


  /**
    Electronic amplitudes in adiabatic representation
    
    Options:
     CMATRIX(nadi, ntraj)
  */
  CMATRIX* ampl_adi; 


  /**
    Cumulative projection matrices for adiabatic states:

    |psi_adi_ordered(t)> = |psi_adi_raw(t)> * proj_adi(t)

    Options:
    vector<ntraj, CMATRIX(nadi, nadi)>
  */
  vector<CMATRIX*> proj_adi;


  /**
    Electronic density matrix in diabatic representation
    
    Options:
     vector<ntraj, CMATRIX(ndia, ndia)>
  */
  vector<CMATRIX*> dm_dia; 


  /**
    Electronic density matrix in adiabatic representation
    
    Options:
     vector<ntraj, CMATRIX(nadi, nadi)>
  */
  vector<CMATRIX*> dm_adi; 


  /**
    Active states for each trajectory
    
    Options:
     vector<int> act_states(ntraj)
  */
  vector<int> act_states;


  ///================= Nuclear variables, for OOP implementation ===================
  /**
    Status of the nuclear vars

    0 - not allocated;
    1 - allocated
  */
  int nuclear_vars_status; 


  /**
    Inverse nuclear masses
    
    Options:
     MATRIX(ndof, 1)
  */
  MATRIX* iM; 


  /**
    Nuclear coordinates 
    
    Options:
     MATRIX(ndof, ntraj)
  */
  MATRIX* q; 


  /**
    Nuclear momenta
    
    Options:
     MATRIX(ndof, ntraj)
  */
  MATRIX* p;


  /**
    Nuclear forces (active)
    
    Options:
      MATRIX(ndof, ntraj)
  */
  MATRIX* f;



  ///================= For A-FSSH ===================
  /**
    Status of the A-FSSH vars

    0 - not allocated;
    1 - allocated
  */
  int afssh_vars_status; 

  /**
    Moments of coordinates in the adiabatic representation for all DOFs and all trajectories
    
    Options:
     CMATRIX(nadi, nadi) x ndof x ntraj, so delta_q[itraj][idof]->get(i,j)

    For Method: A-FSSH
  */
  //CMATRIX*** dR;
  vector< vector<CMATRIX*> > dR;


  /**
    Moments of momenta in the adiabatic representation for all DOFs and all trajectories
    
    Options:
     CMATRIX(nadi, nadi) x ndof x ntraj, so delta_p[itraj][idof]->get(i,j)

    For Method: A-FSSH
    
  */
  //CMATRIX*** dP;
  vector< vector<CMATRIX*> > dP;


  ///================= For BCSH ===================
  /**
    Status of the BCSH vars

    0 - not allocated;
    1 - allocated
  */
  int bcsh_vars_status; 

  /**
    Reversal event matrix
    
    Options:
     MATRIX(nadi, ntraj)

    For Method: BCSH
  */
  MATRIX* reversal_events;

  
  ///================= For DISH ===================
  /**
    Status of the DISH vars

    0 - not allocated;
    1 - allocated
  */
  int dish_vars_status;

  /**
    Coherence times

    Options:
     MATRIX(nadi, ntraj)

    For Method: DISH
  */
  MATRIX* coherence_time;

  ///================= For FSSH2 ===================
  /**
    Status of the FSSH2 vars

    0 - not allocated;
    1 - allocated
  */
  int fssh2_vars_status;

  /**
    Electronic density matrix in diabatic representation, at previos timestep

    Options:
     vector<ntraj, CMATRIX(ndia, ndia)>
  */
  vector<CMATRIX*> dm_dia_prev;


  /**
    Electronic density matrix in adiabatic representation, at previous timestep

    Options:
     vector<ntraj, CMATRIX(nadi, nadi)>
  */
  vector<CMATRIX*> dm_adi_prev;




  ///====================== In dyn_variables.cpp =====================

  void allocate_electronic_vars();
  void allocate_nuclear_vars();
  void allocate_afssh();
  void allocate_bcsh();
  void allocate_dish();
  void allocate_fssh2();

  dyn_variables(int _ndia, int _nadi, int _ndof, int _ntraj);
  dyn_variables(const dyn_variables& x); 
  ~dyn_variables();

  void set_parameters(bp::dict params);


  void set_q(MATRIX& _q){ *q = _q; }
  void set_p(MATRIX& _p){ *p = _p; }
  void set_f(MATRIX& _f){ *f = _f; }

  CMATRIX get_ampl_adi(){ return *ampl_adi; }
  CMATRIX get_ampl_dia(){ return *ampl_dia; }
  CMATRIX get_dm_adi(int i){  return *dm_adi[i]; }
  CMATRIX get_dm_dia(int i){  return *dm_dia[i]; }
  CMATRIX get_dm_adi(int i, int prev_steps);
  CMATRIX get_dm_dia(int i, int prev_steps);
  MATRIX get_imass(){ return *iM; }
  MATRIX get_coords(){ return *q; }
  MATRIX get_momenta(){ return *p; }
  MATRIX get_forces(){ return *f; }
  


  ///====================== In dyn_variables_nuclear.cpp =====================

  void init_nuclear_dyn_var(bp::dict _params, Random& rnd);
  double compute_average_kinetic_energy();
  double compute_average_kinetic_energy(vector<int>& which_dofs);
  double compute_kinetic_energy(int itraj);
  double compute_kinetic_energy(int itraj, vector<int>& which_dofs);
  vector<double> compute_kinetic_energies();
  vector<double> compute_kinetic_energies(vector<int>& which_dofs);


  ///====================== In dyn_variables_electronic.cpp =====================

  void update_amplitudes(dyn_control_params& dyn_params, nHamiltonian& ham);
  void update_amplitudes(bp::dict dyn_params, nHamiltonian& ham);
  void update_amplitudes(dyn_control_params& dyn_params, bp::object compute_model, bp::dict model_params);
  void update_amplitudes(bp::dict dyn_params, bp::object compute_model, bp::dict model_params);

  void update_density_matrix(dyn_control_params& dyn_params, nHamiltonian& ham, int lvl);
  void update_density_matrix(bp::dict dyn_params, nHamiltonian& ham, int lvl);
  void update_density_matrix(dyn_control_params& dyn_params, bp::object compute_model, bp::dict model_params, int lvl);
  void update_density_matrix(bp::dict dyn_params, bp::object compute_model, bp::dict model_params, int lvl);

  void update_active_states();

  void init_amplitudes(bp::dict params, Random& rnd);
  void init_density_matrix(bp::dict _params);
  void init_active_states(bp::dict _params, Random& rnd);

  void init_electronic_dyn_var(bp::dict params, Random& rnd);

  CMATRIX compute_average_dm(int rep);
  vector<double> compute_average_se_pop(int rep);
  vector<double> compute_average_sh_pop();

  void save_curr_dm_into_prev();



  friend bool operator == (const dyn_variables& n1, const dyn_variables& n2){
    return &n1 == &n2;
  }
  friend bool operator != (const dyn_variables& n1, const dyn_variables& n2){
    return !(n1 == n2);  // only compare addresses
  }


};


//vector<int> update_active_states(vector<int>& act_states, vector<CMATRIX*>& T);
CMATRIX orthogonalized_T(CMATRIX& T);


} // libdyn
}// liblibra

#endif // DYN_VARIABLES_H
