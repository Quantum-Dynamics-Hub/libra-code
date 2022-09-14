/*********************************************************************************
* Copyright (C) 2021 Alexey V. Akimov
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

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "../math_linalg/liblinalg.h"


/// liblibra namespace
namespace liblibra{


using namespace liblinalg;

/// libdyn namespace
namespace libdyn{

namespace bp = boost::python;



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


  ///================= General variables, for OOP implementation ===================
  /**
    Status of the general vars

    0 - not allocated;
    1 - allocated
  */
  int gen_vars_status; 

  /**
    Nuclear coordinates 
    
    Options:
     MATRIX(ndof, ntraj)
  */
//  MATRIX* q; 


  /**
    Nuclear momenta
    
    Options:
     MATRIX(ndof, ntraj)
  */
//  MATRIX* p; 


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
//  vector<int> act_states;


  /**
    Projectors transforming dynamically-consistent and raw wavefunctions to each other
    
    Options:
     projectors(nstates, nstates) x ntraj
  */
//  CMATRIX** projectors;



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


  void allocate_gen_vars();
  void allocate_afssh();
  void allocate_bcsh();

  //dyn_variables();
  dyn_variables(int _ndia, int _nadi, int _ndof, int _ntraj);
  dyn_variables(const dyn_variables& x); 
  ~dyn_variables();

  void set_parameters(bp::dict params);


  friend bool operator == (const dyn_variables& n1, const dyn_variables& n2){
    return &n1 == &n2;
  }
  friend bool operator != (const dyn_variables& n1, const dyn_variables& n2){
    return !(n1 == n2);  // only compare addresses
  }


};



} // libdyn
}// liblibra

#endif // DYN_VARIABLES_H
