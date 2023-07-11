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
  \file Hamiltonian_Model.cpp
  \brief The file implements the basic methods of the model Hamiltonian class
    
*/


#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <complex>
#include <cmath>
#endif 

#include "Hamiltonian_Model.h"

/// liblibra namespace
namespace liblibra{

/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_model namespace
namespace libhamiltonian_model{

using namespace libhamiltonian_generic;
using namespace liblinalg;
using namespace libmeigen;
using std::complex;
using std::sin;
using std::cos;
using std::exp;
using std::sqrt;


Hamiltonian_Model::Hamiltonian_Model(int ham_indx_){
/**
  The constructor with parameter

  \param[in] ham_indx_ The parameter that selects the model Hamiltonian
  Possible options: 
  ham_indx      Description     
    0             SAC           
    1             DAC           
    2             ECWR          
    3             Marcus        
    4             SEXCH         
    5             Rabi2         
    6             --
    7             1D sin        
    8             cubic
    9             double_well
   200            2D sin        

*/

  int i;
  ham_indx = ham_indx_;

  // 2-level models
  if(ham_indx==8 || ham_indx==9){  // cubic, double_well
    nelec = 1;
    nnucl = 1;
  }

  else if(ham_indx==0 || ham_indx==1 || ham_indx==2 || ham_indx==3 
  || ham_indx==5 || ham_indx==6 || ham_indx==7 ){  // SAC, DAC, ECWR, Marcus, Rabi2, sin
    nelec = 2;
    nnucl = 1;
  }
  else if(ham_indx==4){    // superexchange
    nelec = 3;
    nnucl = 1;
  }
  else if(ham_indx==200){    // sin_2D
    nelec = 2;
    nnucl = 2;
  }


  ham_dia = new MATRIX(nelec,nelec); *ham_dia = 0.0;
  ham_adi = new MATRIX(nelec,nelec); *ham_adi = 0.0;
  basis_transform = new MATRIX(nelec, nelec);  basis_transform->Init_Unit_Matrix(1.0);  // trivial basis
  basis_overlap = new MATRIX(nelec, nelec); basis_overlap->Init_Unit_Matrix(1.0); // identity matrix
  

  d1ham_dia = vector<MATRIX*>(nnucl);
  d1ham_adi = vector<MATRIX*>(nnucl);
  d1basis_overlap = vector<MATRIX*>(nnucl);

  for(i=0;i<nnucl;i++){  
    d1ham_dia[i] = new MATRIX(nelec,nelec); *d1ham_dia[i] = 0.0; 
    d1ham_adi[i] = new MATRIX(nelec,nelec); *d1ham_adi[i] = 0.0; 
    d1basis_overlap[i] = new MATRIX(nelec,nelec); *d1basis_overlap[i] = 0.0;
  }

  d2ham_dia = vector<MATRIX*>(nnucl*nnucl);
  for(i=0;i<nnucl*nnucl;i++){  d2ham_dia[i] = new MATRIX(nelec,nelec); *d2ham_dia[i] = 0.0;  }


  rep = 0; // default representation is diabatic  
  status_dia = 0;
  status_adi = 0;

}

Hamiltonian_Model::~Hamiltonian_Model(){
/**
  Destructor
*/
  int i;

//  cout<<"Derived Ham_mod. destructor\n";
 
  delete ham_dia;
  delete ham_adi;
  delete basis_transform;
  delete basis_overlap;

  for(i=0;i<nnucl;i++){  
    delete d1ham_dia[i];
    delete d1ham_adi[i];
    delete d1basis_overlap[i];
  }
  for(i=0;i<nnucl*nnucl;i++){  delete d2ham_dia[i]; }

  if(d1ham_dia.size()>0) { d1ham_dia.clear(); }
  if(d2ham_dia.size()>0) { d2ham_dia.clear(); }
  if(d1ham_adi.size()>0) { d1ham_adi.clear(); }
  if(d1basis_overlap.size()>0) { d1basis_overlap.clear(); }


}


void Hamiltonian_Model::set_params(vector<double>& params_){
/**
  \param[in] params_ The array (vector of doubles) of parameters to be used. The expected
  maximal number of parameters is. Smaller # can be provided, but not more:

  ham_indx      Description      max # of parameters
    0             SAC                 4
    1             DAC                 3
    2             ECWR                4
    3             Marcus              3
    4             SEXCH               12
    5             Rabi2               3
    6               --               --
    7             1D sin              5
    8             cubic               3
    9             double_well         9
   200            2D sin              6

  Setting parameters of the model Hamiltonian.
  The number and the meaning of the parameters is different across different models
*/

  int num_params = 0;

  if(ham_indx==0 || ham_indx==2){ // SAC, ECWR
    num_params = 4;
  }
  else if(ham_indx==7){ num_params = 5; } // sin
  else if(ham_indx==5 || ham_indx==6 || ham_indx==8 ){ // Rabi2 or cubic
    num_params = 3;
  }
  else if(ham_indx==1 || ham_indx==3){  // DAC, Marcus
    num_params = 5;
  }
  else if(ham_indx==4){  // SEXCH
    num_params = 12;
  }

  else if(ham_indx==9) { //double_well
    num_params = 1;
  }

  else if(ham_indx==200){  // SEXCH
    num_params = 6;
  }


  params = vector<double>(num_params, 0.0);

//  cout<<"In Hamiltonian_Mode::set_params\n";
  // Now copy input params:
  for(int i=0;i<params_.size() && i<num_params; i++){
    params[i] = params_[i];

//    cout<<"params["<<i<<"] = "<< params_[i]<<endl;
  }

  // Since the parameters have changed - we need to recompute everything
  status_dia = 0;
  status_adi = 0;

}

void Hamiltonian_Model::set_basis_overlap_by_ref(MATRIX& basis_overlap_){
/**
  Sets the internal basis_overlap matrix to the externally-defined value
  This is a setting by reference, so it works in Python. 
  The input object should exist for the entire duratiton of the simulation
  that need using the basis_overlap matrix in the calculations.
  Moreover, if one changes the external object, this will affect the calculations
  since the internal pointer points to the external object.
*/
  if( basis_overlap_.n_cols!=nelec ){
    cout<<"Error in Hamiltonian_Model::set_basis_overlap_by_ref: the number of \
    columns in the input matrix ( "<<basis_overlap_.n_cols<<" ) does not match\
    the number of electronic states used in the present model ( "<<nelec<<") \nExiting now...\n";
    exit(0); 
  }
  if( basis_overlap_.n_rows!=nelec ){
    cout<<"Error in Hamiltonian_Model::set_basis_overlap_by_ref: the number of \
    rows in the input matrix ( "<<basis_overlap_.n_rows<<" ) does not match\
    the number of electronic states used in the present model ( "<<nelec<<") \nExiting now...\n";
    exit(0); 
  }
  
  basis_overlap = &basis_overlap_;

}

void Hamiltonian_Model::set_basis_overlap_by_val(MATRIX& basis_overlap_){
/**
  Sets the internal basis_overlap matrix to the externally-defined value
  This is a setting by value, the input object doesn't have to exist 
  for the entire duratiton of the simulation. 
  Moreover, if one changes the external object, this does not affect the calculations
  unless this value is fed back to the Hamiltonian using this function.
*/

  if( basis_overlap_.n_cols!=nelec ){
    cout<<"Error in Hamiltonian_Model::set_basis_overlap_by_ref: the number of \
    columns in the input matrix ( "<<basis_overlap_.n_cols<<" ) does not match\
    the number of electronic states used in the present model ( "<<nelec<<") \nExiting now...\n";
    exit(0); 
  }
  if( basis_overlap_.n_rows!=nelec ){
    cout<<"Error in Hamiltonian_Model::set_basis_overlap_by_ref: the number of \
    rows in the input matrix ( "<<basis_overlap_.n_rows<<" ) does not match\
    the number of electronic states used in the present model ( "<<nelec<<") \nExiting now...\n";
    exit(0); 
  }
  
  *basis_overlap = basis_overlap_;

}




void Hamiltonian_Model::set_basis_overlap_derivs_by_ref(vector<MATRIX>& basis_overlap_derivs_, vector<int>& indxs_){
/**
*/

}

void Hamiltonian_Model::set_basis_overlap_derivs_by_val(vector<MATRIX>& basis_overlap_derivs_, vector<int>& indxs_){
/**
*/

}





//void Hamiltonian_Model::set_q(vector<double>& q_){ Hamiltonian::set_q(q_); }
//void Hamiltonian_Model::set_q(boost::python::list q_){ Hamiltonian::set_q(q_); }
//void Hamiltonian_Model::set_v(vector<double>& v_){ Hamiltonian::set_v(v_); }
//void Hamiltonian_Model::set_v(boost::python::list v_){ Hamiltonian::set_v(v_); }



void Hamiltonian_Model::compute_diabatic(){
/**
  Compute diabatic Hamiltonian, its 1-st and 2-nd order derivatives w.r.t. nuclear DOF
  Computations are done only if they are not up to date (status_dia = 0)
*/

  if(status_dia == 0){ // only compute this if the result is not up to date
  
    double e;
    double x = q[0];
    double y = 0.0;  if(q.size()>=2){ y = q[1]; }


    //=============== 1D models ========================

    if(ham_indx==0){  SAC_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);     }    // SAC potetnial
    else if(ham_indx==1){ DAC_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params); }    // DAC potential
    else if(ham_indx==2){ ECWR_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);   } // ECWR potential
    else if(ham_indx==3){ Marcus_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params); } // Marcus potential
    else if(ham_indx==4){ SEXCH_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);  } // SEXCH potential
    else if(ham_indx==5){ Rabi2_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);  } // Rabi2 potential
    else if(ham_indx==6){ ;; }// Nothing here
    else if(ham_indx==7){ sin_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);  } // sin potential
    else if(ham_indx==8){ cubic_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);  } // cubic potential
    else if(ham_indx==9){ double_well_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);  } // double_well potential

    //=============== 2D models ========================
    else if(ham_indx==200){  sin_2D_Ham(x, y, ham_dia, d1ham_dia[0], d1ham_dia[1], d2ham_dia[0], d2ham_dia[3], params);     }    // sin potetnial



    // Set status flag 
    status_dia = 1;

  }//   status_dia == 0

}


void Hamiltonian_Model::compute_adiabatic(){
/**
  Compute the adiabatic Hamiltonian, its 1-st and 2-nd order derivatives w.r.t. nuclear DOF
  Computations are done only if they are not up to date (status_dia = 0)
*/

// This function computes adiabatic PESs and derivative couplings

  compute_diabatic();

  if(status_adi == 0){

    // Transformation to adiabatic basis by solving the generalized eigenvalue problem
    // H_dia * C = S * C * H_adi, where C = basis_transform, S = basis_overlap
    // Here, of course, we assume that H_dia does not depend on C!
    solve_eigen(ham_dia, basis_overlap, ham_adi, basis_transform, 0);  


    // Now compute the derivative couplings (off-diagonal, multiplied by energy difference) and adiabatic gradients (diagonal)
    for(int n=0;n<nnucl;n++){

      *d1ham_adi[n] = (*basis_transform).T() * (*d1ham_dia[n] - *ham_adi * *basis_overlap ) * (*basis_transform);

    }// for n

    // Set status flag
    status_adi = 1;

  }// status_adi == 0
  
}


}// namespace libhamiltonian_model
}// namespace libhamiltonian
}// liblibra



