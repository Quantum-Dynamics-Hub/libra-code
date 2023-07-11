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
  \file Hamiltonian_Extern.cpp
  \brief The file implements the external Hamiltonian class - for interface with 3-rd party codes    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <complex>
#include <cmath>
#endif 
#include "Hamiltonian_Extern.h"

/// liblibra namespace
namespace liblibra{

/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_extern namespace
namespace libhamiltonian_extern{

using namespace liblinalg;
using namespace libmeigen;
using std::complex;
using std::sin;
using std::cos;
using std::exp;
using std::sqrt;


Hamiltonian_Extern::Hamiltonian_Extern(int _nelec, int _nnucl){
/**
  \param[in] _nelec The number of electronic basis states 
  \param[in] _nnucl The number of nuclear DOF

  This Hamiltonian_Extern constructor does not allocate any internal memory, so make sure 
  you create external objects for all variables: ham_dia, ham_adi, d1ham_dia, d1ham_adi, d2ham_dia and bind them
  to the internal variables after the Hamiltoninan_Extern object is created
*/

  cout<<"Warning: Hamiltonian_Extern does not allocate any internal memory, so make sure \
  you create external objects for all variables: ham_dia, ham_adi, d1ham_dia, d1ham_adi, d2ham_dia, ham_vib, and bind them \
  to the internal variables\n";

  nelec = _nelec;
  nnucl = _nnucl;

  d1ham_dia = vector<MATRIX*>(nnucl);
  d1ham_adi = vector<MATRIX*>(nnucl);
  d2ham_dia = vector<MATRIX*>(nnucl*nnucl);

/*
  int i;

  ham_dia = new MATRIX(nelec,nelec); *ham_dia = 0.0;
  ham_adi = new MATRIX(nelec,nelec); *ham_adi = 0.0;

  d1ham_dia = vector<MATRIX*>(nnucl);
  d1ham_adi = vector<MATRIX*>(nnucl);
  for(i=0;i<nnucl;i++){  
    d1ham_dia[i] = new MATRIX(nelec,nelec); *d1ham_dia[i] = 0.0; 
    d1ham_adi[i] = new MATRIX(nelec,nelec); *d1ham_adi[i] = 0.0; 
  }

  d2ham_dia = vector<MATRIX*>(nnucl*nnucl);
  for(i=0;i<nnucl*nnucl;i++){  d2ham_dia[i] = new MATRIX(nelec,nelec); *d2ham_dia[i] = 0.0;  }
*/

  rep = 0; // default representation is diabatic  
  status_dia = 0;
  status_adi = 0;


  adiabatic_opt = 0;
  vibronic_opt = 0;

  bs_ham_dia = 0;
  bs_d1ham_dia = 0;
  bs_d2ham_dia = 0;
  bs_ham_adi = 0;
  bs_d1ham_adi = 0;
  bs_ham_vib = 0;

}

void Hamiltonian_Extern::set_adiabatic_opt(int ad_opt){
/**
  \param[in] ad_opt The new 'adiabatic_opt' value to set:
  0 - use provided adiabatic Hamiltonian and its derivatives directly (default)
  1 - use provided diabatic Hamiltonian and its derivatives, and the transformation from diabatic to adiabatic basis

  Sets the adiabatic_opt variable
*/

  if(ad_opt==0 || ad_opt==1){ adiabatic_opt = ad_opt;  }
  else{
    cout<<"Error in Hamiltonian_Extern::set_adiabatic_opt - allowed values are:\n";
    cout<<"  0 - use provided adiabatic Hamiltonian and its derivatives directly (default)\n";
    cout<<"  1 - use provided diabatic Hamiltonian and its derivatives, and the transformation from diabatic to adiabatic basis\n";
    exit(0);
  }

}

void Hamiltonian_Extern::set_vibronic_opt(int vib_opt){
/**
  \param[in] vib_opt The new 'vibronic_opt' value to set:
  0 - use provided vibronic Hamiltonian (default)
  1 - use provided diabatic and/or adiabatic Hamiltonian to compute vibronic Hamiltonian

  Sets the vibronic_opt variable
*/

  if(vib_opt==0 || vib_opt==1){ vibronic_opt = vib_opt;  }
  else{
    cout<<"Error in Hamiltonian_Extern::set_vibronic_opt - allowed values are:\n";
    cout<<"  0 - use provided vibronic Hamiltonian (default)\n";
    cout<<"  1 - use provided diabatic and/or adiabatic Hamiltonian to compute vibronic Hamiltonian\n";
    exit(0);
  }

}


void Hamiltonian_Extern::bind_ham_dia(MATRIX& _ham_dia){
/**
  \param[in] _ham_dia The external matrix containing the diabatic Hamiltonian

  Makes the internal pointer, ham_dia, to point to the external object containing diabatic Hamiltonian matrix
*/

  if(_ham_dia.n_cols!=nelec){
    cout<<"Error in Hamiltonian_Extern::bind_ham_dia\n";
    cout<<"Expected number of electronic DOF = "<<nelec<<" the number of cols in the input matrix is = "<<_ham_dia.n_cols<<"\n";
    exit(0);
  }
  if(_ham_dia.n_rows!=nelec){
    cout<<"Error in Hamiltonian_Extern::bind_ham_dia\n";
    cout<<"Expected number of electronic DOF = "<<nelec<<" the number of rows in the input matrix is = "<<_ham_dia.n_rows<<"\n";
    exit(0);
  }

  // At this point, all is ok

  ham_dia = &_ham_dia;
  bs_ham_dia = 1;
}

void Hamiltonian_Extern::bind_d1ham_dia(vector<MATRIX>& _d1ham_dia){
/**
  \param[in] _d1ham_dia The external vector of matrices containing the 1-st order derivatives of diabatic Hamiltonian
  w.r.t. all nuclear DOF

  Makes the internal pointers, d1ham_dia[i], to point to the external objects containing the 1-st order 
  derivatives of diabatic Hamiltonian matrix w.r.t. nuclear DOF
*/

  int sz = _d1ham_dia.size();

  if(sz!=nnucl){
    cout<<"Error in Hamiltonian_Extern::bind_d1ham_dia\n";
    cout<<"Expected number of nuclear DOF = "<<nnucl<<" the number of provided matrices = "<<sz<<"\n";
    exit(0);
  }
  if(sz>0){

    if(_d1ham_dia[0].n_cols!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d1ham_dia\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of cols in the input matrix is = "<<_d1ham_dia[0].n_cols<<"\n";
      exit(0);
    }
    if(_d1ham_dia[0].n_rows!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d1ham_dia\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of rows in the input matrix is = "<<_d1ham_dia[0].n_rows<<"\n";
      exit(0);
    }
  }

  // At this point, all is ok

  for(int i=0; i<sz; i++){  d1ham_dia[i] = &_d1ham_dia[i];  }
  bs_d1ham_dia = 1;

}

void Hamiltonian_Extern::bind_d2ham_dia(vector<MATRIX>& _d2ham_dia){
/**
  \param[in] _d2ham_dia The external vector of matrices containing the 2-nd order derivatives of diabatic Hamiltonian
  w.r.t. all nuclear DOF

  Makes the internal pointers, d2ham_dia[i], to point to the external objects containing the 2-nd order 
  derivatives of diabatic Hamiltonian matrix w.r.t. nuclear DOF
*/

  int sz = _d2ham_dia.size();

  if(sz!=nnucl*nnucl){
    cout<<"Error in Hamiltonian_Extern::bind_d2ham_dia\n";
    cout<<"Expected number of elements in _d2ham_dia = "<<nnucl*nnucl<<" the number of provided matrices = "<<sz<<"\n";
    exit(0);
  }
  if(sz>0){

    if(_d2ham_dia[0].n_cols!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d2ham_dia\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of cols in the input matrix is = "<<_d2ham_dia[0].n_cols<<"\n";
      exit(0);
    }
    if(_d2ham_dia[0].n_rows!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d2ham_dia\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of rows in the input matrix is = "<<_d2ham_dia[0].n_rows<<"\n";
      exit(0);
    }
  }

  // At this point, all is ok

  for(int i=0; i<sz; i++){  d2ham_dia[i] = &_d2ham_dia[i];  }
  bs_d2ham_dia = 1;

}


void Hamiltonian_Extern::bind_ham_adi(MATRIX& _ham_adi){
/**
  \param[in] _ham_adi The external matrix containing the adiabatic Hamiltonian

  Makes the internal pointer, ham_adi, to point to the external object containing adiabatic Hamiltonian matrix
*/

  if(_ham_adi.n_cols!=nelec){
    cout<<"Error in Hamiltonian_Extern::bind_ham_adi\n";
    cout<<"Expected number of electronic DOF = "<<nelec<<" the number of cols in the input matrix is = "<<_ham_adi.n_cols<<"\n";
    exit(0);
  }
  if(_ham_adi.n_rows!=nelec){
    cout<<"Error in Hamiltonian_Extern::bind_ham_adi\n";
    cout<<"Expected number of electronic DOF = "<<nelec<<" the number of rows in the input matrix is = "<<_ham_adi.n_rows<<"\n";
    exit(0);
  }

  // At this point, all is ok

  ham_adi = &_ham_adi;
  bs_ham_adi = 1;
}

void Hamiltonian_Extern::bind_d1ham_adi(vector<MATRIX>& _d1ham_adi){
/**
  \param[in] _d1ham_adi The external vector of matrices containing the 1-st order derivatives of adiabatic Hamiltonian
  w.r.t. all nuclear DOF

  Makes the internal pointers, d1ham_adi[i], to point to the external objects containing the 1-st order 
  derivatives of adiabatic Hamiltonian matrix w.r.t. nuclear DOF
*/

  int sz = _d1ham_adi.size();

  if(sz!=nnucl){
    cout<<"Error in Hamiltonian_Extern::bind_d1ham_adi\n";
    cout<<"Expected number of nuclear DOF = "<<nnucl<<" the number of provided matrices = "<<sz<<"\n";
    exit(0);
  }
  if(sz>0){

    if(_d1ham_adi[0].n_cols!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d1ham_adi\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of cols in the input matrix is = "<<_d1ham_adi[0].n_cols<<"\n";
      exit(0);
    }
    if(_d1ham_adi[0].n_rows!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d1ham_adi\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of rows in the input matrix is = "<<_d1ham_adi[0].n_rows<<"\n";
      exit(0);
    }
  }

  // At this point, all is ok

  for(int i=0; i<sz; i++){  d1ham_adi[i] = &_d1ham_adi[i];  }
  bs_d1ham_adi = 1;

}


void Hamiltonian_Extern::bind_ham_vib(CMATRIX& _ham_vib){
/**
  \param[in] _ham_vib The external matrix containing the vibronic Hamiltonian

  Makes the internal pointer (defined in the derived Hamiltonian_Extern class), ham_vib, to point to the external 
  object containing the vibronic Hamiltonian matrix
*/

  if(_ham_vib.n_cols!=nelec){
    cout<<"Error in Hamiltonian_Extern::bind_ham_vib\n";
    cout<<"Expected number of electronic DOF = "<<nelec<<" the number of cols in the input matrix is = "<<_ham_vib.n_cols<<"\n";
    exit(0);
  }
  if(_ham_vib.n_rows!=nelec){
    cout<<"Error in Hamiltonian_Extern::bind_ham_vib\n";
    cout<<"Expected number of electronic DOF = "<<nelec<<" the number of rows in the input matrix is = "<<_ham_vib.n_rows<<"\n";
    exit(0);
  }

  // At this point, all is ok

  ham_vib = &_ham_vib;
  bs_ham_vib = 1;
}



Hamiltonian_Extern::~Hamiltonian_Extern(){
/** 
  Destructor - does nothing at all at this point
*/
  int i;

/* 
  delete ham_dia;
  delete ham_adi;

  for(i=0;i<nnucl;i++){  
    delete d1ham_dia[i];
    delete d1ham_adi[i];
  }
  for(i=0;i<nnucl*nnucl;i++){  delete d2ham_dia[i]; }

  if(d1ham_dia.size()>0) { d1ham_dia.clear(); }
  if(d2ham_dia.size()>0) { d2ham_dia.clear(); }
  if(d1ham_adi.size()>0) { d1ham_adi.clear(); }

*/

}



void Hamiltonian_Extern::compute_diabatic(){
/**
  Performs the "computations" of the diabatic Hamiltonian and its derivatives
  Basically, we only check if the ham_dia and d1ham_dia have been bound to the
  external variables, which contain the results.
  If everything is bound, then update the status of diabatic computations,
  if not - there appropriate error message will be printed out
  
*/

  if(status_dia == 0){ // only compute this is the result if not up to date

    // setup ham_dia, d1ham_dia, and d2ham_dia matrices, so below we will basically 
    // check regarding the status of bindings

    if(bs_ham_dia == 0){  
      cout<<"Error in Hamiltonian_Extern::compute_diabatic\n";
      cout<<"Diabatic Hamiltonian has not been bound to the Hamiltonian_Extern object\n";
      cout<<"use \"bind_ham_dia\" function\n";
      exit(0);
    }
    if(bs_d1ham_dia == 0){  
      cout<<"Error in Hamiltonian_Extern::compute_diabatic\n";
      cout<<"Derivatives of diabatic Hamiltonian have not been bound to the Hamiltonian_Extern object\n";
      cout<<"use \"bind_d1ham_dia\" function\n";
      exit(0);
    }

    // Now it is ok, so far we don't care about d2ham_dia matrices!!!  

    // Set status flag 
    status_dia = 1;

  }//   status_dia == 0

}


void Hamiltonian_Extern::compute_adiabatic(){
/** This function "computes" adiabatic PESs (energies and derivatives) and derivative couplings
  Here, we have 2 options:
  - adiabatic Hamiltonian and its derivatives are bound - we simply use them ("adiabatic_opt==0")
  - diabatic Hamiltonian and its derivatives are bound (bud we use "adiabatic_opt ==1" ) -
    - we need to perfrom diabatic -> adiabatic transformation, like in the Hamiltonian_Model case

 To distinguish these two cases, we use additional parameter - "adiabatic_opt"
*/
 
  if(adiabatic_opt==0){ 
    // Everything is done already - don't do anything special, other than check the bindings

    if(status_adi == 0){ // only compute this is the result if not up to date

      // setup ham_adi, d1ham_adi, and d2ham_adi matrices, so below we will basically 
      // check regarding the status of bindings

      if(bs_ham_adi == 0){  
        cout<<"Error in Hamiltonian_Extern::compute_adiabatic (with option adiabatic_opt == 0)\n";
        cout<<"Adiabatic Hamiltonian has not been bound to the Hamiltonian_Extern object\n";
        cout<<"use \"bind_ham_adi\" function\n";
        exit(0);
      }
      if(bs_d1ham_adi == 0){  
        cout<<"Error in Hamiltonian_Extern::compute_adiabatic (with option adiabatic_opt == 0)\n";
        cout<<"Derivatives of adiabatic Hamiltonian have not been bound to the Hamiltonian_Extern object\n";
        cout<<"use \"bind_d1ham_dia\" function\n";
        exit(0);
      }

      // Now it is ok, so far we don't care about d2ham_dia matrices!!!  

      // Set status flag 
      status_adi = 1;

    }//   status_dia == 0

  }// adiabatic_opt == 0

  else if(adiabatic_opt==1){

    compute_diabatic();

    if(status_adi == 0){

      MATRIX* S; S = new MATRIX(nelec, nelec);  S->Init_Unit_Matrix(1.0);
      MATRIX* C; C = new MATRIX(nelec, nelec);  *C = 0.0;

      // Transformation to adiabatic basis
      solve_eigen(ham_dia, S, ham_adi, C, 0);  // H_dia * C = S * C * H_adi

      // Now compute the derivative couplings (off-diagonal, multiplied by energy difference) and adiabatic gradients (diagonal)
      for(int n=0;n<nnucl;n++){

        *d1ham_adi[n] = (*C).T() * (*d1ham_dia[n]) * (*C);

      }// for n

      delete S;
      delete C;

      // Set status flag
      status_adi = 1;

    }// status_adi == 0

  }// adiabatic_opt == 1
  
}

std::complex<double> Hamiltonian_Extern::Hvib(int i,int j){
/**
  Return the vibronic Hamiltonian matrix element

  The returned Hamiltonian depends on the settings of the Hamiltonian_Extern object 
  This function does not invoke actual computation - it only returns whatever exists in the internal variables.

  \param[in] i index of electronic state
  \param[in] j index of electronic state
*/

  if(vibronic_opt == 0){ 
    // Everything is done already - don't do anything special, other than check the bindings

    // check regarding the status of bindings
    if(bs_ham_vib == 0){  
      cout<<"Error in Hamiltonian_Extern::Hvib (with option vibronic_opt == 0)\n";
      cout<<"Vibronic Hamiltonian has not been bound to the Hamiltonian_Extern object\n";
      cout<<"use \"bind_ham_vib\" function\n";
      exit(0);
    }

    return ham_vib->get(i,j);

  }// vibronic_opt == 0

  else if(vibronic_opt == 1){

    return Hamiltonian::Hvib(i,j);

  }

}


}// namespace libhamiltonian_extern
}// namespace libhamiltonian
}// liblibra



