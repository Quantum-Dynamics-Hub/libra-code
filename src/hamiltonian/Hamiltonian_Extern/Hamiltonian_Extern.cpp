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


#include <complex>
#include <cmath>
#include "Hamiltonian_Extern.h"


namespace libhamiltonian{
namespace libhamiltonian_extern{

using namespace libmmath;
using namespace libmmath::libmeigen;
using std::complex;
using std::sin;
using std::cos;
using std::exp;
using std::sqrt;


Hamiltonian_Extern::Hamiltonian_Extern(int _nelec, int _nnucl){

  cout<<"Warning: Hamiltonian_Extern does not allocate any internal memory, so make sure \
  you create external objects for all variables: ham_dia, ham_adi, d1ham_dia, d1ham_adi, d2ham_dia and bind them \
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

  bs_ham_dia = 0;
  bs_d1ham_dia = 0;
  bs_d2ham_dia = 0;
  bs_ham_adi = 0;
  bs_d1ham_adi = 0;


}

void Hamiltonian_Extern::set_adiabatic_opt(int ad_opt){

  if(ad_opt==0 || ad_opt==1){ adiabatic_opt = ad_opt;  }
  else{
    cout<<"Error in Hamiltonian_Extern::set_adiabatic_opt - allowed values are:\n";
    cout<<"  0 - use provided adiabatic Hamiltonian and its derivatives directly (default)\n";
    cout<<"  1 - use provided diabatic Hamiltonian and its derivatives, and the transformation from diabatic to adiabatic basis\n";
    exit(0);
  }

}

void Hamiltonian_Extern::bind_ham_dia(MATRIX& _ham_dia){
  if(_ham_dia.num_of_cols!=nelec){
    cout<<"Error in Hamiltonian_Extern::bind_ham_dia\n";
    cout<<"Expected number of electronic DOF = "<<nelec<<" the number of cols in the input matrix is = "<<_ham_dia.num_of_cols<<"\n";
    exit(0);
  }
  if(_ham_dia.num_of_rows!=nelec){
    cout<<"Error in Hamiltonian_Extern::bind_ham_dia\n";
    cout<<"Expected number of electronic DOF = "<<nelec<<" the number of rows in the input matrix is = "<<_ham_dia.num_of_rows<<"\n";
    exit(0);
  }

  // At this point, all is ok

  ham_dia = &_ham_dia;
  bs_ham_dia = 1;
}

void Hamiltonian_Extern::bind_d1ham_dia(vector<MATRIX>& _d1ham_dia){
  int sz = _d1ham_dia.size();

  if(sz!=nnucl){
    cout<<"Error in Hamiltonian_Extern::bind_d1ham_dia\n";
    cout<<"Expected number of nuclear DOF = "<<nnucl<<" the number of provided matrices = "<<sz<<"\n";
    exit(0);
  }
  if(sz>0){

    if(_d1ham_dia[0].num_of_cols!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d1ham_dia\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of cols in the input matrix is = "<<_d1ham_dia[0].num_of_cols<<"\n";
      exit(0);
    }
    if(_d1ham_dia[0].num_of_rows!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d1ham_dia\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of rows in the input matrix is = "<<_d1ham_dia[0].num_of_rows<<"\n";
      exit(0);
    }
  }

  // At this point, all is ok

  for(int i=0; i<sz; i++){  d1ham_dia[i] = &_d1ham_dia[i];  }
  bs_d1ham_dia = 1;

}

void Hamiltonian_Extern::bind_d2ham_dia(vector<MATRIX>& _d2ham_dia){
  int sz = _d2ham_dia.size();

  if(sz!=nnucl*nnucl){
    cout<<"Error in Hamiltonian_Extern::bind_d2ham_dia\n";
    cout<<"Expected number of elements in _d2ham_dia = "<<nnucl*nnucl<<" the number of provided matrices = "<<sz<<"\n";
    exit(0);
  }
  if(sz>0){

    if(_d2ham_dia[0].num_of_cols!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d2ham_dia\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of cols in the input matrix is = "<<_d2ham_dia[0].num_of_cols<<"\n";
      exit(0);
    }
    if(_d2ham_dia[0].num_of_rows!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d2ham_dia\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of rows in the input matrix is = "<<_d2ham_dia[0].num_of_rows<<"\n";
      exit(0);
    }
  }

  // At this point, all is ok

  for(int i=0; i<sz; i++){  d2ham_dia[i] = &_d2ham_dia[i];  }
  bs_d2ham_dia = 1;

}


void Hamiltonian_Extern::bind_ham_adi(MATRIX& _ham_adi){
  if(_ham_adi.num_of_cols!=nelec){
    cout<<"Error in Hamiltonian_Extern::bind_ham_adi\n";
    cout<<"Expected number of electronic DOF = "<<nelec<<" the number of cols in the input matrix is = "<<_ham_adi.num_of_cols<<"\n";
    exit(0);
  }
  if(_ham_adi.num_of_rows!=nelec){
    cout<<"Error in Hamiltonian_Extern::bind_ham_adi\n";
    cout<<"Expected number of electronic DOF = "<<nelec<<" the number of rows in the input matrix is = "<<_ham_adi.num_of_rows<<"\n";
    exit(0);
  }

  // At this point, all is ok

  ham_adi = &_ham_adi;
  bs_ham_adi = 1;
}

void Hamiltonian_Extern::bind_d1ham_adi(vector<MATRIX>& _d1ham_adi){
  int sz = _d1ham_adi.size();

  if(sz!=nnucl){
    cout<<"Error in Hamiltonian_Extern::bind_d1ham_adi\n";
    cout<<"Expected number of nuclear DOF = "<<nnucl<<" the number of provided matrices = "<<sz<<"\n";
    exit(0);
  }
  if(sz>0){

    if(_d1ham_adi[0].num_of_cols!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d1ham_adi\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of cols in the input matrix is = "<<_d1ham_adi[0].num_of_cols<<"\n";
      exit(0);
    }
    if(_d1ham_adi[0].num_of_rows!=nelec){
      cout<<"Error in Hamiltonian_Extern::bind_d1ham_adi\n";
      cout<<"Expected number of electronic DOF = "<<nelec<<" the number of rows in the input matrix is = "<<_d1ham_adi[0].num_of_rows<<"\n";
      exit(0);
    }
  }

  // At this point, all is ok

  for(int i=0; i<sz; i++){  d1ham_adi[i] = &_d1ham_adi[i];  }
  bs_d1ham_adi = 1;

}




Hamiltonian_Extern::~Hamiltonian_Extern(){
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

  if(status_dia == 0){ // only compute this is the result if not up to date

    // setup ham_dia, d1ham_dia, and d2ham_dia matrices, so below we will basically 
    // check regarding the status of bindings

    if(bs_ham_dia == 0){  
      cout<<"Error in Hamiltonian_Extern::compute_diabatic\n";
      cout<<"Diabatic Hamiltonian has not been bound to the Hamiltonian_Extern object\n";
      cout<<"us \"bind_ham_dia\" function\n";
      exit(0);
    }
    if(bs_d1ham_dia == 0){  
      cout<<"Error in Hamiltonian_Extern::compute_diabatic\n";
      cout<<"Derivatives of diabatic Hamiltonian have not been bound to the Hamiltonian_Extern object\n";
      cout<<"us \"bind_d1ham_dia\" function\n";
      exit(0);
    }

    // Now it is ok, so far we don't care about d2ham_dia matrices!!!  

    // Set status flag 
    status_dia = 1;

  }//   status_dia == 0

}


void Hamiltonian_Extern::compute_adiabatic(){
// This function computes adiabatic PESs and derivative couplings
// Here, we have 2 options:
//  - adiabatic Hamiltonian and its derivatives are bound - we simply use them
//  - diabatic Hamiltonian and its derivatives are bound (and also adiabatic, but they just used a memory arrays) -
//    - we need to perfrom diabatic -> adiabatic transformation, like in the Hamiltonian_Model case
// To distinguish these two cases, we use additional parameter - "adiabatic_opt"
 
  if(adiabatic_opt==0){ 
    // Everything is done already - don't do anything special
  }

  else if(adiabatic_opt==1){

    compute_diabatic();

    if(status_adi == 0){

      MATRIX* S; S = new MATRIX(nelec, nelec);  S->Init_Unit_Matrix(1.0);
      MATRIX* C; C = new MATRIX(nelec, nelec);  *C = 0.0;

      // Transformation to adiabatic basis
      solve_eigen(nelec, ham_dia, S, ham_adi, C);  // H_dia * C = S * C * H_adi

      // Now compute the derivative couplings (off-diagonal, multiplied by energy difference) and adiabatic gradients (diagonal)
      for(int n=0;n<nnucl;n++){

        *d1ham_adi[n] = (*C).T() * (*d1ham_dia[n]) * (*C);

      }// for n

      delete S;
      delete C;

      // Set status flag
      status_adi = 1;

    }// status_adi == 0

  }// adiabatic_opt
  
}


}// namespace libhamiltonian_extern
}// namespace libhamiltonian




