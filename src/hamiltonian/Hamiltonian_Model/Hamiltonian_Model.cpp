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


#include <complex>
#include <cmath>
#include "Hamiltonian_Model.h"


namespace libhamiltonian{
namespace libhamiltonian_model{

using namespace libmmath;
using namespace libmmath::libmeigen;
using std::complex;
using std::sin;
using std::cos;
using std::exp;
using std::sqrt;


Hamiltonian_Model::Hamiltonian_Model(int ham_indx_){

//cout<<"Derived Ham_mod. constructor\n";

  int i;
  ham_indx = ham_indx_;

  // 2-level models
  if(ham_indx==0 || ham_indx==1 || ham_indx==2 || ham_indx==3 || ham_indx==5){  // SAC, DAC, ECWR, Marcus, Rabi2
    nelec = 2;
    nnucl = 1;
  }
  else if(ham_indx==4){    // superexchange
    nelec = 3;
    nnucl = 1;
  }

  ham_dia = new MATRIX(nelec,nelec); *ham_dia = 0.0;
  ham_adi = new MATRIX(nelec,nelec); *ham_adi = 0.0;

  d1ham_dia = vector<MATRIX*>(nnucl);
  d2ham_dia = vector<MATRIX*>(nnucl);
  d1ham_adi = vector<MATRIX*>(nnucl);
  for(i=0;i<nnucl;i++){  
    d1ham_dia[i] = new MATRIX(nelec,nelec); *d1ham_dia[i] = 0.0; 
    d1ham_adi[i] = new MATRIX(nelec,nelec); *d1ham_adi[i] = 0.0; 
  }
  for(i=0;i<nnucl*nnucl;i++){  d2ham_dia[i] = new MATRIX(nelec,nelec); *d2ham_dia[i] = 0.0;  }


  rep = 0; // default representation is diabatic  
  status_dia = 0;
  status_adi = 0;

}

Hamiltonian_Model::~Hamiltonian_Model(){
  int i;

//  cout<<"Derived Ham_mod. destructor\n";
 
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


}


void Hamiltonian_Model::set_params(vector<double>& params_){

  int num_params = 0;

  if(ham_indx==0 || ham_indx==2){ // SAC, ECWR
    num_params = 4;
  }
  else if(ham_indx==5){ // Rabi2
    num_params = 3;
  }
  else if(ham_indx==1 || ham_indx==3){  // DAC, Marcus
    num_params = 5;
  }
  else if(ham_indx==4){  // SEXCH
    num_params = 12;
  }

  params = vector<double>(num_params, 0.0);

  // Now copy input params:
  for(int i=0;i<params_.size() && i<num_params; i++){
    params[i] = params_[i];
  }

  // Since the parameters have changed - we need to recompute everything
  status_dia = 0;
  status_adi = 0;

}


void Hamiltonian_Model::compute_diabatic(){

  if(status_dia == 0){ // only compute this is the result is not up to date
  
    double e;

    //=============== 1D models ========================
    double x = q[0];

    if(ham_indx==0){  SAC_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);     }    // SAC potetnial
    else if(ham_indx==1){ DAC_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params); }    // DAC potential
    else if(ham_indx==2){ ECWR_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);   } // ECWR potential
    else if(ham_indx==3){ Marcus_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params); } // Marcus potential
    else if(ham_indx==4){ SEXCH_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);  } // SEXCH potential
    else if(ham_indx==5){ Rabi2_Ham(x, ham_dia,d1ham_dia[0],d2ham_dia[0], params);  } // Rabi2 potential


    // Set status flag 
    status_dia = 1;

  }//   status_dia == 0

}


void Hamiltonian_Model::compute_adiabatic(){
// This function computes adiabatic PESs and derivative couplings

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
  
}


}// namespace libhamiltonian_model
}// namespace libhamiltonian




