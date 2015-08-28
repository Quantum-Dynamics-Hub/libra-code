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

#include "Electronic.h"

/****************************************************************************
  This file contains following functions:

  void propagate_electronic(double dt,Electronic& el,Hamiltonian& pot)

****************************************************************************/

namespace libdyn{
namespace libelectronic{


using libmmath::liboperators::rotate;


void Electronic::propagate_electronic(double dt,Hamiltonian* ham){
  libdyn::libelectronic::propagate_electronic(dt,this,ham);
}

void Electronic::propagate_electronic(double dt,Hamiltonian& ham){
  libdyn::libelectronic::propagate_electronic(dt,this,&ham);
}


void propagate_electronic(double dt,Electronic* el,Hamiltonian* ham){
// This version is based on the Hamiltonian formulation of TD-SE
// This propagator is good for general Hamiltonian - diabatic of adiabatic
//
// iL = iL_qp + iL_qq + iL_pp
//
// iL_qp = sum_i,j {}
//
// ...


  int i,j;

  double dt_half = 0.5*dt;

  //------------- Phase evolution (adiabatic) ---------------- 
  // exp(iL_qp * dt/2)
  for(i=0;i<el->nstates;i++){
    for(j=0;j<el->nstates;j++){

      rotate(el->p[j],el->q[i], dt_half*ham->Hvib(i,j).real());

    }// for j
  }// for i

  //------------- Population transfer (adiabatic) ----------------
  // exp((iL_qq + iL_pp) * dt/2)
  for(i=0;i<el->nstates;i++){
    for(j=i+1;j<el->nstates;j++){

      rotate(el->q[j],el->q[i], dt_half*ham->Hvib(i,j).imag());
      rotate(el->p[j],el->p[i], dt_half*ham->Hvib(i,j).imag());

    }// for j
  }// for i

  // exp((iL_qq + iL_pp) * dt/2)
  for(i=el->nstates-1;i>=0;i--){
    for(j=el->nstates-1;j>i;j--){


      rotate(el->q[j],el->q[i], dt_half*ham->Hvib(i,j).imag());
      rotate(el->p[j],el->p[i], dt_half*ham->Hvib(i,j).imag());

    }// for j
  }// for i


  //------------- Phase evolution (adiabatic) ----------------
  // exp(iL_qp * dt/2)
  for(i=el->nstates-1;i>=0;i--){
    for(j=el->nstates-1;j>=0;j--){

      rotate(el->p[j],el->q[i], dt_half*ham->Hvib(i,j).real());

    }// for j
  }// for i


}// propagate_electronic


}// namespace libelectronic
}// namespace libdyn

