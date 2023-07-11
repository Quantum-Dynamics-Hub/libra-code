/*********************************************************************************
* Copyright (C) 2017-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian_compute_nac.cpp
  \brief The file implements computations of the Nonadiabatic couplings (NACs) and
  vibronic Hamiltonians
    
*/


#if defined(USING_PCH)
#include "../pch.h"
#else
#include <stdlib.h>
#endif 

#include "nHamiltonian.h"
#include "../math_meigen/libmeigen.h"


/// liblibra namespace
namespace liblibra{

/// libnhamiltonian namespace 
namespace libnhamiltonian{


using namespace liblinalg;
using namespace libmeigen;




void nHamiltonian::compute_nac_dia(MATRIX& p, const MATRIX& invM){
/**  
  Function to compute nonadiabatic coupling (nac) in the diabatic representation
  This coupiling is also known as the time-derivative coupling and it is a scalar

  Assumed dimensions: 
  p - Ndof x 1 matrix of nuclear momenta
  invM - Ndof x 1 matrix of inverse masses for all the nuclear DOFs

  This is going to be <i|d/dt|j> = sum_n { p[n]/mass[n] * dc1_adi[n] }
*/

  if(nac_dia_mem_status==0){ cout<<"Error in compute_nac_dia(): the memory is not allocated for \
  nac_dia but is needed for the calculations \n"; exit(0); }

  nac_dia->set(-1, std::complex<double>(0.0, 0.0));

  for(int n=0;n<nnucl;n++){ 

    double v_n = p.get(n,0) * invM.get(n,0);

    if(dc1_dia_mem_status[n]==0){ cout<<"Error in compute_nac_dia(): the derivatives couplings matrix in the adiabatic \
    basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }


    for(int i=0;i<ndia;i++){
      for(int j=0;j<ndia;j++){
        nac_dia->add(i,j, v_n * dc1_dia[n]->get(i,j));
      }// for j
    }// for i

  }// for n

}


void nHamiltonian::compute_nac_dia(MATRIX& p, const MATRIX& invM, vector<int>& id_){
/**
  See the description of the compute_nac_dia(MATRIX& p, MATRIX& invM) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   compute_nac_dia(p, invM);    }
    else{ cout<<"ERROR in compute_nac_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    children[id_[1]]->compute_nac_dia(p, invM, next);
  }
}


void nHamiltonian::compute_nac_dia(MATRIX& p, const MATRIX& invM, int lvl, int split){
/**
  \brief Function to compute nonadiabatic coupling (nac) in the diabatic representation
   This coupiling is also known as the time-derivative coupling and it is a scalar.


  \param[in] p [ndof x 1] or [ndof x ntraj] matrix of nuclear momenta
  \param[in] invM [ndof x 1] matrix of inverse nuclear masses
  \param[in] lvl The level in the hierarchy of Hamiltonians at which we will perform the transformation
  \param[in] split The flag tells the function to perform the calculations 
  column by column, with column being used by a children Hamiltonian (next level 
  relative to the input parameter lvl). This is only allowed if the number of columns in p
  is equal to the number of children of the present-level Hamiltonian.

  There are 4 possible use cases:
  a) p is a ndof x 1 matrix (one trajectory) and is meant to be used by the current
  Hamiltonian - this is already implemented above
  b) p is a ndof x 1 matrix (one trajectory) but is meant to be used by the Hamiltonian
  of a different level - call the appropriate level children.
  c) p is a ndof x ntraj matrix containing momenta for multiple trajectories and they are
  meant to be used by the children sub-Hamiltonians of the current level Hamiltonian 
  one by one (this is activated by the split=1 option)
  d) Same as "c", but relative to a different level of the Hamiltonians hierarchy.

*/
  int i;

  if(lvl==level){  // Case "a" 
    if(split==0){  // Transform all the columns using the present level Hamiltonian
      compute_nac_dia(p, invM);
    }
    else if(split==1){  // Case "c"
      // Check whether we have enough sub-Hamiltonians
      if(children.size()!=p.n_cols){
        cout<<"ERROR in void nHamiltonian::compute_nac_dia(const MATRIX& p, const MATRIX& invM, int lvl, int split):\n";
        cout<<"The number of columns of the p ("<<p.n_cols<<")";
        cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
        cout<<"Exiting...\n";
        exit(0);
      }

      MATRIX p_tmp(p.n_rows, 1);
      vector<int> stenc_p(p.n_rows, 0);
      vector<int> stenc_col(1, 0);
      for(i=0;i<p.n_rows;i++){ stenc_p[i] = i;}


      for(i=0;i<children.size();i++){
        stenc_col[0] = i;

        pop_submatrix(p, p_tmp, stenc_p, stenc_col);
        children[i]->compute_nac_dia(p_tmp, invM);
 
      }// for all children

    }// split==1
    else{
      cout<<"ERROR in void nHamiltonian::compute_nac_dia(const MATRIX& p, const MATRIX& invM, int lvl, int split):\n";
      cout<<"The parameters split = "<<split<<" is not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }
  }// lvl == level

  else if(lvl>level){  // Cases "b" or "d"
  
    for(int i=0;i<children.size();i++){
      children[i]->compute_nac_dia(p, invM, lvl, split);
    }

  }// lvl >level

  else{
    cout<<"WARNING in nHamiltonian::compute_nac_dia\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }
}



void nHamiltonian::compute_nac_adi(double dt, int method){
/***
  One of the ways to update scalar NAC matrix, using the time-overlap info

*/ 

  if(time_overlap_adi_mem_status==0){ cout<<"Error in compute_nac_adi(): the memory is not allocated for \
  time_overlap_adi but is needed for the calculations \n"; exit(0); }

  if(nac_adi_mem_status==0){ cout<<"Error in compute_nac_adi(): the memory is not allocated for \
  nac_adi but is needed for the calculations \n"; exit(0); }

  nac_adi->set(-1, std::complex<double>(0.0, 0.0));


  // HST formula
  if(method==0){  

 //   nac_adi->
  }

}


void nHamiltonian::compute_nac_adi(MATRIX& p, const MATRIX& invM){
/**  
  Function to compute nonadiabatic coupling (nac) in the adiabatic representation
  This coupiling is also known as the time-derivative coupling and it is a scalar

  This is going to be <i|d/dt|j> = sum_n { p[n]/mass[n] * dc1_adi[n] }
*/

  if(nac_adi_mem_status==0){ cout<<"Error in compute_nac_adi(): the memory is not allocated for \
  nac_adi but is needed for the calculations \n"; exit(0); }

  nac_adi->set(-1, std::complex<double>(0.0, 0.0));

  for(int n=0;n<nnucl;n++){ 

    if(dc1_adi_mem_status[n]==0){ cout<<"Error in compute_nac_dia(): the derivatives couplings matrix in the adiabatic \
    basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }

    double v_n = p.get(n,0) * invM.get(n,0);

    for(int i=0;i<nadi;i++){
      for(int j=0;j<nadi;j++){
        nac_adi->add(i,j, v_n * dc1_adi[n]->get(i,j));
      }// for j
    }// for i

  }// for n


}


void nHamiltonian::compute_nac_adi(MATRIX& p, const MATRIX& invM, vector<int>& id_){
/**
  See the description of the compute_nac_adi(MATRIX& p, MATRIX& invM) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   compute_nac_adi(p, invM);    }
    else{ cout<<"ERROR in compute_nac_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    children[id_[1]]->compute_nac_adi(p, invM, next);
  }
}


void nHamiltonian::compute_nac_adi(MATRIX& p, const MATRIX& invM, int lvl, int split){
/**
  \brief Function to compute nonadiabatic coupling (nac) in the adiabatic representation
   This coupiling is also known as the time-derivative coupling and it is a scalar.

  \param[in] p [ndof x 1] or [ndof x ntraj] matrix of nuclear momenta
  \param[in] invM [ndof x 1] matrix of inverse nuclear masses
  \param[in] lvl The level in the hierarchy of Hamiltonians at which we will perform the transformation
  \param[in] split The flag tells the function to perform the calculations 
  column by column, with column being used by a children Hamiltonian (next level 
  relative to the input parameter lvl). This is only allowed if the number of columns in p
  is equal to the number of children of the present-level Hamiltonian.

  There are 4 possible use cases:
  a) p is a ndof x 1 matrix (one trajectory) and is meant to be used by the current
  Hamiltonian - this is already implemented above
  b) p is a ndof x 1 matrix (one trajectory) but is meant to be used by the Hamiltonian
  of a different level - call the appropriate level children.
  c) p is a ndof x ntraj matrix containing momenta for multiple trajectories and they are
  meant to be used by the children sub-Hamiltonians of the current level Hamiltonian 
  one by one (this is activated by the split=1 option)
  d) Same as "c", but relative to a different level of the Hamiltonians hierarchy.

*/
  int i;

  if(lvl==level){  // Case "a" 
    if(split==0){  // Transform all the columns using the present level Hamiltonian
      compute_nac_adi(p, invM);
    }
    else if(split==1){  // Case "c"
      // Check whether we have enough sub-Hamiltonians
      if(children.size()!=p.n_cols){
        cout<<"ERROR in void nHamiltonian::compute_nac_adi(const MATRIX& p, const MATRIX& invM, int lvl, int split):\n";
        cout<<"The number of columns of the p ("<<p.n_cols<<")";
        cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
        cout<<"Exiting...\n";
        exit(0);
      }

      MATRIX p_tmp(p.n_rows, 1);
      vector<int> stenc_p(p.n_rows, 0);
      vector<int> stenc_col(1, 0);
      for(i=0;i<p.n_rows;i++){ stenc_p[i] = i;}


      for(i=0;i<children.size();i++){
        stenc_col[0] = i;

        pop_submatrix(p, p_tmp, stenc_p, stenc_col);
        children[i]->compute_nac_adi(p_tmp, invM);
 
      }// for all children

    }// split==1
    else{
      cout<<"ERROR in void nHamiltonian::compute_nac_adi(const MATRIX& p, const MATRIX& invM, int lvl, int split):\n";
      cout<<"The parameters split = "<<split<<" is not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }
  }// lvl == level

  else if(lvl>level){  // Cases "b" or "d"
  
    for(int i=0;i<children.size();i++){
      children[i]->compute_nac_adi(p, invM, lvl, split);
    }

  }// lvl >level

  else{
    cout<<"WARNING in nHamiltonian::compute_nac_adi\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }
}






void nHamiltonian::compute_hvib_dia(){
/**  
  Function to compute vibronic Hamiltonian in the diabatic representation
  Assume the NAC is updated!

  This is going to be Hvib = H_el - i*hbar * NAC
*/
  complex<double> ihbar(0.0, 1.0); // working in atomic units

  if(nac_dia_mem_status==0){ cout<<"Error in compute_hvib_dia(): the memory is not allocated for \
  nac_dia but is needed for the calculations \n"; exit(0); }

  if(ham_dia_mem_status==0){ cout<<"Error in compute_hvib_dia(): the memory is not allocated for \
  ham_dia but is needed for the calculations \n"; exit(0); }

  if(hvib_dia_mem_status==0){ cout<<"Error in compute_hvib_dia(): the memory is not allocated for \
  hvib_dia but is needed for the calculations \n"; exit(0); }


  for(int i=0;i<ndia;i++){
    for(int j=0;j<ndia;j++){
      hvib_dia->set(i,j, ham_dia->get(i,j) - ihbar*nac_dia->get(i,j) );
    }// for j
  }// for i

}

void nHamiltonian::compute_hvib_dia(vector<int>& id_){
/**
  See the description of the compute_hvib_dia() function
*/
  if(id_.size()==1){
    if(id_[0]==id){   compute_hvib_dia();    }
    else{ cout<<"ERROR in compute_hvib_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    children[id_[1]]->compute_hvib_dia(next);
  }
}



void nHamiltonian::compute_hvib_dia(int lvl){
/**
  Re-evaluates the diabatic vibronic Hamiltonians for all 
  children Hamiltonians of a given level lvl.
  The top-most level is 0. 

  \param[in] lvl The level of the sub-Hamiltonians to be re-evaluated
    
*/

  if(level==lvl){   compute_hvib_dia();   }// level == lvl
  else if(lvl>level){
  
    for(int i=0;i<children.size();i++){   children[i]->compute_hvib_dia(lvl);   }

  }// lvl >level
  else{
    cout<<"WARNING in nHamiltonian::compute_hvib_dia\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }
  
}





void nHamiltonian::compute_hvib_adi(){
/**  
  Function to compute vibronic Hamiltonian in the adiabatic representation
  Assume the NAC is updated!

  This is going to be Hvib = H_el - i*hbar * NAC
*/
  complex<double> ihbar(0.0, 1.0); // working in atomic units

  if(nac_adi_mem_status==0){ cout<<"Error in compute_hvib_adi(): the memory is not allocated for \
  nac_adi but is needed for the calculations \n"; exit(0); }

  if(ham_adi_mem_status==0){ cout<<"Error in compute_hvib_adi(): the memory is not allocated for \
  ham_adi but is needed for the calculations \n"; exit(0); }

  if(hvib_adi_mem_status==0){ cout<<"Error in compute_hvib_adi(): the memory is not allocated for \
  hvib_adi but is needed for the calculations \n"; exit(0); }


  for(int i=0;i<nadi;i++){
    for(int j=0;j<nadi;j++){
      hvib_adi->set(i,j, ham_adi->get(i,j) - ihbar*nac_adi->get(i,j) );
    }// for j
  }// for i

}

void nHamiltonian::compute_hvib_adi(vector<int>& id_){
/**
  See the description of the compute_hvib_adi() function
*/
  if(id_.size()==1){
    if(id_[0]==id){   compute_hvib_adi();    }
    else{ cout<<"ERROR in compute_hvib_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    children[id_[1]]->compute_hvib_adi(next);
  }
}


void nHamiltonian::compute_hvib_adi(int lvl){
/**
  Re-evaluates the adiabatic vibronic Hamiltonians for all 
  children Hamiltonians of a given level lvl.
  The top-most level is 0. 

  \param[in] lvl The level of the sub-Hamiltonians to be re-evaluated
    
*/

  if(level==lvl){   compute_hvib_adi();   }// level == lvl
  else if(lvl>level){
  
    for(int i=0;i<children.size();i++){   children[i]->compute_hvib_adi(lvl);   }

  }// lvl >level
  else{
    cout<<"WARNING in nHamiltonian::compute_hvib_adi\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }
  
}



}// namespace libnhamiltonian
}// liblibra

