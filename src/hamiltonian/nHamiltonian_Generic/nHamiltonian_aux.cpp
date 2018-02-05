/*********************************************************************************
* Copyright (C) 2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian_aux.cpp
  \brief The file implements some auxiliary functions
    
*/


#include <stdlib.h>

#include "nHamiltonian.h"


/// liblibra namespace
namespace liblibra{

/// libhamiltonian namespace 
namespace libhamiltonian{

/// libhamiltonian_generic namespace 
namespace libhamiltonian_generic{

void check_mat_dimensions(CMATRIX* x, int nrows, int ncols){

  /// Check that the size of the external objects
  /// corresponds to the internal dimensions of the Hamiltonian
  if (x->n_rows!=nrows){  cout<<"Error in check_mat_dimensions: The number of rows of the matrix object "<<x->n_rows<<"\
    does not match the expected dimensionality of "<<nrows<<endl;
    exit(0);
  }

  if (x->n_cols!=ncols){  cout<<"Error in check_mat_dimensions: The number of cols of the matrix object "<<x->n_cols<<"\
    does not match the expected dimensionality of "<<ncols<<endl;
    exit(0);
  }

}


void check_vector_dimensions(vector<CMATRIX*> ptx, vector<CMATRIX>& x_, vector<int>& x_mem_status, int nnucl){
/**
  This is an auxiliary function to set a pointer to a matrix object to
  make reference to the existing object.
*/
  if(ptx.size()!=x_.size()){
    cout<<"Error in check_vector_dimensions: the size of the ptx ("<<ptx.size()<<") is not equal to the \
    size of x_ ("<<x_.size()<<")\n"; exit(0);
  }
  if(ptx.size()!=x_mem_status.size()){
    cout<<"Error in check_vector_dimensions: the size of the ptx ("<<ptx.size()<<") is not equal to the \
    size of x_mem_status ("<<x_mem_status.size()<<")\n"; exit(0);
  }
  if(ptx.size()!=nnucl){
    cout<<"Error in check_vector_dimensions: the size of the ptx ("<<ptx.size()<<") is not equal to the \
    nnucl in the input ("<<nnucl<<")\n"; exit(0);
  }

}


void set_X1_by_ref(CMATRIX* ptx, CMATRIX& x_, int& x_mem_status, int nrows, int ncols){
/**
  This is an auxiliary function to set a pointer to a matrix object to
  make reference to the existing object.
*/
  check_mat_dimensions(&x_,nrows,ncols);

  if(x_mem_status==0){ ptx = &x_; } /// Not allocated - not a problem
  else if(x_mem_status==1){ delete ptx; ptx = &x_;   }
  else if(x_mem_status==2){ ptx = &x_; }

  x_mem_status = 2; // Allocated externally

}


void set_X1_by_val(CMATRIX* ptx, CMATRIX& x_, int& x_mem_status, int nrows, int ncols){
/**
  This is an auxiliary function to set a pointer to a matrix object to
  make reference to the existing object.
*/

  check_mat_dimensions(&x_,nrows,ncols);

  /// Not allocated - not a problem - allocate
  if(x_mem_status==0){  ptx = new CMATRIX(nrows,ncols);  } 
  /// Allocated internally - gret! - just make sure it matches
  else if(x_mem_status==1){ check_mat_dimensions(ptx,nrows,ncols);  }
  /// Allocated externally - allocate memory and make the pointer
  /// to point to the internal location
  else if(x_mem_status==2){  ptx = new CMATRIX(nrows,ncols);    }

  *ptx = x_;
  x_mem_status = 1; // Allocated internally


}



void init_X1(CMATRIX* ptx, int& x_mem_status, int nrows, int ncols){
/**
  This is an auxiliary function that checks the memory allocation and allocates it if needed
*/

  if(x_mem_status==0){ // Not allocated
    ptx = new CMATRIX(nrows,ncols);
    x_mem_status = 1;
  }
  else if(x_mem_status==1){ // Allocated internally

    /// Check that the size of the internal storage corresponds to what we need
    if(ptx->n_rows==nrows && ptx->n_cols==ncols){ ;; }
    else{
      delete ptx; ptx = new CMATRIX(nrows, ncols);
    }
  }
  else if(x_mem_status==2){ // Allocated externally

    if(ptx->n_rows!=nrows){  cout<<"Error: The number of rows of the external object "<<ptx->n_rows<<"\
      does not match the expected dimensionality of "<<nrows<<endl;
      exit(0);
    }

    if(ptx->n_cols!=ncols){  cout<<"Error: The number of cols of the external object "<<ptx->n_cols<<"\
      does not match the expected dimensionality of "<<ncols<<endl;
      exit(0);
    }
  }
}



void set_X2_by_ref(vector<CMATRIX*> ptx, vector<CMATRIX>& x_, vector<int>& x_mem_status, int nrows, int ncols, int nnucl){
/**
  This is an auxiliary function to set a pointer to a matrix object to
  make reference to the existing object.
*/
  if(ptx.size()!=x_.size()){
    cout<<"Error in set_X2_by_ref: the size of the ptx ("<<ptx.size()<<") is not equal to the \
    size of x_ ("<<x_.size()<<")\n"; exit(0);
  }
  if(ptx.size()!=x_mem_status.size()){
    cout<<"Error in set_X2_by_ref: the size of the ptx ("<<ptx.size()<<") is not equal to the \
    size of x_mem_status ("<<x_mem_status.size()<<")\n"; exit(0);
  }
  if(ptx.size()!=nnucl){
    cout<<"Error in set_X2_by_ref: the size of the ptx ("<<ptx.size()<<") is not equal to the \
    nnucl in the input ("<<nnucl<<")\n"; exit(0);
  }

  // If we are here - we are good:
  for(int n=0;n<nnucl;n++){
    set_X1_by_ref(ptx[n], x_[n], x_mem_status[n], nrows, ncols);
  }

}


void set_X2_by_val(vector<CMATRIX*> ptx, vector<CMATRIX>& x_, vector<int>& x_mem_status, int nrows, int ncols, int nnucl){
/**
  This is an auxiliary function to set a pointer to a matrix object to
  make reference to the existing object.
*/
  if(ptx.size()!=x_.size()){
    cout<<"Error in set_X2_by_val: the size of the ptx ("<<ptx.size()<<") is not equal to the \
    size of x_ ("<<x_.size()<<")\n"; exit(0);
  }
  if(ptx.size()!=x_mem_status.size()){
    cout<<"Error in set_X2_by_val: the size of the ptx ("<<ptx.size()<<") is not equal to the \
    size of x_mem_status ("<<x_mem_status.size()<<")\n"; exit(0);
  }
  if(ptx.size()!=nnucl){
    cout<<"Error in set_X2_by_val: the size of the ptx ("<<ptx.size()<<") is not equal to the \
    nnucl in the input ("<<nnucl<<")\n"; exit(0);
  }

  // If we are here - we are good:
  for(int n=0;n<nnucl;n++){
    set_X1_by_val(ptx[n], x_[n], x_mem_status[n], nrows, ncols);
  }

}


void init_X2(vector<CMATRIX*> ptx, vector<int>& x_mem_status, int nrows, int ncols, int nnucl){


  if(ptx.size()!=x_mem_status.size()){
    cout<<"Error in init_X2: the size of the ptx ("<<ptx.size()<<") is not equal to the \
    size of x_mem_status ("<<x_mem_status.size()<<")\n"; exit(0);
  }
  if(ptx.size()!=nnucl){
    cout<<"Error in init_X2: the size of the ptx ("<<ptx.size()<<") is not equal to the \
    nnucl in the input ("<<nnucl<<")\n"; exit(0);
  }


  // If we are here - we are good:
  for(int n=0;n<nnucl;n++){
    init_X1(ptx[n], x_mem_status[n], nrows, ncols);
  }

}






}// namespace libhamiltonian_generic
}// namespace libhamiltonian
}// liblibra

