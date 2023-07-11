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
  \file nHamiltonian_aux.cpp
  \brief The file implements some auxiliary functions
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <stdlib.h>
#endif

#include "nHamiltonian.h"


/// liblibra namespace
namespace liblibra{

namespace bp = boost::python;

/// libnhamiltonian namespace 
namespace libnhamiltonian{


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



void nHamiltonian::check_cmatrix(bp::object obj, std::string matrix_name, int nrows, int ncols){
/**
  This function checks whether the dimensions of the matrix named <matrix_name>
  extracted from a Python object <obj> agree with the expected number of rows <nrows>
  and number of columns <ncols>
*/

    CMATRIX& tmp = extract<CMATRIX&>(obj.attr(matrix_name.c_str()));

    if(tmp.n_rows != nrows){
      cout<<"ERROR in nHamiltonian...\n"; 
      cout<<"The number of rows in the extracted "<<matrix_name<<" object ("<<tmp.n_rows<<")"
          <<" is not equal to the expected number of rows = "<<nrows<<"\nExiting...\n";
      exit(0);
    }
    if(tmp.n_cols != ncols){
      cout<<"ERROR in nHamiltonian...\n"; 
      cout<<"The number of rows in the extracted "<<matrix_name<<" object ("<<tmp.n_cols<<")"
          <<" is not equal to the expected number of columns = "<<ncols<<"\nExiting...\n";
      exit(0);
    }

}


void nHamiltonian::check_cmatrix_list(bp::object obj, std::string matrix_name, int nrows, int ncols, int length){
/**
  This function checks whether the dimensions of the vector of matrices named <matrix_name>
  extracted from a Python object <obj> agree with the expected number of rows <nrows>
  and number of columns <ncols> for each element and whether the number of elements of the list is 
  equal to the expected value of <length>
*/

    vector<CMATRIX>& tmp = extract<vector<CMATRIX>&>(obj.attr(matrix_name.c_str()));

    if(tmp.size() != length){
      cout<<"ERROR in nHamiltonian...\n"; 
      cout<<"The number of elements in the extracted "<<matrix_name<<" object ("<<tmp.size()<<")"
          <<"is not equal to the expected number of elements = "<<length<<"\nExiting...\n";
      exit(0);
    }

    // For speed, check the dimensions of only the first element
    if(tmp.size()>0){

      if(tmp[0].n_rows != nrows){
        cout<<"ERROR in nHamiltonian...\n"; 
        cout<<"The number of rows in the extracted "<<matrix_name<<"[0] object ("<<tmp[0].n_rows<<")"
            <<" is not equal to the expected number of rows = "<<nrows<<"\nExiting...\n";
        exit(0);
      }
      if(tmp[0].n_cols != ncols){
        cout<<"ERROR in nHamiltonian...\n"; 
        cout<<"The number of rows in the extracted "<<matrix_name<<"[0] object ("<<tmp[0].n_cols<<")"
            <<" is not equal to the expected number of columns = "<<ncols<<"\nExiting...\n";
        exit(0);
      }
    }


}




void nHamiltonian::add_branches(int target_level, vector<nHamiltonian*>& res){
/**
*/ 

  if(level==target_level){  res.push_back(this);  }
  else{

    int sz = children.size();
    for(int i; i<sz; i++){   children[i]->add_branches(target_level, res);   }

  }

}

vector<nHamiltonian*> nHamiltonian::get_branches(int target_level){
/**
  An auxiliary function to traverse the entire tree to return the pointers to the 
  sub-Hamiltonian of the given level.

  The tree traversal is done according to: the leftmost-first principle

*/ 

  vector<nHamiltonian*> res;

  add_branches(target_level, res);

  return res;

}



}// namespace libnhamiltonian
}// liblibra

