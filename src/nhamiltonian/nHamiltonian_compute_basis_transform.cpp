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
  \file nHamiltonian_compute_basis_transform.cpp
  \brief The file implements various versions of functions to transform between diabatic
  and adiabatic bases and vice versa
    
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





void nHamiltonian::ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi){
/**

  |PSI> = |psi_dia> * C_dia = |psi_adi> * C_adi

  |psi_adi> = |psi_dia> * U, so:

  C_dia = U * C_adi

  U^H * S * U = I ==> U^-1 = U^H() * S, and C_adi = U^H() * S * C_dia

*/

  if(ovlp_dia_mem_status==0){ cout<<"Error in ampl_dia2adi(): the overlap matrix in the diabatic basis is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(basis_transform_mem_status==0){ cout<<"Error in ampl_dia2adi(): the transformation basis matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

//  if(ampl_dia_mem_status==0){ cout<<"Error in ampl_dia2adi(): the amplitudes of the diabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }

//  if(ampl_adi_mem_status==0){ cout<<"Error in ampl_dia2adi(): the amplitudes of the adiabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }



  ampl_adi = (*basis_transform).H() * (*ovlp_dia) * ampl_dia; 

}


void nHamiltonian::ampl_dia2adi(CMATRIX* ampl_dia, CMATRIX* ampl_adi){
/**

  |PSI> = |psi_dia> * C_dia

  C_dia = U * C_adi

  U^H * S * U = I ==> U^-1 = U^H() * S, and C_adi = U^H() * S * C_dia

*/

  if(ovlp_dia_mem_status==0){ cout<<"Error in ampl_dia2adi(): the overlap matrix in the diabatic basis is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(basis_transform_mem_status==0){ cout<<"Error in ampl_dia2adi(): the transformation basis matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

//  if(ampl_dia_mem_status==0){ cout<<"Error in ampl_dia2adi(): the amplitudes of the diabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }

//  if(ampl_adi_mem_status==0){ cout<<"Error in ampl_dia2adi(): the amplitudes of the adiabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }



  (*ampl_adi) = (*basis_transform).H() * (*ovlp_dia) * (*ampl_dia); 

}


void nHamiltonian::ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi, vector<int>& id_){
/**
  See the description of the ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return ampl_dia2adi(ampl_dia, ampl_adi);    }
    else{ cout<<"ERROR in ampl_dia2adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->ampl_dia2adi(ampl_dia, ampl_adi, next);
  }
}

void nHamiltonian::ampl_dia2adi(CMATRIX* ampl_dia, CMATRIX* ampl_adi, vector<int>& id_){
/**
  See the description of the ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return ampl_dia2adi(ampl_dia, ampl_adi);    }
    else{ cout<<"ERROR in ampl_dia2adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->ampl_dia2adi(ampl_dia, ampl_adi, next);
  }
}




void nHamiltonian::ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split){
/**
  \brief Transforms diabatic amplitudes to the adiabatic:  ampl_dia -> ampl_adi
  \param[in] ampl_dia A [ndia x 1] or [ndia x ntraj] matrix of diabatic amplitudes for
  one of many trajectories
  \param[in, out] ampl_adi A [nadi x 1] or [nadi x ntraj] matrix of adiabatic amplitudes for
  one of many trajectories
  \param[in] lvl The level in the hierarchy of Hamiltonians at which we will perform the transformation
  \param[in] split The flag telling the function to perform the transformations of the matrix
  column by column, with each transformation being done by a children Hamiltonian (next level 
  relative to the input parameter lvl). This is only allowed if the number of columns in both
  ampl_dia and ampl_adi are equal to the number of children of the present-level Hamiltonian.

  There are 4 possible use cases:
  a) Both ampl_dia and ampl_adi are N x 1 matrices (one trajectory) and are meant to 
  be handled by the current Hamiltonian - this is already implemented above
  b) Both ampl_dia and ampl_adi are N x Ntraj matrices (Ntraj trajectories) and are still
  meant to be transformed by the same Hamiltonian - the current one. This case fits the 
  use case "a" due to the underlying matrix transformation in the called function.
  c, d) These use cases are analogous to "a" and "b", except for the fact that we need to
  call a lower-level Hamiltonian.
  Any of the above case can be combined with the "split" option to do the calculations 
  with the sub-Hamiltonians.

*/
  int i;

  if(lvl==level){  // Cases "a" and "b"
    if(split==0){  // Transform all the columns using the present level Hamiltonian
      ampl_dia2adi(ampl_dia, ampl_adi);
    }
    else if(split==1){
      // Check whether we have enough sub-Hamiltonians
      if(children.size()!=ampl_dia.n_cols){
        cout<<"ERROR in void nHamiltonian::ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
        cout<<"The number of columns of the ampl_dia ("<<ampl_dia.n_cols<<")";
        cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
        cout<<"Exiting...\n";
        exit(0);
      }
      if(children.size()!=ampl_adi.n_cols){
        cout<<"ERROR in void nHamiltonian::ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
        cout<<"The number of columns of the ampl_adi ("<<ampl_adi.n_cols<<")";
        cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
        cout<<"Exiting...\n";
        exit(0);
      }

      CMATRIX dia_tmp(ampl_dia.n_rows, 1);
      CMATRIX adi_tmp(ampl_adi.n_rows, 1);

      vector<int> stenc_dia1(ampl_dia.n_rows, 0);
      vector<int> stenc_adi1(ampl_dia.n_rows, 0);      
      vector<int> stenc_col(1, 0);

      for(i=0;i<ampl_dia.n_rows;i++){ stenc_dia1[i] = i;}
      for(i=0;i<ampl_adi.n_rows;i++){ stenc_adi1[i] = i;}

      for(i=0;i<children.size();i++){
        stenc_col[0] = i;

        pop_submatrix(ampl_dia, dia_tmp, stenc_dia1, stenc_col);

        children[i]->ampl_dia2adi(dia_tmp, adi_tmp);

        push_submatrix(ampl_adi, adi_tmp, stenc_adi1, stenc_col);
 
      }// for all children

    }// split==1
    else{
      cout<<"ERROR in void nHamiltonian::ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
      cout<<"The parameters split = "<<split<<" is not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }
  }// lvl == level

  else if(lvl>level){  // Cases "c" and "d"
  
    for(int i=0;i<children.size();i++){
      children[i]->ampl_dia2adi(ampl_dia, ampl_adi, lvl, split);
    }

  }// lvl >level

  else{
    cout<<"WARNING in nHamiltonian::ampl_dia2adi\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }
}


void nHamiltonian::ampl_dia2adi(CMATRIX* ampl_dia, CMATRIX* ampl_adi, int lvl, int split){
/**
  \brief Transforms diabatic amplitudes to the adiabatic:  ampl_dia -> ampl_adi
  \param[in] ampl_dia A [ndia x 1] or [ndia x ntraj] matrix of diabatic amplitudes for
  one of many trajectories
  \param[in, out] ampl_adi A [nadi x 1] or [nadi x ntraj] matrix of adiabatic amplitudes for
  one of many trajectories
  \param[in] lvl The level in the hierarchy of Hamiltonians at which we will perform the transformation
  \param[in] split The flag telling the function to perform the transformations of the matrix
  column by column, with each transformation being done by a children Hamiltonian (next level 
  relative to the input parameter lvl). This is only allowed if the number of columns in both
  ampl_dia and ampl_adi are equal to the number of children of the present-level Hamiltonian.

  There are 4 possible use cases:
  a) Both ampl_dia and ampl_adi are N x 1 matrices (one trajectory) and are meant to 
  be handled by the current Hamiltonian - this is already implemented above
  b) Both ampl_dia and ampl_adi are N x Ntraj matrices (Ntraj trajectories) and are still
  meant to be transformed by the same Hamiltonian - the current one. This case fits the 
  use case "a" due to the underlying matrix transformation in the called function.
  c, d) These use cases are analogous to "a" and "b", except for the fact that we need to
  call a lower-level Hamiltonian.
  Any of the above case can be combined with the "split" option to do the calculations 
  with the sub-Hamiltonians.

*/
  int i;

  if(lvl==level){  // Cases "a" and "b"
    if(split==0){  // Transform all the columns using the present level Hamiltonian
      ampl_dia2adi(ampl_dia, ampl_adi);
    }
    else if(split==1){
      // Check whether we have enough sub-Hamiltonians
      if(children.size()!=ampl_dia->n_cols){
        cout<<"ERROR in void nHamiltonian::ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
        cout<<"The number of columns of the ampl_dia ("<<ampl_dia->n_cols<<")";
        cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
        cout<<"Exiting...\n";
        exit(0);
      }
      if(children.size()!=ampl_adi->n_cols){
        cout<<"ERROR in void nHamiltonian::ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
        cout<<"The number of columns of the ampl_adi ("<<ampl_adi->n_cols<<")";
        cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
        cout<<"Exiting...\n";
        exit(0);
      }

      CMATRIX dia_tmp(ampl_dia->n_rows, 1);
      CMATRIX adi_tmp(ampl_adi->n_rows, 1);

      vector<int> stenc_dia1(ampl_dia->n_rows, 0);
      vector<int> stenc_adi1(ampl_dia->n_rows, 0);      
      vector<int> stenc_col(1, 0);

      for(i=0;i<ampl_dia->n_rows;i++){ stenc_dia1[i] = i;}
      for(i=0;i<ampl_adi->n_rows;i++){ stenc_adi1[i] = i;}

      for(i=0;i<children.size();i++){
        stenc_col[0] = i;

        pop_submatrix(ampl_dia, &dia_tmp, stenc_dia1, stenc_col);

        children[i]->ampl_dia2adi(dia_tmp, adi_tmp);

        push_submatrix(ampl_adi, &adi_tmp, stenc_adi1, stenc_col);
 
      }// for all children

    }// split==1
    else{
      cout<<"ERROR in void nHamiltonian::ampl_dia2adi(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
      cout<<"The parameters split = "<<split<<" is not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }
  }// lvl == level

  else if(lvl>level){  // Cases "c" and "d"
  
    for(int i=0;i<children.size();i++){
      children[i]->ampl_dia2adi(ampl_dia, ampl_adi, lvl, split);
    }

  }// lvl >level

  else{
    cout<<"WARNING in nHamiltonian::ampl_dia2adi\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }
}



void nHamiltonian::ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi){
/**

  |PSI> = |psi_dia> * C_dia

  C_dia = U * C_adi


*/

  if(basis_transform_mem_status==0){ cout<<"Error in ampl_dia2adi(): the transformation basis matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

//  if(ampl_dia_mem_status==0){ cout<<"Error in ampl_dia2adi(): the amplitudes of the diabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }

//  if(ampl_adi_mem_status==0){ cout<<"Error in ampl_dia2adi(): the amplitudes of the adiabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }


  ampl_dia = (*basis_transform) * ampl_adi; 

}



void nHamiltonian::ampl_adi2dia(CMATRIX* ampl_dia, CMATRIX* ampl_adi){
/**

  |PSI> = |psi_dia> * C_dia

  C_dia = U * C_adi


*/

  if(basis_transform_mem_status==0){ cout<<"Error in ampl_dia2adi(): the transformation basis matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

//  if(ampl_dia_mem_status==0){ cout<<"Error in ampl_dia2adi(): the amplitudes of the diabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }

//  if(ampl_adi_mem_status==0){ cout<<"Error in ampl_dia2adi(): the amplitudes of the adiabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }


  (*ampl_dia) = (*basis_transform) * (*ampl_adi); 

}



void nHamiltonian::ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi, vector<int>& id_){
/**
  See the description of the ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return ampl_adi2dia(ampl_dia, ampl_adi);    }
    else{ cout<<"ERROR in ampl_adi2dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->ampl_adi2dia(ampl_dia, ampl_adi, next);
  }
}


void nHamiltonian::ampl_adi2dia(CMATRIX* ampl_dia, CMATRIX* ampl_adi, vector<int>& id_){
/**
  See the description of the ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return ampl_adi2dia(ampl_dia, ampl_adi);    }
    else{ cout<<"ERROR in ampl_adi2dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->ampl_adi2dia(ampl_dia, ampl_adi, next);
  }
}




void nHamiltonian::ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split){
/**
  \brief Transforms adiabatic amplitudes to the diabatic:  ampl_adi -> ampl_dia
  \param[in, out] ampl_dia A [ndia x 1] or [ndia x ntraj] matrix of diabatic amplitudes for
  one of many trajectories
  \param[in] ampl_adi A [nadi x 1] or [nadi x ntraj] matrix of adiabatic amplitudes for
  one of many trajectories
  \param[in] lvl The level in the hierarchy of Hamiltonians at which we will perform the transformation
  \param[in] split The flag telling the function to perform the transformations of the matrix
  column by column, with each transformation being done by a children Hamiltonian (next level 
  relative to the input parameter lvl). This is only allowed if the number of columns in both
  ampl_dia and ampl_adi are equal to the number of children of the present-level Hamiltonian.

  There are 4 possible use cases:
  a) Both ampl_dia and ampl_adi are N x 1 matrices (one trajectory) and are meant to 
  be handled by the current Hamiltonian - this is already implemented above
  b) Both ampl_dia and ampl_adi are N x Ntraj matrices (Ntraj trajectories) and are still
  meant to be transformed by the same Hamiltonian - the current one. This case fits the 
  use case "a" due to the underlying matrix transformation in the called function.
  c, d) These use cases are analogous to "a" and "b", except for the fact that we need to
  call a lower-level Hamiltonian.
  Any of the above case can be combined with the "split" option to do the calculations 
  with the sub-Hamiltonians.

*/
  int i;

  if(lvl==level){  // Cases "a" and "b"
    if(split==0){  // Transform all the columns using the present level Hamiltonian
      ampl_adi2dia(ampl_dia, ampl_adi);
    }
    else if(split==1){
      // Check whether we have enough sub-Hamiltonians
      if(children.size()!=ampl_dia.n_cols){
        cout<<"ERROR in void nHamiltonian::ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
        cout<<"The number of columns of the ampl_dia ("<<ampl_dia.n_cols<<")";
        cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
        cout<<"Exiting...\n";
        exit(0);
      }
      if(children.size()!=ampl_adi.n_cols){
        cout<<"ERROR in void nHamiltonian::ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
        cout<<"The number of columns of the ampl_adi ("<<ampl_adi.n_cols<<")";
        cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
        cout<<"Exiting...\n";
        exit(0);
      }

      CMATRIX dia_tmp(ampl_dia.n_rows, 1);
      CMATRIX adi_tmp(ampl_adi.n_rows, 1);

      vector<int> stenc_dia1(ampl_dia.n_rows, 0);
      vector<int> stenc_adi1(ampl_dia.n_rows, 0);      
      vector<int> stenc_col(1, 0);

      for(i=0;i<ampl_dia.n_rows;i++){ stenc_dia1[i] = i;}
      for(i=0;i<ampl_adi.n_rows;i++){ stenc_adi1[i] = i;}

      for(i=0;i<children.size();i++){
        stenc_col[0] = i;

        pop_submatrix(ampl_adi, adi_tmp, stenc_adi1, stenc_col);

        children[i]->ampl_adi2dia(dia_tmp, adi_tmp);

        push_submatrix(ampl_dia, dia_tmp, stenc_dia1, stenc_col);
 
      }// for all children

    }// split==1
    else{
      cout<<"ERROR in void nHamiltonian::ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
      cout<<"The parameters split = "<<split<<" is not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }
  }// lvl == level

  else if(lvl>level){  // Cases "c" and "d"
  
    for(int i=0;i<children.size();i++){
      children[i]->ampl_adi2dia(ampl_dia, ampl_adi, lvl, split);
    }

  }// lvl >level

  else{
    cout<<"WARNING in nHamiltonian::ampl_dia2adi\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }
}



void nHamiltonian::ampl_adi2dia(CMATRIX* ampl_dia, CMATRIX* ampl_adi, int lvl, int split){
/**
  \brief Transforms adiabatic amplitudes to the diabatic:  ampl_adi -> ampl_dia
  \param[in, out] ampl_dia A [ndia x 1] or [ndia x ntraj] matrix of diabatic amplitudes for
  one of many trajectories
  \param[in] ampl_adi A [nadi x 1] or [nadi x ntraj] matrix of adiabatic amplitudes for
  one of many trajectories
  \param[in] lvl The level in the hierarchy of Hamiltonians at which we will perform the transformation
  \param[in] split The flag telling the function to perform the transformations of the matrix
  column by column, with each transformation being done by a children Hamiltonian (next level 
  relative to the input parameter lvl). This is only allowed if the number of columns in both
  ampl_dia and ampl_adi are equal to the number of children of the present-level Hamiltonian.

  There are 4 possible use cases:
  a) Both ampl_dia and ampl_adi are N x 1 matrices (one trajectory) and are meant to 
  be handled by the current Hamiltonian - this is already implemented above
  b) Both ampl_dia and ampl_adi are N x Ntraj matrices (Ntraj trajectories) and are still
  meant to be transformed by the same Hamiltonian - the current one. This case fits the 
  use case "a" due to the underlying matrix transformation in the called function.
  c, d) These use cases are analogous to "a" and "b", except for the fact that we need to
  call a lower-level Hamiltonian.
  Any of the above case can be combined with the "split" option to do the calculations 
  with the sub-Hamiltonians.

*/
  int i;

  if(lvl==level){  // Cases "a" and "b"
    if(split==0){  // Transform all the columns using the present level Hamiltonian
      ampl_adi2dia(ampl_dia, ampl_adi);
    }
    else if(split==1){
      // Check whether we have enough sub-Hamiltonians
      if(children.size()!=ampl_dia->n_cols){
        cout<<"ERROR in void nHamiltonian::ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
        cout<<"The number of columns of the ampl_dia ("<<ampl_dia->n_cols<<")";
        cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
        cout<<"Exiting...\n";
        exit(0);
      }
      if(children.size()!=ampl_adi->n_cols){
        cout<<"ERROR in void nHamiltonian::ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
        cout<<"The number of columns of the ampl_adi ("<<ampl_adi->n_cols<<")";
        cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
        cout<<"Exiting...\n";
        exit(0);
      }

      CMATRIX dia_tmp(ampl_dia->n_rows, 1);
      CMATRIX adi_tmp(ampl_adi->n_rows, 1);

      vector<int> stenc_dia1(ampl_dia->n_rows, 0);
      vector<int> stenc_adi1(ampl_dia->n_rows, 0);      
      vector<int> stenc_col(1, 0);

      for(i=0;i<ampl_dia->n_rows;i++){ stenc_dia1[i] = i;}
      for(i=0;i<ampl_adi->n_rows;i++){ stenc_adi1[i] = i;}

      for(i=0;i<children.size();i++){
        stenc_col[0] = i;

        pop_submatrix(ampl_adi, &adi_tmp, stenc_adi1, stenc_col);

        children[i]->ampl_adi2dia(dia_tmp, adi_tmp);

        push_submatrix(ampl_dia, &dia_tmp, stenc_dia1, stenc_col);
 
      }// for all children

    }// split==1
    else{
      cout<<"ERROR in void nHamiltonian::ampl_adi2dia(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split):\n";
      cout<<"The parameters split = "<<split<<" is not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }
  }// lvl == level

  else if(lvl>level){  // Cases "c" and "d"
  
    for(int i=0;i<children.size();i++){
      children[i]->ampl_adi2dia(ampl_dia, ampl_adi, lvl, split);
    }

  }// lvl >level

  else{
    cout<<"WARNING in nHamiltonian::ampl_dia2adi\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }
}



}// namespace libnhamiltonian
}// liblibra

