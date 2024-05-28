/*********************************************************************************
* Copyright (C) 2017-2023 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian_compute_QTSH_forces.cpp
  \brief The file implements the calculations of various variants of the QTSH forces 
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




CMATRIX nHamiltonian::QTSH_forces_dia_unit(CMATRIX& ampl_dia, int option){
/**
  \param[in] ampl_dia: MATRIX(ndia, 1) diabatic amplitudes for one trajectory

  Returns:
  MATRIX(ndof, 1) - QTSH nonclassical forces in diabatic representation, for a single trajectory
 
  The nonclassical force in QTSH has the form of the off-diagonal terms in the Ehrenfest force.

*/

  if(ovlp_dia_mem_status==0){ cout<<"Error in QTSH_forces_dia_unit(): the overlap matrix in the diabatic basis is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(ham_dia_mem_status==0){ cout<<"Error in QTSH_forces_dia_unit(): the diabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }


  CMATRIX res(nnucl,1);
  CMATRIX* dtilda; dtilda = new CMATRIX(nadi,nadi);
  CMATRIX* invS; invS = new CMATRIX(nadi, nadi); 

  FullPivLU_inverse(*ovlp_dia, *invS);

  complex<double> norm = ( ampl_dia.H() * (*ovlp_dia) * ampl_dia ).M[0]; 

  CMATRIX* temp_diag = new CMATRIX(nadi, nadi);
  
  for(int n=0;n<nnucl;n++){

      if(d1ham_dia_mem_status[n]==0){ cout<<"Error in QTSH_forces_dia_unit(): the derivatives of the Hamiltonian matrix in the \
      diabatic basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }

      if(dc1_dia_mem_status[n]==0){ cout<<"Error in QTSH_forces_dia_unit(): the derivatives couplings matrix in the diabatic \
      basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }


      if(option==0){
        *dtilda = ( (*dc1_dia[n]).H() * (*invS) * (*ham_dia) + (*ham_dia) * (*invS) * (*dc1_dia[n])  );

        for(int i=0;i<nadi;i++){
          temp_diag->set(i,i, (*d1ham_dia[n] - *dtilda ).get(i,i)); 
        }
        res.M[n] = -( ampl_dia.H() * (*d1ham_dia[n] - *dtilda - *temp_diag ) * ampl_dia ).M[0];
      }
      else if(option==1){
        for(int i=0;i<nadi;i++){
          temp_diag->set(i,i, d1ham_dia[n]->get(i,i)); 
        }
        res.M[n] = -( ampl_dia.H() * (*d1ham_dia[n] - *temp_diag ) * ampl_dia ).M[0];
      }
      

  }// for n

  res /= norm; 


  delete dtilda;
  delete invS;
  delete temp_diag;

  return res;
 
}



CMATRIX nHamiltonian::QTSH_forces_dia_unit(CMATRIX& ampl_dia){
  return QTSH_forces_dia_unit(ampl_dia, 0);
}



CMATRIX nHamiltonian::QTSH_forces_dia(CMATRIX& ampl_dia, int lvl, int option){
/**
  \brief Computes the QTSH forces in the diabatic basis

  \param[in] ampl_dia [ndia x ntraj] matrix of diabatic amplitudes for
  one or many (ntraj) trajectories
  \param[in] lvl [0 or 1] - 0 - use the present level for all trajectories, 1 - use the next 
  level for the trajectories (one sub-Hamiltonian per each trajectory) 

  There are 2 possible use cases:
  a) lvl = 0 ampl_dia is a ndia x ntraj matrix and is meant to 
  be handled by the current Hamiltonian (e.g. like the NBRA) 

  b) lvl = 1 ampl_adi is a ndia x ntraj matrix (ntraj trajectories) and 
  each trajectory is meant to be handled by a separate sub-Hamiltonian 

  Returns:
  MATRIX(ndof, ntraj) - QTSH forces in diabatic representation, for multiple trajectories


*/
  int i;

  if(lvl==1){
    // Check whether we have enough sub-Hamiltonians
    if(children.size()!=ampl_dia.n_cols){
      cout<<"ERROR in CMATRIX nHamiltonian::QTSH_forces_dia(CMATRIX& ampl_dia):\n";
      cout<<"The number of columns of the ampl_dia ("<<ampl_dia.n_cols<<")";
      cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
      cout<<"Exiting...\n";
      exit(0);
    }
  }

  CMATRIX ampl_tmp(ampl_dia.n_rows, 1);
  CMATRIX frc_tmp(nnucl, 1);
  CMATRIX F(nnucl, ampl_dia.n_cols);

  vector<int> stenc_ampl(ampl_dia.n_rows, 0);
  vector<int> stenc_frc(nnucl, 0);
  vector<int> stenc_col(1, 0);

  // The indicies in stenc_frc reflects which indicies are to be updated in the final F matrix
  for(i=0;i<nnucl;i++) { stenc_frc[i] = i;}    

  for(i=0;i<ampl_dia.n_rows;i++){ stenc_ampl[i] = i;}

  for(i=0;i<ampl_dia.n_cols;i++){
    stenc_col[0] = i;

    pop_submatrix(ampl_dia, ampl_tmp, stenc_ampl, stenc_col);

    if(lvl==0){
        frc_tmp = QTSH_forces_dia_unit(ampl_tmp, option);
    }
    if(lvl==1){
        frc_tmp = children[i]->QTSH_forces_dia_unit(ampl_tmp, option);
    }

    push_submatrix(F, frc_tmp, stenc_frc, stenc_col);
 
  }// for all children

  return F;
  
}


CMATRIX nHamiltonian::QTSH_forces_dia(CMATRIX& ampl_dia, int lvl){
  return QTSH_forces_dia(ampl_dia, lvl, 0);
}




/*
CMATRIX nHamiltonian::QTSH_forces_dia(CMATRIX& ampl_dia, vector<int>& id_){
//
//  See the description of the QTSH_forces_dia(CMATRIX& ampl_dia) function
//
  if(id_.size()==1){
    if(id_[0]==id){   return QTSH_forces_dia(ampl_dia);    }
    else{ cout<<"ERROR in QTSH_forces_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->QTSH_forces_dia(ampl_dia, next);
  }
}
*/



CMATRIX nHamiltonian::QTSH_forces_adi_unit(CMATRIX& ampl_adi, int option, CMATRIX& transform){
/**
  \param[in] ampl_adi: MATRIX(nadi, 1) diabatic amplitudes for one trajectory

  \params[in] option [0 or 1] - option 0 keeps all the terms in the Ehrenfest-like force expression, including NAC
  option 1 removes all the derivative NACs - this is to enforce the local diabatization approximation, to be
  consistent with it
  
  The nonclassical force in QTSH has the form of the off-diagonal terms in the Ehrenfest force.

  Returns:
  MATRIX(ndof, 1) - QTSH forces in adiabatic representation, for single trajectory


*/

  if(ham_adi_mem_status==0){ cout<<"Error in QTSH_forces_adi(): the adiabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }


  complex<double> norm = (ampl_adi.H() * ampl_adi).M[0];


  CMATRIX res(nnucl,1);
  CMATRIX tmp(nadi, nadi);

//  CMATRIX& T = transform;
  CMATRIX T(transform);  T.identity();

  CMATRIX temp_diag(nadi, nadi);

  for(int n=0;n<nnucl;n++){

    if(d1ham_adi_mem_status[n]==0){ cout<<"Error in QTSH_forces_adi_unit(): the derivatives of the Hamiltonian matrix in the \
    adiabatic basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }

    if(dc1_adi_mem_status[n]==0){ cout<<"Error in QTSH_forces_adi_unit(): the derivatives couplings matrix in the adiabatic \
    basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }


    if(option==0){ // Original formulation with NACs - for non-LD integrators

      tmp = (T.H() * (*dc1_adi[n]) * T ).H() *  (T.H() * (*ham_adi) * T);      
      tmp = tmp + tmp.H();
      
      for(int i=0;i<nadi;i++){
        temp_diag.set(i,i, ( T.H() * (*d1ham_adi[n]) *T - tmp ).get(i,i)); 
      }
      res.M[n] = -( ampl_adi.H() * ( T.H() * (*d1ham_adi[n]) * T - tmp - temp_diag ) * ampl_adi ).M[0];
    }
    else if(option==1){ // Options that disregard the NACs - appropriate for the LD integrators
      for(int i=0;i<nadi;i++){
        temp_diag.set(i,i, d1ham_adi[n]->get(i,i) ); 
      }
      res.M[n] = -( ampl_adi.H() * T.H() * ( *d1ham_adi[n] - temp_diag ) * T * ampl_adi ).M[0];
    }

  }// for n

  res /= norm; 
  //delete tmp;

  return res;
}

CMATRIX nHamiltonian::QTSH_forces_adi_unit(CMATRIX& ampl_adi, int option){
  CMATRIX I(nadi, nadi); I.identity();
  return QTSH_forces_adi_unit(ampl_adi,option, I);
}

CMATRIX nHamiltonian::QTSH_forces_adi_unit(CMATRIX& ampl_adi){
  CMATRIX I(nadi, nadi); I.identity();
  return QTSH_forces_adi_unit(ampl_adi, 0, I);
}


CMATRIX nHamiltonian::QTSH_forces_adi(CMATRIX& ampl_adi, int lvl, int option, vector<CMATRIX*>& transforms){
/**
  \brief Computes the QTSH forces in the adiabatic basis

  \param[in] ampl_adi [nadi x ntraj] matrix of adiabatic amplitudes for
  one of many (ntraj) trajectories
  \param[in] lvl [0 or 1] - 0 - use the present level for all trajectories, 1 - use the next 
  level for the trajectories (one sub-Hamiltonian per each trajectory) 

  There are 2 possible use cases:
  a) lvl = 0 ampl_adi is a nadi x ntraj matrix and is meant to 
  be handled by the current Hamiltonian (e.g. like the NBRA) 

  b) lvl = 1 ampl_adi is a nadi x ntraj matrix (ntraj trajectories) and 
  each trajectory is meant to be handled by a separate sub-Hamiltonian 

  \params[in] option [0 or 1] - option 0 keeps all the terms in the QTSH force expression, including NAC
  option 1 removes all the derivative NACs - this is to enforce the local diabatization approximation, to be
  consistent with it

  Returns:
  MATRIX(ndof, ntraj) - QTSH forces in adiabatic representation, for multiple trajectories

*/
  int i;

  if(lvl==1){
    // Check whether we have enough sub-Hamiltonians
    if(children.size()!=ampl_adi.n_cols){
      cout<<"ERROR in CMATRIX nHamiltonian::QTSH_forces_adi(CMATRIX& ampl_adi):\n";
      cout<<"The number of columns of the ampl_adi ("<<ampl_adi.n_cols<<")";
      cout<<" should be equal to the number of children Hamiltonians ("<<children.size()<<")\n";
      cout<<"Exiting...\n";
      exit(0);
    }
  }


  CMATRIX ampl_tmp(ampl_adi.n_rows, 1);
  CMATRIX frc_tmp(nnucl, 1);
  CMATRIX F(nnucl, ampl_adi.n_cols);

  vector<int> stenc_ampl(ampl_adi.n_rows, 0);
  vector<int> stenc_frc(nnucl, 0);
  vector<int> stenc_col(1, 0);

  // The indicies in stenc_frc reflects which indicies are to be updated in the final F matrix
  for(i=0;i<nnucl;i++) { stenc_frc[i] = i;}   

  for(i=0;i<ampl_adi.n_rows;i++){ stenc_ampl[i] = i;}

  for(i=0;i<ampl_adi.n_cols;i++){
    stenc_col[0] = i;

    pop_submatrix(ampl_adi, ampl_tmp, stenc_ampl, stenc_col);

    if(lvl==0){
        frc_tmp = QTSH_forces_adi_unit(ampl_tmp, option);
    }
    if(lvl==1){
        frc_tmp = children[i]->QTSH_forces_adi_unit(ampl_tmp, option, *transforms[i] );
    }

    push_submatrix(F, frc_tmp, stenc_frc, stenc_col);
 
  }// for all children

  return F;
}

CMATRIX nHamiltonian::QTSH_forces_adi(CMATRIX& ampl_adi, int lvl, int option){
  int ntraj = ampl_adi.n_cols;
  vector<CMATRIX*> I(ntraj);

  for(int itraj=0; itraj<ntraj; itraj++){
    I[itraj] = new CMATRIX(nadi, nadi);
    I[itraj]->load_identity();
  }      

  CMATRIX F(nnucl, ampl_adi.n_cols);

  F = QTSH_forces_adi(ampl_adi, lvl, option, I);

  for(int itraj=0; itraj<ntraj; itraj++){  delete I[itraj];  }  I.clear();

  return F;

}

CMATRIX nHamiltonian::QTSH_forces_adi(CMATRIX& ampl_adi, int lvl){
  int ntraj = ampl_adi.n_cols;
  vector<CMATRIX*> I(ntraj);

  for(int itraj=0; itraj<ntraj; itraj++){
    I[itraj] = new CMATRIX(nadi, nadi);
    I[itraj]->load_identity();
  }  

  CMATRIX F(nnucl, ampl_adi.n_cols);

  F = QTSH_forces_adi(ampl_adi, lvl, 0, I);

  for(int itraj=0; itraj<ntraj; itraj++){  delete I[itraj];  }  I.clear();
  return F;
}


/*
CMATRIX nHamiltonian::QTSH_forces_adi(CMATRIX& ampl_adi, vector<int>& id_){
//
//  See the description of the QTSH_forces_adi(CMATRIX& ampl_adi) function
//
  if(id_.size()==1){
    if(id_[0]==id){   return QTSH_forces_adi(ampl_adi);    }
    else{ cout<<"ERROR in QTSH_forces_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->QTSH_forces_adi(ampl_adi, next);
  }
}
*/




}// namespace libnhamiltonian
}// liblibra

