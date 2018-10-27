/*********************************************************************************
* Copyright (C) 2017-2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian_compute_forces.cpp
  \brief The file implements computations of various types of forces and force tensors
    
*/


#include <stdlib.h>

#include "nHamiltonian.h"
#include "../../math_meigen/libmeigen.h"


/// liblibra namespace
namespace liblibra{

/// libhamiltonian namespace 
namespace libhamiltonian{

/// libhamiltonian_generic namespace 
namespace libhamiltonian_generic{

using namespace liblinalg;
using namespace libmeigen;





vector<CMATRIX> nHamiltonian::forces_tens_adi(CMATRIX& ampl_adi){
/**

  The dependence of U on R is neglected

  The elements of the returned list are the matrixes F_adi that, when multiplied
  by the adiabatic amplitudes give the derivative of the total energy in the adiabatic basis

  f_adi.M[n] = Cadi.H() * (F_adi[n]) * Cadi   = -dE/dR

  The normalization factor is embedded in F_adi[n]

*/

  vector<CMATRIX> res; res = vector<CMATRIX>(nnucl, CMATRIX(nadi,nadi));

//  if(ampl_adi_mem_status==0){ cout<<"Error in forces_tens_adi(): the amplitudes of the adiabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }

  complex<double> norm = (ampl_adi.H() * ampl_adi).M[0]; 


  for(int n=0;n<nnucl;n++){

    if(d1ham_adi_mem_status[n]==0){ cout<<"Error in forces_tens_adi(): the derivatives of the Hamiltonian matrix in the \
    adiabatic basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }

    res[n] = (-1.0/norm)*(*d1ham_adi[n]); 

  }// for n

  return res;

}

vector<CMATRIX> nHamiltonian::forces_tens_adi(CMATRIX& ampl_adi, vector<int>& id_){
/**
  See the description of the forces_tens_adi(CMATRIX& ampl_adi) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return forces_tens_adi(ampl_adi);    }
    else{ cout<<"ERROR in forces_tens_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->forces_tens_adi(ampl_adi, next);
  }
}




CMATRIX nHamiltonian::forces_adi(CMATRIX& ampl_adi){
/**
  The dependence of U on R is neglected

  The elements of the returned matrix f_adi are the components 
  of the vector of the adiabatic force allong all DOFs

  Return : f_adi.M[n] = Cadi.H() * (F_adi[n]) * Cadi   = -dE/dR

*/
  CMATRIX res(nnucl, 1);
  vector<CMATRIX> dEdR;   dEdR = forces_tens_adi(ampl_adi);

//  if(ampl_adi_mem_status==0){ cout<<"Error in forces_adi(): the amplitudes of the adiabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }

  for(int n=0;n<nnucl;n++){  res.M[n] = (ampl_adi.H() * dEdR[n] * ampl_adi).M[0];   }

  return res;
}

CMATRIX nHamiltonian::forces_adi(CMATRIX& ampl_adi, vector<int>& id_){
/**
  See the description of the forces_adi(CMATRIX& ampl_adi) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return forces_adi(ampl_adi);    }
    else{ cout<<"ERROR in forces_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->forces_adi(ampl_adi, next);
  }
}


CMATRIX nHamiltonian::forces_adi(vector<int>& act_states){
/**
  Get a matrix of adiabatic Forces:  ndof x ntraj

  Each column contains the force for given trajectory (in the adiabatic state given by the "act_states")
  Each raw corresponds to the nuclear DOF

  This is a convenience function: meant to extract the children forces

*/

  int ntraj = act_states.size();
  int ndof = nnucl;

  vector<int> t1(ndof, 0); for(int dof=0;dof<ndof;dof++){  t1[dof] = dof; }
  vector<int> t2(1,0);
  vector<int> traj_id(2,0);

  CMATRIX F(ndof, ntraj);
  CMATRIX f(ndof, 1);
  CMATRIX Cadi(nadi,1); 

  for(int traj=0; traj<ntraj; traj++){
    t2[0] = traj;  
    traj_id[0] = 0;  traj_id[1] = traj;

    Cadi.set(act_states[traj],0, 1.0,0.0);
    f = forces_adi(Cadi, traj_id);

    push_submatrix(F, f, t1, t2);
  }

  return F;

}


CMATRIX nHamiltonian::forces_adi(int act_state){
/**
  Get a matrix of adiabatic Forces:  ndof x ntraj

  Each column contains the force for given trajectory (in the adiabatic state given by the "act_state")
  Each raw corresponds to the nuclear DOF

  This is a convenience function: meant to extract the children forces

*/

  int ntraj = children.size();
  int ndof = nnucl;

  vector<int> t1(ndof, 0); for(int dof=0;dof<ndof;dof++){  t1[dof] = dof; }
  vector<int> t2(1,0);
  vector<int> traj_id(2,0);

  CMATRIX F(ndof, ntraj);
  CMATRIX f(ndof, 1);
  CMATRIX Cadi(nadi,1); 

  for(int traj=0; traj<ntraj; traj++){
    t2[0] = traj;  
    traj_id[0] = 0;  traj_id[1] = traj;

    Cadi.set(act_state, 0, 1.0,0.0);
    f = forces_adi(Cadi, traj_id);

    push_submatrix(F, f, t1, t2);
  }

  return F;

}




vector<CMATRIX> nHamiltonian::forces_tens_dia(CMATRIX& ampl_dia){
/**

  The elements of the returned list are the matrixes F_dia that, when multiplied
  by the diabatic amplitudes give the derivative of the total energy in the diabatic basis

  f_dia.M[n] = Cdia.H() * (F_dia[n]) * Cdia   = -dE/dR

  The normalization factor is embedded in F_dia[n]

*/


  vector<CMATRIX> res; res = vector<CMATRIX>(nnucl, CMATRIX(ndia,ndia));

  if(ovlp_dia_mem_status==0){ cout<<"Error in forces_tens_dia(): the overlap matrix in the diabatic basis is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(basis_transform_mem_status==0){ cout<<"Error in forces_tens_dia(): the transformation basis matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

//  if(ampl_dia_mem_status==0){ cout<<"Error in forces_dia(): the amplitudes of the diabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }


  complex<double> norm = ( ampl_dia.H() * (*ovlp_dia) * ampl_dia).M[0]; 

  complex<double> Etot = (ampl_dia.H() * (*ham_dia) * ampl_dia).M[0]/norm;



  for(int n=0;n<nnucl;n++){

    if(d1ham_dia_mem_status[n]==0){ cout<<"Error in forces_tens_dia(): the derivatives of the Hamiltonian matrix in the \
    adiabatic basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }

    if(dc1_dia_mem_status[n]==0){ cout<<"Error in forces_tens_dia(): the derivatives couplings matrix in the adiabatic \
    basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }


    res[n] = -1.0 * (*d1ham_dia[n] - Etot * (*dc1_dia[n] + (*dc1_dia[n]).H()) ); 

  }// for n


  return res;
}

vector<CMATRIX> nHamiltonian::forces_tens_dia(CMATRIX& ampl_dia, vector<int>& id_){
/**
  See the description of the forces_tens_dia(CMATRIX& ampl_dia) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return forces_tens_dia(ampl_dia);    }
    else{ cout<<"ERROR in force_tens_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->forces_tens_dia(ampl_dia, next);
  }
}




CMATRIX nHamiltonian::forces_dia(CMATRIX& ampl_dia){
/**

  The dependence of U on R is neglected

  The elements of the returned matrix f_dia are the components 
  of the vector of the diabatic force allong all DOFs

  Return : f_dia.M[n] = Cdia.H() * (F_dia[n]) * Cdia   = -dE/dR

*/


  CMATRIX res(nnucl, 1);
  vector<CMATRIX> dEdR;   dEdR = forces_tens_dia(ampl_dia);

//  if(ampl_dia_mem_status==0){ cout<<"Error in forces_dia(): the amplitudes of the diabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }

  for(int n=0;n<nnucl;n++){  res.M[n] = ( ampl_dia.H() * dEdR[n] *  ampl_dia).M[0];   }

  return res;
}


CMATRIX nHamiltonian::forces_dia(CMATRIX& ampl_dia, vector<int>& id_){
/**
  See the description of the forces_dia(CMATRIX& ampl_dia) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return forces_dia(ampl_dia);    }
    else{ cout<<"ERROR in force_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->forces_dia(ampl_dia, next);
  }
}






}// namespace libhamiltonian_generic
}// namespace libhamiltonian
}// liblibra

