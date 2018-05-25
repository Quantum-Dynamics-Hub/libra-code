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
  \file nHamiltonian_compute_Ehrenfest_forces.cpp
  \brief The file implements the calculations of various variants of the Ehrenfest forces 
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




CMATRIX nHamiltonian::Ehrenfest_forces_dia_unit(CMATRIX& ampl_dia){
/**

  These are the Ehrenfest forces derived such the EOMs derived from the
  quntum-classical energy would conserve:

  H_qc = sum_i {p_i^2/2m_i} + <PSI|H|PSI>/<PSI|PSI>

  The wavefunction is expressed in the diabatic basis:

  |PSI> = |psi_dia> * C_dia

  and evolves according to the TD-SE:

  i * hbar * d|PSI>/dt = H |PSI> 

  for a systematic derivations, look here: 
  https://github.com/alexvakimov/Derivatory/blob/master/Ehrenfest.pdf

*/

  if(ovlp_dia_mem_status==0){ cout<<"Error in Ehrenfest_forces_dia_unit(): the overlap matrix in the diabatic basis is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(ham_dia_mem_status==0){ cout<<"Error in Ehrenfest_forces_dia_unit(): the diabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }


  CMATRIX res(nnucl,1);
  CMATRIX* dtilda; dtilda = new CMATRIX(nadi,nadi);
  CMATRIX* invS; invS = new CMATRIX(nadi, nadi); 

  FullPivLU_inverse(*ovlp_dia, *invS);

  complex<double> norm = ( ampl_dia.H() * (*ovlp_dia) * ampl_dia ).M[0]; 

  
  for(int n=0;n<nnucl;n++){

      if(d1ham_dia_mem_status[n]==0){ cout<<"Error in Ehrenfest_forces_dia_unit(): the derivatives of the Hamiltonian matrix in the \
      diabatic basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }

      if(dc1_dia_mem_status[n]==0){ cout<<"Error in Ehrenfest_forces_dia_unit(): the derivatives couplings matrix in the diabatic \
      basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }


      *dtilda = (*dc1_dia[n]).H() * (*invS) * (*ham_dia);
      *dtilda = (*dtilda + (*dtilda).H() ) ;

      res.M[n] = -( ampl_dia.H() * (*d1ham_dia[n] - *dtilda ) * ampl_dia ).M[0];

  }// for n

  res /= norm; 


  delete dtilda;
  delete invS;

  return res;
 
}



CMATRIX nHamiltonian::Ehrenfest_forces_dia(CMATRIX& ampl_dia, int lvl){
/**
  \brief Computes the Ehrenfest forces in the diabatic basis

  \param[in] ampl_dia A [ndia x ntraj] matrix of diabatic amplitudes for
  one of many (ntraj) trajectories
  \param[in] lvl [0 or 1] - 0 - use the present level for all trajectories, 1 - use the next 
  level for the trajectories (one sub-Hamiltonian per each trajectory) 

  There are 2 possible use cases:
  a) lvl = 0 ampl_dia is a ndia x ntraj matrix and is meant to 
  be handled by the current Hamiltonian (e.g. like the NBRA) 

  b) lvl = 1 ampl_adi is a ndia x ntraj matrix (ntraj trajectories) and 
  each trajectory is meant to be handled by a separate sub-Hamiltonian 

*/
  int i;

  if(lvl==1){
    // Check whether we have enough sub-Hamiltonians
    if(children.size()!=ampl_dia.n_cols){
      cout<<"ERROR in CMATRIX nHamiltonian::Ehrenfest_forces_dia(CMATRIX& ampl_dia):\n";
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

  for(i=0;i<ampl_dia.n_rows;i++){ stenc_ampl[i] = i;}

  for(i=0;i<ampl_dia.n_cols;i++){
    stenc_col[0] = i;

    pop_submatrix(ampl_dia, ampl_tmp, stenc_ampl, stenc_col);

    if(lvl==0){
        frc_tmp = Ehrenfest_forces_dia_unit(ampl_tmp);
    }
    if(lvl==1){
        frc_tmp = children[i]->Ehrenfest_forces_dia_unit(ampl_tmp);
    }

    push_submatrix(F, frc_tmp, stenc_frc, stenc_col);
 
  }// for all children

  return F;
  
}

/*
CMATRIX nHamiltonian::Ehrenfest_forces_dia(CMATRIX& ampl_dia, vector<int>& id_){
//
//  See the description of the Ehrenfest_forces_dia(CMATRIX& ampl_dia) function
//
  if(id_.size()==1){
    if(id_[0]==id){   return Ehrenfest_forces_dia(ampl_dia);    }
    else{ cout<<"ERROR in Ehrenfest_forces_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->Ehrenfest_forces_dia(ampl_dia, next);
  }
}
*/



CMATRIX nHamiltonian::Ehrenfest_forces_adi_unit(CMATRIX& ampl_adi){
/**

  These are the Ehrenfest forces derived such the EOMs derived from the
  quntum-classical energy would conserve:

  H_qc = sum_i {p_i^2/2m_i} + <PSI|H|PSI>/<PSI|PSI>

  The wavefunction is expressed in the adiabatic basis:

  |PSI> = |psi_adi> * C_adi

  and evolves according to the TD-SE:

  i * hbar * d|PSI>/dt = H |PSI> 


  Some useful theory can be found here: 
  http://www.theochem.ruhr-uni-bochum.de/~nikos.doltsinis/nic_10_doltsinis.pdf

  for a systematic derivations, look here: 
  https://github.com/alexvakimov/Derivatory/blob/master/Ehrenfest.pdf

*/

  if(ham_adi_mem_status==0){ cout<<"Error in Ehrenfest_forces_adi(): the adiabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }


  complex<double> norm = (ampl_adi.H() * ampl_adi).M[0]; 

  CMATRIX res(nnucl,1);

  CMATRIX* tmp; tmp = new CMATRIX(nadi, nadi);


  for(int n=0;n<nnucl;n++){

    if(d1ham_adi_mem_status[n]==0){ cout<<"Error in Ehrenfest_forces_adi_unit(): the derivatives of the Hamiltonian matrix in the \
    adiabatic basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }

    if(dc1_adi_mem_status[n]==0){ cout<<"Error in Ehrenfest_forces_adi_unit(): the derivatives couplings matrix in the adiabatic \
    basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }


    *tmp = (*dc1_adi[n]).H() * (*ham_adi);
    *tmp = (*tmp + (*tmp).H() );

    res.M[n] = -( ampl_adi.H() * (*d1ham_adi[n] - *tmp ) * ampl_adi ).M[0];

  }// for n

  res /= norm; 
  delete tmp;

  return res;
}



CMATRIX nHamiltonian::Ehrenfest_forces_adi(CMATRIX& ampl_adi, int lvl){
/**
  \brief Computes the Ehrenfest forces in the adiabatic basis

  \param[in] ampl_adi A [nadi x ntraj] matrix of adiabatic amplitudes for
  one of many (ntraj) trajectories
  \param[in] lvl [0 or 1] - 0 - use the present level for all trajectories, 1 - use the next 
  level for the trajectories (one sub-Hamiltonian per each trajectory) 

  There are 2 possible use cases:
  a) lvl = 0 ampl_adi is a nadi x ntraj matrix and is meant to 
  be handled by the current Hamiltonian (e.g. like the NBRA) 

  b) lvl = 1 ampl_adi is a nadi x ntraj matrix (ntraj trajectories) and 
  each trajectory is meant to be handled by a separate sub-Hamiltonian 
*/
  int i;

  if(lvl==1){
    // Check whether we have enough sub-Hamiltonians
    if(children.size()!=ampl_adi.n_cols){
      cout<<"ERROR in CMATRIX nHamiltonian::Ehrenfest_forces_adi(CMATRIX& ampl_adi):\n";
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

  for(i=0;i<ampl_adi.n_rows;i++){ stenc_ampl[i] = i;}

  for(i=0;i<ampl_adi.n_cols;i++){
    stenc_col[0] = i;

    pop_submatrix(ampl_adi, ampl_tmp, stenc_ampl, stenc_col);

    if(lvl==0){
        frc_tmp = Ehrenfest_forces_adi_unit(ampl_tmp);
    }
    if(lvl==1){
        frc_tmp = children[i]->Ehrenfest_forces_adi_unit(ampl_tmp);
    }

    push_submatrix(F, frc_tmp, stenc_frc, stenc_col);
 
  }// for all children

  return F;
}


/*
CMATRIX nHamiltonian::Ehrenfest_forces_adi(CMATRIX& ampl_adi, vector<int>& id_){
//
//  See the description of the Ehrenfest_forces_adi(CMATRIX& ampl_adi) function
//
  if(id_.size()==1){
    if(id_[0]==id){   return Ehrenfest_forces_adi(ampl_adi);    }
    else{ cout<<"ERROR in Ehrenfest_forces_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->Ehrenfest_forces_adi(ampl_adi, next);
  }
}
*/




}// namespace libhamiltonian_generic
}// namespace libhamiltonian
}// liblibra

