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
  \file nHamiltonian_compute_Ehrenfest.cpp
  \brief The file implements the calculations of the Ehrenfest-related properties:
  energies, forces and state-resolved force contributions
    
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



complex<double> nHamiltonian::Ehrenfest_energy_dia(CMATRIX& ampl_dia){
/**
  Compute the expectation value of the Hamiltonian in the diabatic basis:

  This is an Ehrenfest energy for the superposition:

  |PSI> = |psi_dia> * C_dia 

  Here C_adi are the: 

  E = <PSI|H|PSI>/<PSI|PSI> = C_dia.H() * <psi_dia|H|psi_dia> * C_dia / C_dia.H() * S * C_dia

*/

  if(ovlp_dia_mem_status==0){ cout<<"Error in Ehrenfest_energy_dia(): the overlap matrix in the diabatic basis is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(ham_dia_mem_status==0){ cout<<"Error in Ehrenfest_energy_dia(): the diabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

//  if(ampl_dia_mem_status==0){ cout<<"Error in Ehrenfest_energy_dia(): the amplitudes of the diabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }


  complex<double> norm = (ampl_dia.H() * (*ovlp_dia) * ampl_dia).M[0]; 

  return (ampl_dia.H() * (*ham_dia) * ampl_dia).M[0]/norm;

}



complex<double> nHamiltonian::Ehrenfest_energy_dia(CMATRIX& ampl_dia, vector<int>& id_){
/**
  See the description of the Ehrenfest_energy_dia(CMATRIX& ampl_dia) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return Ehrenfest_energy_dia(ampl_dia);    }
    else{ cout<<"ERROR in force_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->Ehrenfest_energy_dia(ampl_dia, next);
  }
}




complex<double> nHamiltonian::Ehrenfest_energy_adi(CMATRIX& ampl_adi){
/**
  Compute the expectation value of the Hamiltonian in the adiabatic basis:

  This is an Ehrenfest energy for the superposition:

  |PSI> = |psi_adi> * C_adi 

  Here C_adi are the: 

  E = <PSI|H|PSI> = C_adi.H() * <psi_adi|H|psi_adi> * C_adi

*/

  if(ham_adi_mem_status==0){ cout<<"Error in Ehrenfest_energy_adi(): the adiabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

//  if(ampl_adi_mem_status==0){ cout<<"Error in Ehrenfest_energy_adi(): the amplitudes of the adiabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }


  complex<double> norm = (ampl_adi.H() * ampl_adi).M[0]; 
  
  return (ampl_adi.H() * (*ham_adi) * ampl_adi).M[0] / norm;

}


complex<double> nHamiltonian::Ehrenfest_energy_adi(CMATRIX& ampl_adi, vector<int>& id_){
/**
  See the description of the Ehrenfest_energy_adi(CMATRIX& ampl_adi) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return Ehrenfest_energy_adi(ampl_adi);    }
    else{ cout<<"ERROR in force_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->Ehrenfest_energy_adi(ampl_adi, next);
  }
}






vector<CMATRIX> nHamiltonian::Ehrenfest_forces_tens_adi(CMATRIX& ampl_adi){
/**
  The elements of the returned list are the matrixes F_adi^MF that, when multiplied
  by the adiabatic amplitudes give the Ehrenfest forces in the adiabatic basis

  f_adi^MF.M[n] = Cadi.H() * (F_adi^MF[n]) * Cadi  

  The normalization factor is embedded in F_adi^MF[n]
  
*/

  vector<CMATRIX> res; res = vector<CMATRIX>(nnucl, CMATRIX(nadi,nadi));

  if(ham_adi_mem_status==0){ cout<<"Error in Ehrenfest_forces_tens_adi(): the adiabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

//  if(ampl_adi_mem_status==0){ cout<<"Error in Ehrenfest_forces_tens_adi(): the amplitudes of the adiabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }



  complex<double> norm = (ampl_adi.H() * ampl_adi).M[0]; 

  CMATRIX* tmp; tmp = new CMATRIX(nadi, nadi);


  for(int n=0;n<nnucl;n++){

    if(d1ham_adi_mem_status[n]==0){ cout<<"Error in Ehrenfest_forces_tens_adi(): the derivatives of the Hamiltonian matrix in the \
    adiabatic basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }

    if(dc1_adi_mem_status[n]==0){ cout<<"Error in Ehrenfest_forces_tens_adi(): the derivatives couplings matrix in the adiabatic \
    basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }

    *tmp = (*dc1_adi[n]).H() * (*ham_adi);
    *tmp = (*tmp + (*tmp).H() );

    res[n] = (-1.0/norm) * (*d1ham_adi[n] - *tmp );

  }// for n

  delete tmp;

  return res;

}

vector<CMATRIX> nHamiltonian::Ehrenfest_forces_tens_adi(CMATRIX& ampl_adi, vector<int>& id_){
/**
  See the description of the Ehrenfest_forces_tens_adi(CMATRIX& ampl_adi) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return Ehrenfest_forces_tens_adi(ampl_adi);    }
    else{ cout<<"ERROR in Ehrenfest_forces_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->Ehrenfest_forces_tens_adi(ampl_adi, next);
  }
}




vector<CMATRIX> nHamiltonian::Ehrenfest_forces_tens_dia(CMATRIX& ampl_dia){
/**
  The elements of the returned list are the matrixes F_dia^MF that, when multiplied
  by the diabatic amplitudes give the Ehrenfest forces in the diabatic basis

  f_dia^MF.M[n] = Cdia.H() * (F_dia^MF[n]) * Cdia  

  The normalization factor is embedded in F_dia^MF[n]

*/

  vector<CMATRIX> res; res = vector<CMATRIX>(nnucl, CMATRIX(ndia,ndia));

  if(ovlp_dia_mem_status==0){ cout<<"Error in Ehrenfest_forces_tens_dia(): the overlap matrix in the diabatic basis is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(ham_dia_mem_status==0){ cout<<"Error in Ehrenfest_forces_tens_dia(): the diabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

//  if(ampl_dia_mem_status==0){ cout<<"Error in Ehrenfest_forces_tens_dia(): the amplitudes of the diabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }


  CMATRIX* dtilda; dtilda = new CMATRIX(nadi,nadi);
  CMATRIX* invS; invS = new CMATRIX(nadi, nadi); 

  FullPivLU_inverse(*ovlp_dia, *invS);
  complex<double> norm = (ampl_dia.H() * (*ovlp_dia) * ampl_dia).M[0]; 

  
  for(int n=0;n<nnucl;n++){

      if(d1ham_dia_mem_status[n]==0){ cout<<"Error in Ehrenfest_forces_tens_dia(): the derivatives of the Hamiltonian matrix in the \
      diabatic basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }

      if(dc1_dia_mem_status[n]==0){ cout<<"Error in Ehrenfest_forces_tens_dia(): the derivatives couplings matrix in the diabatic \
      basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for the calculations \n"; exit(0); }


      *dtilda = (*dc1_dia[n]).H() * (*invS) * (*ham_dia);
      *dtilda = (*dtilda + (*dtilda).H() ) ;

      res[n] = (-1.0/norm) * (*d1ham_dia[n] - *dtilda );

  }// for n


  delete dtilda;
  delete invS;

  return res;
 
}

vector<CMATRIX> nHamiltonian::Ehrenfest_forces_tens_dia(CMATRIX& ampl_dia, vector<int>& id_){
/**
  See the description of the Ehrenfest_forces_tens_dia(CMATRIX& ampl_dia) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return Ehrenfest_forces_tens_dia(ampl_dia);    }
    else{ cout<<"ERROR in Ehrenfest_forces_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->Ehrenfest_forces_tens_dia(ampl_dia, next);
  }
}







}// namespace libnhamiltonian
}// liblibra

