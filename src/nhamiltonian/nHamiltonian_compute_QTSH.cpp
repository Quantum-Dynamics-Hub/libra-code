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
  \file nHamiltonian_compute_QTSH.cpp
  \brief The file implements the calculations of the QTSH-related properties:
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



complex<double> nHamiltonian::QTSH_energy_dia(CMATRIX& ampl_dia){
/**
  Compute the coherence energy of QTSH in the diabatic representation.
  The coherence energy in QTSH has the form of the off-diagonal terms of the Ehrenfest energy.
*/

  if(ovlp_dia_mem_status==0){ cout<<"Error in QTSH_energy_dia(): the overlap matrix in the diabatic basis is not allocated \
  but it is needed for the calculations\n"; exit(0); }

  if(ham_dia_mem_status==0){ cout<<"Error in QTSH_energy_dia(): the diabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

//  if(ampl_dia_mem_status==0){ cout<<"Error in QTSH_energy_dia(): the amplitudes of the diabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }

  complex<double> norm = (ampl_dia.H() * (*ovlp_dia) * ampl_dia).M[0]; 

  return (ampl_dia.H() * complex<double>(0.0, -1.0)* (*nac_dia) * ampl_dia).M[0] / norm;

}



complex<double> nHamiltonian::QTSH_energy_dia(CMATRIX& ampl_dia, vector<int>& id_){
/**
  See the description of the QTSH_energy_dia(CMATRIX& ampl_dia) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return QTSH_energy_dia(ampl_dia);    }
    else{ cout<<"ERROR in force_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->QTSH_energy_dia(ampl_dia, next);
  }
}




complex<double> nHamiltonian::QTSH_energy_adi(CMATRIX& ampl_adi, CMATRIX& transform){
/**
  Compute the coherence energy of QTSH in the adiabatic representation.
  The coherence energy in QTSH has the form of the off-diagonal terms of the Ehrenfest energy.
*/

  if(ham_adi_mem_status==0){ cout<<"Error in QTSH_energy_adi(): the adiabatic Hamiltonian matrix is not allocated \
  but it is needed for the calculations\n"; exit(0); }

//  if(ampl_adi_mem_status==0){ cout<<"Error in QTSH_energy_adi(): the amplitudes of the adiabatic states are\
//  not allocated, but they are needed for the calculations\n"; exit(0); }

//  CMATRIX& T = transform;
  CMATRIX T(transform); T.identity();
  
  complex<double> norm = (ampl_adi.H() * ampl_adi).M[0]; 
  
  return (ampl_adi.H() * ( T.H() * complex<double>(0.0, -1.0)* (*nac_adi) * T) * ampl_adi).M[0] / norm;

}

complex<double> nHamiltonian::QTSH_energy_adi(CMATRIX& ampl_adi){
  CMATRIX I(nadi, nadi); I.identity();
  return QTSH_energy_adi(ampl_adi, I);
}


complex<double> nHamiltonian::QTSH_energy_adi(CMATRIX& ampl_adi, vector<int>& id_, CMATRIX& transform){
/**
  See the description of the QTSH_energy_adi(CMATRIX& ampl_adi) function
*/
  if(id_.size()==1){
    if(id_[0]==id){   return QTSH_energy_adi(ampl_adi);    }
    else{ cout<<"ERROR in force_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->QTSH_energy_adi(ampl_adi, next, transform);
  }
}

complex<double> nHamiltonian::QTSH_energy_adi(CMATRIX& ampl_adi, vector<int>& id_){
  CMATRIX I(nadi, nadi); I.identity();
  return QTSH_energy_adi(ampl_adi, id_, I);
}

}// namespace libnhamiltonian
}// liblibra

