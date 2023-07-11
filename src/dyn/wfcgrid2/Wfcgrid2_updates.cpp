/*********************************************************************************
* Copyright (C) 2019-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Wfcgrid2_updates.cpp
  \brief The file implements functions various types of computations and updates
    
*/

#include "Wfcgrid2.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

using namespace libwfcgrid;

/// libwfcgrid namespace
namespace libwfcgrid2{


void Wfcgrid2::update_Hamiltonian(bp::object py_funct, bp::object params, int rep){
/**
  \brief Update the Hamiltonian for nd-D grid
  \param[in] py_funct the Python function to be called to compute the Hamiltonian matrix
  in a given representation
  \param[in] params the parameters fed into the Python function that performs the calculations
  \param[in] rep representation: 0 - diabatic, 1 - adiabatic

  This function recomputes the Hamiltonian for all points
  working in atomic units: hbar = 1
*/

  int idof, ipt;
  int has_attr;

  double res; res = 0.0;
  double nrm; nrm = 0.0;

  for(int npt1=0; npt1<Npts; npt1++){

    MATRIX q(ndof, 1);
    for(idof=0; idof<ndof; idof++){
      ipt = gmap[npt1][idof];
      q.set(idof, rgrid[idof]->get(ipt) );
    }

    // Call the Python function with such arguments
    bp::object obj = py_funct(bp::object(q), params);  


    // Try extract the adiabatic Hamiltonian
    has_attr = (int)hasattr(obj,"v_complex");
    if(has_attr){   Vcomplex[npt1] = extract<CMATRIX>(obj.attr("v_complex"));    } 


    if(rep==0){

      // Try to extract the diabatic Hamiltonian
      has_attr = (int)hasattr(obj,"ham_dia");
      if(has_attr){  Hdia[npt1] = extract<CMATRIX>(obj.attr("ham_dia"));    } 

    }// rep == 0

    else if(rep==1){

      // Try extract the adiabatic Hamiltonian
      has_attr = (int)hasattr(obj,"ham_adi");
      if(has_attr){   Hadi[npt1] = extract<CMATRIX>(obj.attr("ham_adi"));    } 

      // Try extract the derivative couplings in adiabatic representation
      has_attr = (int)hasattr(obj,"dc1_adi");        
      if(has_attr){
        vector<CMATRIX> _dc1_adi(ndof, CMATRIX(nstates,nstates));
        _dc1_adi = extract<CMATRIXList>(obj.attr("dc1_adi"));    

        for(idof=0; idof<ndof; idof++){  NAC1[npt1][idof] = _dc1_adi[idof];  }
      }

    }// rep == 1


  }// for allgrid points

}// update_Hamiltonian





void Wfcgrid2::update_adiabatic(){
/**
  Update the adiabatic wfc from the diabatic: DIA -> ADI

  Should have called ```update_Hamiltonian``` and then ```update_propagator_H``` prior to this

  Hdia * U = S * U * Hadi

  |psi_adi> = |psi_dia> * U  <=>  phi_i^adi = sum_{a} [ U_ai * phi_a^dia ]

  Then: <adi|H|adi> = U.H() * <dia|H|dia> * U or:

  Hdia * U = U * Hadi

*/

  for(int npt1=0; npt1<Npts; npt1++){ PSI_adi[npt1] = U[npt1].T() * PSI_dia[npt1];  }
  
}

void Wfcgrid2::update_diabatic(){
/**
  Update the diabatic wfc from the adiabatic: ADI -> DIA

  Should have called ```update_Hamiltonian``` and then ```update_propagator_H``` prior to this

  Hdia * U = S * U * Hadi

  |psi_adi> = |psi_dia> * U  <=>  phi_i^adi = sum_{a} [ U_ai * phi_a^dia ]

  Then: <adi|H|adi> = U.H() * <dia|H|dia> * U or:

  Hdia * U = U * Hadi

*/

  for(int npt1=0; npt1<Npts; npt1++){ PSI_dia[npt1] = U[npt1].conj() * PSI_adi[npt1];  }
  
}





}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

