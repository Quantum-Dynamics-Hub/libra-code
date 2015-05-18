#include "Control_Parameters.h"
#include "Model_Parameters.h"
#include "units.h"

/*********************************************************************************
  This file contains the following functions:


  void set_parameters_hf(Control_Parameters& prms,Model_Parameters& modprms,vector<AO>& basis_ao)  

*********************************************************************************/

void set_parameters_hf(Control_Parameters& prms,Model_Parameters& modprms,vector<AO>& basis_ao){


  int a,b,c,d,A,B,C,D;
  int Norb = basis_ao.size(); // total number of AOs 


  // Formation of the Fock matrix: add Coulomb and Exchange parts
  for(a=0;a<Norb;a++){
    for(b=0;b<Norb;b++){
      for(c=0;c<Norb;c++){
        for(d=0;d<Norb;d++){

          //  (P_cd * (ab|cd) - P_alp_cd*(ad|cb))
          double J_abcd = ELECTRON_REPULSION_INTEGRAL(basis_ao[a],basis_ao[b],basis_ao[c],basis_ao[d]);
          double K_adcb = ELECTRON_REPULSION_INTEGRAL(basis_ao[a],basis_ao[d],basis_ao[c],basis_ao[b]);


          modprms.hf_int.set_JK_values(a,b,c,d,J_abcd,K_adcb);


        }// for d
      }// for c
    }// for b
  }// for a



}
