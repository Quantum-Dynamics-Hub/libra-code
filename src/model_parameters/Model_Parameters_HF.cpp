/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Model_Parameters.h"


/// liblibra namespace
namespace liblibra{


namespace libmodel_parameters{


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
          double J_abcd = electron_repulsion_integral(basis_ao[a],basis_ao[b],basis_ao[c],basis_ao[d]);
          double K_adcb = electron_repulsion_integral(basis_ao[a],basis_ao[d],basis_ao[c],basis_ao[b]);


          modprms.hf_int.set_JK_values(a,b,c,d,J_abcd,K_adcb);


        }// for d
      }// for c
    }// for b
  }// for a



}


}// namespace libmodel_parameters
}// namespace liblibra

