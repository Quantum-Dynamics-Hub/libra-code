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

#include "ForceField.h"

/// liblibra namespace
namespace liblibra{

//using namespace liblinalg;
//using namespace libio;


namespace libforcefield{


void ForceField::oop_rule(std::string ff_type1, std::string ff_type2, std::string ff_type3, std::string ff_type4,
                          double& K, int& is_K, 
                          double& C0, int& is_C0,  double& C1, int& is_C1,
                          double& C2, int& is_C2){
// ff_type2 - the central atom, connected to ff_type1, ff_type3, and ff_type4

  K = 0.0;
  C0 = C1 = C2 = 0.0;
  is_K = is_C0 = is_C1 = is_C2 = 0;
  
  if(ForceField_Name=="UFF"){
    std::string elt1, elt2, elt3, elt4;
    elt1 = ff_type1.c_str()[0] + ff_type1.c_str()[1];
    elt2 = ff_type2.c_str()[0] + ff_type2.c_str()[1];
    elt3 = ff_type3.c_str()[0] + ff_type3.c_str()[1];
    elt4 = ff_type4.c_str()[0] + ff_type4.c_str()[1];
 
    if(ff_type2=="C_2" || ff_type2=="C_R" ){  C0 = 1.0; C1 = -1.0; C2 = 0.0;  K = 6.0;
    is_K = is_C0 = is_C1 = is_C2 = 1;
    }
    else if(elt2=="N_" ){
      K = 0.0; 
      
      // If at least two other atoms are C_2 or C_R - treat this Nitrogen as aromatic
      int a1 = (elt1=="C_2" || elt1=="C_R");
      int a3 = (elt3=="C_2" || elt3=="C_R");
      int a4 = (elt4=="C_2" || elt4=="C_R");

      if(a1+a3+a4>=2){
        C0 = 1.0; C1 = -1.0; C2 = 0.0;  K = 6.0;
        is_K = is_C0 = is_C1 = is_C2 = 1;
      }
    }

    // NEED TO DEFINE P, As, Sb, Bi !!!

  }// UFF

}



int ForceField::get_oop_parameters(string ff_type1, string ff_type2, /*Inputs*/
                                   string ff_type3, string ff_type4,
                                   map<string,double>& prms          /*Outputs*/
                                  ){
/******************************************************************
  Function returns the status of parameters assigned:
  1 - if the parameters were successfully obtained
  0 - otherwise 
******************************************************************/

  double Vphi,Vphi1,Vphi2,Vphi3,phi0,K,C0,C1,C2;
  int n,opt;
  int is_Vphi,is_Vphi1,is_Vphi2,is_Vphi3,is_phi0,is_n,is_opt, is_K, is_C0, is_C1, is_C2;
  is_Vphi = is_Vphi1 = is_Vphi2 = is_Vphi3 = is_phi0 = is_n = is_opt = is_K = is_C0 = is_C1 = is_C2 = 0;
  int status = 0;
  
  //-------------- Start with looking the whole record ---------------------
  // Find index of Dihedral_Record corresponding to the bond formed by atom types
  // ff_type1, ff_type2, ff_type3 and ff_type4
  int ff_dihedral_indx = Dihedral_Record_Index(ff_type1,ff_type2,ff_type3,ff_type4);

  //------------- If not found - apply rules  ------------------------------
  if(ff_dihedral_indx>-1){ // if the record has been found
    // Force constant
    if(Dihedral_Records[ff_dihedral_indx].is_Dihedral_vphi){
      Vphi = Dihedral_Records[ff_dihedral_indx].Dihedral_vphi; is_Vphi = 1;
    }
    // Force constant
    if(Dihedral_Records[ff_dihedral_indx].is_Dihedral_vphi1){
      Vphi1 = Dihedral_Records[ff_dihedral_indx].Dihedral_vphi1; is_Vphi1 = 1;
    }
    // Force constant
    if(Dihedral_Records[ff_dihedral_indx].is_Dihedral_vphi2){
      Vphi2 = Dihedral_Records[ff_dihedral_indx].Dihedral_vphi2; is_Vphi2 = 1;
    }
    // Force constant
    if(Dihedral_Records[ff_dihedral_indx].is_Dihedral_vphi3){
      Vphi3 = Dihedral_Records[ff_dihedral_indx].Dihedral_vphi3; is_Vphi3 = 1;
    }
    // Phase shift - equilibrium dihedral/torsion angle
    if(Dihedral_Records[ff_dihedral_indx].is_Dihedral_phase){
      phi0 = Dihedral_Records[ff_dihedral_indx].Dihedral_phase; is_phi0 = 1;
    }
    // Potntial periodicity/multiplicity
    if(Dihedral_Records[ff_dihedral_indx].is_Dihedral_mult){
      n = Dihedral_Records[ff_dihedral_indx].Dihedral_mult; is_n = 1;
    }

  }

  //------------- Apply rules to calculate missing parameters --------------
  if(1){  
    oop_rule(ff_type1, ff_type2, ff_type3, ff_type4, K, is_K, C0, is_C0, C1, is_C1, C2, is_C2); 
  }


  // Convert to atomic units
  K *= (1.0/hartree);

  // Assign parameters according to the force field and the potential used 
  // Do necessary scaling of the pre-computed variables
  if(oop_functional=="Fourier"){
    prms["K"] = K;
    prms["C0"] = C0;
    prms["C1"] = C1;
    prms["C2"] = C2;
    prms["opt"] = -1;   // TORSION
    status = is_K * is_C0 * is_C1 * is_C2;
  }

  return status;

}


}// namespace libforcefield
}// namespace liblibra

