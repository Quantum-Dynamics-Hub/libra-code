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


namespace libforcefield{



int ForceField::get_elec_parameters(string ff_type1, string ff_type2,std::string excl_pair, /*Inputs*/
                                    double q1,double q2,int is_q1,int is_q2,
                                    map<string,double>& prms                                /*Outputs*/
                                   ){
/******************************************************************
  Function returns the status of parameters assigned:
  1 - if the parameters were successfully obtained
  0 - otherwise
******************************************************************/
  double qi,qj,delta,eps;
  int is_delta,is_eps;  is_delta = is_eps = 0;
  int is_qi, is_qj; is_qi = is_q1; is_qj = is_q2;
  int status = 0;

  //-------------- Start with looking the whole record ---------------------
  // Find index of Atom_Record corresponding to the angle formed by atom types ff_type1 and ff_type2
  int ff_atom1_indx = Atom_Record_Index(ff_type1);
  int ff_atom2_indx = Atom_Record_Index(ff_type2);

  //------------- If not found - apply rules  ------------------------------
  if(!is_qi){
    if(ff_atom1_indx>-1){ // if the record has been found
      if(Atom_Records[ff_atom1_indx].is_Atom_partial_charge){ 
      qi = Atom_Records[ff_atom1_indx].Atom_partial_charge; is_qi = 1;    }
    }
  }else{ qi = q1; }
  if(!is_qj){
    if(ff_atom2_indx>-1){ // if the record has been found
      if(Atom_Records[ff_atom2_indx].is_Atom_partial_charge){  
      qj = Atom_Records[ff_atom2_indx].Atom_partial_charge; is_qj = 1;    }
    }
  }else{ qj = q2; }

  delta = 0.0; is_delta = 1;
  eps = 1.0; is_eps = 1;

  // Assign parameters according to the force field and the potential used
  // Do necessary scaling of the pre-computed variables
  prms["scale"] = 1.0;
  if(excl_pair=="12"){ if(is_elec_scale12){  prms["scale"] = elec_scale12; } }
  else if(excl_pair=="13"){ if(is_elec_scale13){  prms["scale"] = elec_scale13;  } }
  else if(excl_pair=="14"){ if(is_elec_scale14){  prms["scale"] = elec_scale14;  } }

  // Using cutoff:
  int is_cut = 0;
  prms["is_cutoff"] = 0;
  double R_on,R_off;
  if((is_R_elec_off==1)&&(is_R_elec_on==1)){ is_cut = 1; R_off = R_elec_off; R_on = R_elec_on; }
  else if((is_R_elec_off==0)&&(is_R_elec_on==0)){}
  else{
    if(is_R_elec_off==1){// Hence R_on not defined
      is_cut = 1;
      R_off  = R_elec_off;
      R_on   = R_off - 1.0;
    }else{            // Hence R_off not defined
      is_cut = 1;
      R_on   = R_elec_on;
      R_off  = R_on + 1.0;
    }
  }


  // Convert to atomic units
  R_on *= Angst;
  R_off *= Angst;
  delta *= Angst;
  //eps *= hartree; // eps is the permittivity, so in denominator


  if(is_cut){
    prms["is_cutoff"] = 1;
    prms["R_on"]  = R_on;
    prms["R_off"] = R_off;
    prms["R_on2"] = R_on*R_on;
    prms["R_off2"]= R_off*R_off;
  }

  if(elec_functional=="Coulomb"){     
    prms["q1"] = qi;
    prms["q2"] = qj; 
    prms["eps"]= eps;
    prms["delta"]=delta;
    status = is_qi * is_qj * is_eps * is_delta ;
    if(qi==0.0 && qj==0.0){ status = 0; }
  }
  return status;
}


}// namespace libforcefield
}// namespace liblibra

