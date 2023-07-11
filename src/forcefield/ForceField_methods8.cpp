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




int ForceField::get_mb_parameters(map<string,double>& prms                           /*Outputs*/
                                   ){
/******************************************************************
  Function returns the status of parameters assigned:
  1 - if the parameters were successfully obtained
  0 - otherwise
******************************************************************/
  int status = 0;

  // Using cutoff:
  int is_cut = 0;
  prms["is_cutoff"] = 0;
  double R_on,R_off;

  if(mb_functional=="Ewald_3D"){
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
    elec_etha *= (1.0/Angst);

    if(is_cut){
      prms["is_cutoff"] = 1;
      prms["R_on"]  = R_on;
      prms["R_off"] = R_off;
      prms["R_on2"] = R_on*R_on;
      prms["R_off2"]= R_off*R_off;
      status = 1;
    }
    // Default values:
    prms["elec_etha"] = 3.0 * (1.0/Angst); 
    if(is_elec_etha){  prms["elec_etha"] = elec_etha; }

  }// Ewald_3D

  else if(mb_functional=="vdw_LJ"||mb_functional=="vdw_LJ1"||mb_functional=="LJ_Coulomb"){

    if((is_R_vdw_off==1)&&(is_R_vdw_on==1)){ is_cut = 1; R_off = R_vdw_off; R_on = R_vdw_on; }
    else if((is_R_vdw_off==0)&&(is_R_vdw_on==0)){}
    else{
      if(is_R_vdw_off==1){// Hence R_on not defined
        is_cut = 1;
        R_off  = R_vdw_off;
        R_on   = R_off - 1.0;
      }else{            // Hence R_off not defined
        is_cut = 1;
        R_on   = R_vdw_on;
        R_off  = R_on + 1.0;
      }
    }

    R_on *= Angst;
    R_off *= Angst;

    if(is_cut){
      prms["is_cutoff"] = 1;
      prms["R_on"]  = R_on;
      prms["R_off"] = R_off;
      prms["R_on2"] = R_on*R_on;
      prms["R_off2"]= R_off*R_off;
      status = 1;
    }

  }// vdw_LJ

  prms["time"] = 0;

  return status;
}


int ForceField::get_mb_excl_parameters(map<string,double>& prms                           /*Outputs*/
                                   ){
/******************************************************************
  Function returns the status of parameters assigned:
  1 - if the parameters were successfully obtained
  0 - otherwise
******************************************************************/
  int status = 0;

  // Using cutoff:
  int is_cut = 0;
  prms["is_cutoff"] = 0;
  double R_on,R_off;

  if(mb_excl_functional=="Ewald_3D"){
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
    elec_etha *= (1.0/Angst);

    if(is_cut){
      prms["is_cutoff"] = 1;
      prms["R_on"]  = R_on;
      prms["R_off"] = R_off;
      prms["R_on2"] = R_on*R_on;
      prms["R_off2"]= R_off*R_off;
      status = 1;
    }
    // Default values:
    prms["elec_etha"] = 3.0 * (1.0/Angst);
    if(is_elec_etha){  prms["elec_etha"] = elec_etha; }

  }// Ewald_3D

  else if(mb_excl_functional=="vdw_LJ"||mb_excl_functional=="vdw_LJ1"){

    if((is_R_vdw_off==1)&&(is_R_vdw_on==1)){ is_cut = 1; R_off = R_vdw_off; R_on = R_vdw_on; }
    else if((is_R_vdw_off==0)&&(is_R_vdw_on==0)){}
    else{
      if(is_R_vdw_off==1){// Hence R_on not defined
        is_cut = 1;
        R_off  = R_vdw_off;
        R_on   = R_off - 1.0;
      }else{            // Hence R_off not defined
        is_cut = 1;
        R_on   = R_vdw_on;
        R_off  = R_on + 1.0;
      }
    }

    // Convert to atomic units
    R_on *= Angst;
    R_off *= Angst;

    if(is_cut){
      prms["is_cutoff"] = 1;
      prms["R_on"]  = R_on;
      prms["R_off"] = R_off;
      prms["R_on2"] = R_on*R_on;
      prms["R_off2"]= R_off*R_off;
      status = 1;
    }

  }// vdw_LJ

  prms["time"] = 0;

  return status;
}


}// namespace libforcefield
}// namespace liblibra

