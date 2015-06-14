#include "../ForceField.h"


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
    if(is_cut){
      prms["is_cutoff"] = 1;
      prms["R_on"]  = R_on;
      prms["R_off"] = R_off;
      prms["R_on2"] = R_on*R_on;
      prms["R_off2"]= R_off*R_off;
      status = 1;
    }
    // Default values:
    prms["elec_etha"] = 3.0; 
    if(is_elec_etha){  prms["elec_etha"] = elec_etha; }

  }// Ewald_3D

  else if(mb_functional=="vdw_LJ"||mb_functional=="vdw_LJ1"){

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
    if(is_cut){
      prms["is_cutoff"] = 1;
      prms["R_on"]  = R_on;
      prms["R_off"] = R_off;
      prms["R_on2"] = R_on*R_on;
      prms["R_off2"]= R_off*R_off;
      status = 1;
    }
    // Default values:
    prms["elec_etha"] = 3.0;
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


