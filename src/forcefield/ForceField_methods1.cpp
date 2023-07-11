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
#include "../math_specialfunctions/libspecialfunctions.h"

/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace libio;

namespace libforcefield{


void ForceField::vdw_sigma_rule(string ff_type1,string ff_type2,double& sigma,int& is_sigma){

  int ff_type1_indx,ff_type2_indx;
  double R_i,R_j,A_i,A_j,a_i,a_j,gamma_ij;
  is_sigma = 1;

  // Find indexes of Atom_Record corresponding to atom types ff_type1 and ff_type2
  ff_type1_indx = Atom_Record_Index(ff_type1);
  ff_type2_indx = Atom_Record_Index(ff_type2);

  if(ff_type1_indx==-1 || ff_type2_indx==-1){ is_sigma = 0; }
  else{

    if(sigma_comb_rule=="SLATER-KIRKWOOD"){
      if(Atom_Records[ff_type1_indx].is_Atom_A_scale){ A_i = Atom_Records[ff_type1_indx].Atom_A_scale;} else{ is_sigma = 0;}
      if(Atom_Records[ff_type2_indx].is_Atom_A_scale){ A_j = Atom_Records[ff_type2_indx].Atom_A_scale;} else{ is_sigma = 0;}
      if(Atom_Records[ff_type1_indx].is_Atom_alpha){ a_i = Atom_Records[ff_type1_indx].Atom_alpha;} else{ is_sigma = 0;}
      if(Atom_Records[ff_type2_indx].is_Atom_alpha){ a_j = Atom_Records[ff_type2_indx].Atom_alpha;} else{ is_sigma = 0;}
    }
    else{
        if(Atom_Records[ff_type1_indx].is_Atom_sigma){ R_i = Atom_Records[ff_type1_indx].Atom_sigma; } else{ is_sigma = 0;}
        if(Atom_Records[ff_type2_indx].is_Atom_sigma){ R_j = Atom_Records[ff_type2_indx].Atom_sigma; } else{ is_sigma = 0;}
    }

    // calculations of auxiliary parameters
    if(is_sigma){
      if(is_sigma_comb_rule){
        if(sigma_comb_rule=="SUM"){ sigma = (R_i + R_j); }
        else if(sigma_comb_rule=="ARITHMETIC") {  sigma = 0.5*(R_i + R_j); }
        else if(sigma_comb_rule=="GEOMETRIC")  {  sigma = sqrt(R_i * R_j); }
        else if(sigma_comb_rule=="SLATER-KIRKWOOD"){
          // For this option we need other parameters
          R_i = A_i * POW(a_i,0.25);
          R_j = A_j * POW(a_j,0.25);
          gamma_ij =(R_i - R_j)/(R_i + R_j);
          sigma = 0.5*(R_i + R_j)*(1.0 + 0.2*(1.0 - exp(-12.0*gamma_ij * gamma_ij)));
        }
        else{
            std::cout<<"Warning: Unknown ForceField_sigma_comb_rule. Setting it to GEOMETRIC\n";
            sigma_comb_rule="GEOMETRIC";  is_sigma_comb_rule=1;
            sigma = sqrt(R_i * R_j);
        }
      }
    }// is_sigma
  }
}


void ForceField::vdw_epsilon_rule(string ff_type1,string ff_type2,double& epsilon,int& is_epsilon){

  int ff_type1_indx,ff_type2_indx;
  double R_i,R_j,D_i,D_j,R_ij,r_0,r_0_3,r_0_6;
  double A_i,A_j,a_i,a_j,N_i,N_j,G_i,G_j,gamma_ij;
  is_epsilon = 1;

  // Find indexes of Atom_Record corresponding to atom types ff_type1 and ff_type2
  ff_type1_indx = Atom_Record_Index(ff_type1);
  ff_type2_indx = Atom_Record_Index(ff_type2);

  if(epsilon_comb_rule=="SLATER-KIRKWOOD"){
    if(Atom_Records[ff_type1_indx].is_Atom_A_scale){ A_i = Atom_Records[ff_type1_indx].Atom_A_scale; } else{ is_epsilon = 0; }
    if(Atom_Records[ff_type2_indx].is_Atom_A_scale){ A_j = Atom_Records[ff_type2_indx].Atom_A_scale; } else{ is_epsilon = 0; }
    if(Atom_Records[ff_type1_indx].is_Atom_G_scale){ G_i = Atom_Records[ff_type1_indx].Atom_G_scale; } else{ is_epsilon = 0; }
    if(Atom_Records[ff_type2_indx].is_Atom_G_scale){ G_j = Atom_Records[ff_type2_indx].Atom_G_scale; } else{ is_epsilon = 0; }
    if(Atom_Records[ff_type1_indx].is_Atom_alpha){ a_i = Atom_Records[ff_type1_indx].Atom_alpha;} else{ is_epsilon = 0; }
    if(Atom_Records[ff_type2_indx].is_Atom_alpha){ a_j = Atom_Records[ff_type2_indx].Atom_alpha;} else{ is_epsilon = 0; }
    if(Atom_Records[ff_type1_indx].is_Atom_N_eff){ N_i = Atom_Records[ff_type1_indx].Atom_N_eff;} else{ is_epsilon = 0; }
    if(Atom_Records[ff_type2_indx].is_Atom_N_eff){ N_j = Atom_Records[ff_type2_indx].Atom_N_eff;} else{ is_epsilon = 0; }
  }
  else{
    if(Atom_Records[ff_type1_indx].is_Atom_epsilon){ D_i = Atom_Records[ff_type1_indx].Atom_epsilon;} else{ is_epsilon = 0; }
    if(Atom_Records[ff_type2_indx].is_Atom_epsilon){ D_j = Atom_Records[ff_type2_indx].Atom_epsilon;} else{ is_epsilon = 0; }
  }

  if(is_epsilon){
    if(is_epsilon_comb_rule){
      if(epsilon_comb_rule=="PROD"){   epsilon = (D_i * D_j);    }
      else if(epsilon_comb_rule=="GEOMETRIC")  {  epsilon = sqrt(D_i * D_j); }
      else if(epsilon_comb_rule=="SLATER-KIRKWOOD"){
        // For this option we need other parameters
        R_i = A_i * POW(a_i,0.25);
        R_j = A_j * POW(a_j,0.25);
        gamma_ij =(R_i - R_j)/(R_i + R_j);
        R_ij = 0.5*(R_i + R_j)*(1.0 + 0.2*(1.0 - exp(-12.0*gamma_ij * gamma_ij)));
        r_0 = (1.0/R_ij);
        r_0_3 = r_0 * r_0 * r_0;
        r_0_6 = r_0_3 * r_0_3;
        epsilon = 181.16 * G_i * G_j * a_i * a_j /(sqrt(a_i/N_i) + sqrt(a_j/N_j))*r_0_6;
      }
      else{
        std::cout<<"Warning: Unknown ForceField_epsilon_comb_rule. Setting it to GEOMETRIC\n";
        epsilon_comb_rule="GEOMETRIC";   is_epsilon_comb_rule=1;
        epsilon = sqrt(D_i * D_j);
      }
    }// if rule
  }// is_epsilon

}

int ForceField::get_vdw_parameters(string ff_type1, string ff_type2,string excl_pair, /*Inputs*/
                                   map<string,double>& prms                           /*Outputs*/
                                   ){
/******************************************************************
  Function returns the status of parameters assigned:
  1 - if the parameters were successfully obtained
  0 - otherwise
******************************************************************/
  double sigma,epsilon,D,r0,alpha,K;
  int is_sigma,is_epsilon,is_D,is_r0,is_alpha,is_K;
  is_sigma = is_epsilon = is_D = is_r0 = is_alpha = is_K = 0;
  int status = 0;

  //------------- Apply rules to calculate missing parameters --------------
  if(!is_sigma){ vdw_sigma_rule(ff_type1,ff_type2,sigma,is_sigma); }
  if(!is_epsilon){ vdw_epsilon_rule(ff_type1,ff_type2,epsilon,is_epsilon); }

  //------------ Additional relations --------------------------------------
  if(!is_alpha){ if(is_D && is_K){ alpha = sqrt(0.5*K/D); is_alpha = 1;}  }
  if(!is_D)    { if(is_K && is_alpha){ D = 0.5*K/(alpha*alpha); is_D = 1; }  }
  if(!is_K)    { if(is_D && is_alpha){ K = 2.0*D*alpha*alpha; is_K = 1; } }

  // Assign parameters according to the force field and the potential used
  // Do necessary scaling of the pre-computed variables
  prms["scale"] = 1.0;
  if(excl_pair=="12"){ if(is_vdw_scale12){  prms["scale"] = vdw_scale12; } }
  else if(excl_pair=="13"){ if(is_vdw_scale13){  prms["scale"] = vdw_scale13;  } }
  else if(excl_pair=="14"){ if(is_vdw_scale14){  prms["scale"] = vdw_scale14;  } }

  // Using cutoff:
  int is_cut = 0;
  prms["is_cutoff"] = 0; // Default!
  double R_on,R_off;
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
  sigma *= Angst;
  r0 *= Angst;
  epsilon *= (1.0/hartree);
  D *= (1.0/hartree);
  alpha *= (1.0/Angst);


  if(is_cut){  
    prms["is_cutoff"] = 1;
    prms["R_on"]  = R_on;
    prms["R_off"] = R_off;
    prms["R_on2"] = R_on*R_on;
    prms["R_off2"]= R_off*R_off;
  }

  if(vdw_functional=="LJ"||vdw_functional=="Buffered14_7"||mb_functional=="vdw_LJ"||mb_functional=="vdw_LJ1"||mb_functional=="LJ_Coulomb"){
    prms["sigma"] = sigma;
    prms["epsilon"]= epsilon;
    status = is_sigma * is_epsilon ;
    if(epsilon==0.0){ status = 0; }
  }
  else if(vdw_functional=="Morse"){
    prms["D"] = D;
    prms["r0"]= r0;
    prms["alpha"]=alpha;
    status = is_K * is_r0 * is_alpha;
    if(D==0.0){ status = 0; }
  }
  prms["time"] = 0; // time since last re-calculation

  return status;

}

int ForceField::get_vdw_parameters(int sz,vector<string> types,double** epsilon, double** sigma){
// sz - number of atoms
// types - type of each atom
// epsilon - epsilon parameter of each atom
// sigma - sigma parameter of each atom

cout<<"In get_vdw_parameters:\n";

  int res = 1;
  for(int i=0;i<sz;i++){
    //-------------- Start with looking the whole record ---------------------
    // Find index of Atom_Record corresponding to the force field type types[i]
    int ff_atom1_indx = Atom_Record_Index(types[i]);
    if(ff_atom1_indx==-1){
      *epsilon[i] = 0.0; *sigma[i] = 0.0;
      cout<<"Warning: In ForceField::get_vdw_parameters : Can not fing the atom type ff_atom type in force field "<<ForceField_Name<<"\n";
      cout<<"Setting A and B constants to to zero\n";
    }
    else{
      map<string,double> d;
      get_vdw_parameters(types[i],types[i],"none",d);
      *epsilon[i] = sqrt(d["epsilon"]);
      *sigma[i] = sqrt(d["sigma"]);
    }
    cout<<"i = "<<i<<"type = "<<types[i]<<" epsilon = "<<*epsilon[i]<<" sigma = "<<*sigma[i]<<endl;
  }// for i

  return res;

}

}// namespace libforcefield
}// namespace liblibra


