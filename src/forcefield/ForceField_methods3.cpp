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



void ForceField::angle_theta_0_rule(std::string ff_type1,std::string ff_type2,std::string ff_type3,
                                    double bond_order12,double bond_order23,int coordination,
                                    double& theta_0,int& is_theta_0){

  int ff_type2_indx;
  is_theta_0 = 1;
  // Find indexe of Atom_Record corresponding to atom types ff_type2
  ff_type2_indx = Atom_Record_Index(ff_type2);

  if(ff_type2_indx==-1){ is_theta_0 = 0; }
  else{
    // If we are at this point - than all atom types
    // are correct - now we collect information about
    // properties of atom types and use them to calculate
    // all necessary force field parameters
    if(ForceField_Name=="UFF"||ForceField_Name=="DREIDING"){
      if(Atom_Records[ff_type2_indx].is_Atom_theta){ theta_0 = deg_to_rad*Atom_Records[ff_type2_indx].Atom_theta; }else{ is_theta_0 = 0; }
    }
    else if(ForceField_Name=="GAFF"){ is_theta_0 = 0; }
    else if(ForceField_Name=="MMFF94"){is_theta_0 = 0; }
    else{ is_theta_0 = 0; }

  }// if ( ff_type2_indx>-1)

}

void ForceField::angle_k_theta_rule(std::string ff_type1,std::string ff_type2,std::string ff_type3,
                                    double bond_order12,double bond_order23,int coordination,
                                    double& k_theta,int& is_k_theta){

  int ff_type1_indx,ff_type2_indx,ff_type3_indx;
  is_k_theta = 1; k_theta = 0.0;
  if(ForceField_Name=="DREIDING"){ k_theta = 100.0; }
  else{
    // Find indexe of Atom_Record corresponding to atom types ff_type2
    ff_type1_indx = Atom_Record_Index(ff_type1);
    ff_type2_indx = Atom_Record_Index(ff_type2);
    ff_type3_indx = Atom_Record_Index(ff_type3);

    if(ff_type1_indx==-1 || ff_type2_indx==-1 || ff_type3_indx==-1){ is_k_theta = 0; }
    else{
      if(ForceField_Name=="UFF"){
        double theta_0,rij,rjk,Zi,Zk;
        int is_theta_0,is_rij,is_rjk,is_Zi,is_Zk;
        is_theta_0 = is_rij = is_rjk = is_Zi = is_Zk = 0;

        if(Atom_Records[ff_type2_indx].is_Atom_theta){ theta_0 = deg_to_rad*Atom_Records[ff_type2_indx].Atom_theta; is_theta_0 = 1;}
        if(Atom_Records[ff_type1_indx].is_Atom_Z_star){ Zi = Atom_Records[ff_type1_indx].Atom_Z_star; is_Zi = 1;}
        if(Atom_Records[ff_type3_indx].is_Atom_Z_star){ Zk = Atom_Records[ff_type3_indx].Atom_Z_star; is_Zk = 1;}
        bond_r0_rule(ff_type1,ff_type2,bond_order12,rij,is_rij);
        bond_r0_rule(ff_type2,ff_type3,bond_order23,rjk,is_rjk);

        if(is_rij && is_rjk && is_theta_0 && is_Zi && is_Zk){
          double rik2,rik,cs;
          rik2 = rij*rij + rjk*rjk - 2.0*rij*rjk*cos(theta_0); rik = sqrt(rik2); 
          cs = cos(theta_0);
          k_theta = 664.12*(Zi*Zk/(rik2*rik2*rik))*(3.0*rij*rjk*(1.0-cs*cs)-rik2*cs);
        }
        else{is_k_theta = 0;}
      }
      else{ is_k_theta = 0; }
    }// found all indexes
  }// all force fields except DREIDING

}

int ForceField::get_angle_parameters(string ff_type1, string ff_type2,string ff_type3,
                                     double bond_order12,double bond_order23,int coordination, /*Inputs*/
                                     map<string,double>& prms                                 /*Outputs*/
                                    ){
/******************************************************************
  Function returns the status of parameters assigned:
  1 - if the parameters were successfully obtained
  0 - otherwise
******************************************************************/

  cout<<"In ForceField::get_angle_parameters()...\n";
  cout<<"ff_type1 = "<<ff_type1<<"ff_type2 = "<<ff_type2<<"ff_type3 = "<<ff_type3<<endl;

  double k_theta,theta_0,cos_theta_0,C0,C1,C2;
  int is_k_theta,is_theta_0,is_cos_theta_0,is_C0,is_C1,is_C2;
  is_k_theta = is_theta_0 = is_cos_theta_0 = is_C0 = is_C1 = is_C2 = 0;
  int status = 0;
  //-------------- Start with looking the whole record ---------------------
  // Find index of Angle_Record corresponding to the angle formed by atom types ff_type1, ff_type2 and ff_type3
  int ff_angle_indx = Angle_Record_Index(ff_type1,ff_type2,ff_type3);

  //------------- If not found - apply rules  ------------------------------
  if(ff_angle_indx>-1){ // if the record has been found
    cout<<"Angle record has been found in the Force Field used\n";
    // Equilibrium angle
    if(Angle_Records[ff_angle_indx].is_Angle_theta_eq){
      theta_0 = deg_to_rad*Angle_Records[ff_angle_indx].Angle_theta_eq; is_theta_0 = 1;
    }
    // Angle bending force constant
    if(Angle_Records[ff_angle_indx].is_Angle_k_angle){
      k_theta = Angle_Records[ff_angle_indx].Angle_k_angle; is_k_theta = 1;
    }
//    // Angle coordination
//    if(Angle_Records[ff_angle_indx].is_Angle_k_angle){
//      k_theta = Angle_Records[ff_angle_indx].Angle_k_angle; is_k_theta = 1;
//    }


  }

  //------------- Apply rules to calculate missing parameters --------------
  if(!is_theta_0){  
    cout<<"Angle theta_0 is not defined, so using the rule\n";
    angle_theta_0_rule(ff_type1,ff_type2,ff_type3,bond_order12,bond_order23,coordination,theta_0,is_theta_0); 
  } 
  if(!is_k_theta){
    cout<<"Angle k_theta is not defined, so using the rule\n";
    angle_k_theta_rule(ff_type1,ff_type2,ff_type3,bond_order12,bond_order23,coordination,k_theta,is_k_theta);
  }

  //------------ Additional relations --------------------------------------
  if(!is_theta_0){  if(is_cos_theta_0) { theta_0 = acos(cos_theta_0); is_theta_0 = 1; } }
  if(!is_cos_theta_0){ if(is_theta_0) { cos_theta_0 = cos(theta_0); is_cos_theta_0 = 1; } }
  if(!is_C0){  if(is_theta_0){   double cs,sn; cs = cos(theta_0); sn = sin(theta_0);  C0 = (2.0*cs*cs + 1.0)/(4.0*sn*sn); is_C0 = 1;  }   }
  if(!is_C1){  if(is_theta_0){   double cs,sn; cs = cos(theta_0); sn = sin(theta_0);  C1 = -cs/(sn*sn); is_C1 = 1;    }  }
  if(!is_C2){  if(is_theta_0){   double sn; sn = sin(theta_0);                        C2 = 1.0/(4.0*sn*sn); is_C2 = 1;  }  }

  cout<<"Here is what we have so far (in conventional units)...\n";
  cout<<" is_theta_0 = "<<is_theta_0<<" is_k_theta = "<<is_k_theta<<" is_cos_theta_0 = "<<is_cos_theta_0
      <<" is_C0 = "<<is_C0<<" is_C1 = "<<is_C1<<" is_C2 = "<<is_C2<<endl;
  cout<<" theta_0 = "<<theta_0<<" k_theta = "<<k_theta<<" is_cos_theta_0 = "<<cos_theta_0
      <<" C0 = "<<C0<<" C1 = "<<C1<<" C2 = "<<C2<<endl;


  // Convert to atomic units
  k_theta *= (1.0/hartree);


  // Assign parameters according to the force field and the potential used 
  // Do necessary scaling of the pre-computed variables
  cout<<"The parameters determined will be used with the angle bending functional = "<<angle_functional<<endl;
  
  if(angle_functional=="Harmonic"){
    prms["k_theta"] = 0.5*k_theta;
    prms["theta_0"]= theta_0;
    prms["coordination"] = coordination;
    status = is_k_theta * is_theta_0 ;
    if(k_theta==0.0){ status = 0; }
  }
  else if(angle_functional=="Fourier"){
    prms["coordination"] = coordination;
    prms["k_theta"]= k_theta;
    prms["C0"]=C0;
    prms["C1"]=C1;
    prms["C2"]=C2; 
    prms["theta_0"] = theta_0;
    prms["cos_theta_0"] = cos_theta_0;
    status = is_k_theta * is_C0 * is_C1 * is_C2;
    if(k_theta==0.0){ status = 0; }
  }
  else if(angle_functional=="Harmonic_Cos"){
    prms["coordination"] = coordination;
    prms["theta_0"] = theta_0;
    prms["cos_theta_0"] = cos_theta_0;
    prms["k_theta"]= k_theta;
    status = is_k_theta * is_cos_theta_0;
    if(k_theta==0.0){ status = 0; }
  }
  else if(angle_functional=="Harmonic_Cos_General"){
    prms["coordination"] = coordination;
    prms["theta_0"] = theta_0;
    prms["cos_theta_0"] = cos_theta_0;
    prms["k_theta"]= 0.5*k_theta/(sin(theta_0)*sin(theta_0));
    status = is_k_theta * is_cos_theta_0;
    if(k_theta==0.0){ status = 0; }
  }
  cout<<"The outcome is: the parameters are ";
  if(status){ cout<<" sufficient\n"; }
  else{ cout<<" insufficient (the interaction will not be added)\n"; }

  return status;

}


}// namespace libforcefield
}// namespace liblibra

