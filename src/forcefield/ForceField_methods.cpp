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



int ForceField::is_valid_atom_type(std::string ff_type){
  int res = 1;
  if(is_ForceField_Name){
    int indx = Atom_Record_Index(ff_type);
    if(indx==-1){ cout<<"Error: Force field "<<ForceField_Name<<" does not have atom type "<<ff_type<<endl; res = 0; }
  }
  else{  cout<<"Warning!: Force field name is not defined. But this allows you to manually define any atom type!\n";  res = 1;  }
  return 1;
}

int ForceField::is_valid_fragment_type(std::string ff_type){
  return 1;
}


void ForceField::bond_r0_rule(string ff_type1,string ff_type2,double bond_order,double& r0,int& is_r0){

  int ff_type1_indx, ff_type2_indx;
  double r1,r2;
  is_r0 = 1; 
  // Find indexes of Atom_Record corresponding to atom types ff_type1 and ff_type2
  ff_type1_indx = Atom_Record_Index(ff_type1);
  ff_type2_indx = Atom_Record_Index(ff_type2);

  if(ff_type1_indx==-1 || ff_type2_indx==-1){ is_r0 = 0; }
  else{
    // If we are at this point - than all atom types
    // are correct - now we collect information about
    // properties of atom types and use them to calculate
    // all necessary force field parameters
    if(Atom_Records[ff_type1_indx].is_Atom_radius){ r1 = Atom_Records[ff_type1_indx].Atom_radius; }else{ is_r0 = 0; }
    if(Atom_Records[ff_type2_indx].is_Atom_radius){ r2 = Atom_Records[ff_type2_indx].Atom_radius; }else{ is_r0 = 0; }

    if(is_r0){
      // Equilibrium bond length
      r0 = r1 + r2;

      if(ForceField_Name=="UFF"){
        // Bond order correction
        r0 += (-0.1332*(r0)*log(bond_order));

        // Electronegativity correction
        if( Atom_Records[ff_type1_indx].is_Atom_GMP &&
            Atom_Records[ff_type2_indx].is_Atom_GMP){
            double& xi1 = Atom_Records[ff_type1_indx].Atom_GMP;
            double& xi2 = Atom_Records[ff_type2_indx].Atom_GMP;
            r0 -= r1*r2*(xi1-2.0*sqrt(xi1*xi2)+xi2)/(xi1*r1+xi2*r2); 
        }

      }// UFF
      else if(ForceField_Name=="DREIDING"){   r0 -= 0.01;    }
      else if(ForceField_Name=="GAFF"){ is_r0 = 0;     }
      else if(ForceField_Name=="MMFF94"){is_r0 = 0; } 
      else{ is_r0 = 0; }

    }// if is_r0
  }// if (ff_type1_indx>=-1 && ff_type2_indx>-1)

 
}

void ForceField::bond_K_rule(string ff_type1,string ff_type2,double bond_order,double& K,int& is_K){

  int ff_type1_indx, ff_type2_indx;
  double r1,r2,z1,z2,r_tot;
  is_K = 1;
  // Find indexes of Atom_Record corresponding to atom types ff_type1 and ff_type2
  ff_type1_indx = Atom_Record_Index(ff_type1);
  ff_type2_indx = Atom_Record_Index(ff_type2);

  if(ff_type1_indx==-1 || ff_type2_indx==-1){ is_K = 0; }
  else{
    if(ForceField_Name=="UFF"){
      if(Atom_Records[ff_type1_indx].is_Atom_radius){ r1 = Atom_Records[ff_type1_indx].Atom_radius; }else{ is_K = 0; }
      if(Atom_Records[ff_type2_indx].is_Atom_radius){ r2 = Atom_Records[ff_type2_indx].Atom_radius; }else{ is_K = 0; }
      if(Atom_Records[ff_type1_indx].is_Atom_Z_star){ z1 = Atom_Records[ff_type1_indx].Atom_Z_star; }else{ is_K = 0; }
      if(Atom_Records[ff_type2_indx].is_Atom_Z_star){ z2 = Atom_Records[ff_type2_indx].Atom_Z_star; }else{ is_K = 0; }

      if(is_K){
        // Equilibrium bond length
        r_tot = r1 + r2;
        // Bond order correction
        double bo_corr = (-0.1332*(r_tot)*log(bond_order));
        r_tot += bo_corr;
        // Electronegativity correction
        if(Atom_Records[ff_type1_indx].is_Atom_GMP &&
          Atom_Records[ff_type2_indx].is_Atom_GMP){
          double xi1 = Atom_Records[ff_type1_indx].Atom_GMP;
          double xi2 = Atom_Records[ff_type2_indx].Atom_GMP; 
          double el_corr = r1*r2*(xi1-2.0*sqrt(xi1*xi2)+xi2)/(xi1*r1+xi2*r2);
          r_tot -= el_corr;
        }
        K = 664.12*(z1*z2)/(r_tot*r_tot*r_tot);
      }
    }// UFF
    else if(ForceField_Name=="DREIDING"){ K = bond_order*700.0;    }
    else if(ForceField_Name=="GAFF"){ is_K = 0;     }
    else if(ForceField_Name=="MMFF94"){is_K = 0; }
    else{ is_K = 0; }

  }// if (ff_type1_indx>=-1 && ff_type2_indx>-1)

}

void ForceField::bond_D_rule(string ff_type1,string ff_type2,double bond_order,double& D,int& is_D){

  is_D = 1;
  if(ForceField_Name=="UFF"){is_D = 0; }
  else if(ForceField_Name=="DREIDING"){ D = bond_order*70.0;    }
  else if(ForceField_Name=="GAFF"){ is_D = 0;     }
  else if(ForceField_Name=="MMFF94"){is_D = 0; }
  else{ is_D = 0; }

}

void ForceField::bond_alpha_rule(string ff_type1,string ff_type2,double bond_order,double& alpha,int& is_alpha){

  is_alpha = 1;
  if(ForceField_Name=="UFF"){is_alpha = 0; }
  else if(ForceField_Name=="DREIDING"){ alpha = sqrt(5.0);   }
  else if(ForceField_Name=="GAFF"){ is_alpha = 0;     }
  else if(ForceField_Name=="MMFF94"){is_alpha = 0; }
  else{ is_alpha = 0; }

}



int ForceField::get_bond_parameters(string ff_type1, string ff_type2, /*Inputs*/
                                    double bond_order,
                                    map<string,double>& prms          /*Outputs*/
                                    ){
/******************************************************************
  Function returns the status of parameters assigned:
  1 - if the parameters were successfully obtained
  0 - otherwise
******************************************************************/
  double K,D,r0,alpha;
  int is_K,is_D,is_r0,is_alpha;
  is_K = is_D = is_r0 = is_alpha = 0;
  int status = 0;
  //-------------- Start with looking the whole record ---------------------
  // Find index of Bond_Record corresponding to the bond formed by atom types ff_type1 and ff_type2
  int ff_bond_indx = Bond_Record_Index(ff_type1,ff_type2);

  //------------- If not found - apply rules  ------------------------------
  if(ff_bond_indx>-1){ // if the record has been found
    // Equilibrium bond length
    if(Bond_Records[ff_bond_indx].is_Bond_r_eq){
      r0 = Bond_Records[ff_bond_indx].Bond_r_eq; is_r0 = 1;
    }       
    // Force constants
    if(Bond_Records[ff_bond_indx].is_Bond_k_bond){
      K = Bond_Records[ff_bond_indx].Bond_k_bond; is_K = 1;
    }
    // Energy depth
    if(Bond_Records[ff_bond_indx].is_Bond_D_bond){
      D = Bond_Records[ff_bond_indx].Bond_D_bond; is_D = 1;
    }
    // Morse shape parameter
    if(Bond_Records[ff_bond_indx].is_Bond_alpha){
      alpha = Bond_Records[ff_bond_indx].Bond_alpha; is_alpha = 1;
    }
  }

  //------------- Apply rules to calculate missing parameters --------------
  if(!is_r0){  bond_r0_rule(ff_type1,ff_type2,bond_order,r0,is_r0); } 
  if(!is_K) {  bond_K_rule(ff_type1,ff_type2,bond_order,K,is_K); }
  if(!is_D) {  bond_D_rule(ff_type1,ff_type2,bond_order,D,is_D); }
  if(!is_alpha){ bond_alpha_rule(ff_type1,ff_type2,bond_order,alpha,is_alpha); }

  //------------ Additional relations --------------------------------------
  if(!is_alpha){ if(is_D && is_K){ alpha = sqrt(0.5*K/D); is_alpha = 1;}  }
  if(!is_D)    { if(is_K && is_alpha){ D = 0.5*K/(alpha*alpha); is_D = 1; }  }
  if(!is_K)    { if(is_D && is_alpha){ K = 2.0*D*alpha*alpha; is_K = 1; } }


  // Convert to atomic units
  r0 *= Angst;
  K *= (1.0/hartree)/(Angst*Angst);
  D *= (1.0/hartree);
  alpha *= (1.0/Angst);



  // Assign parameters according to the force field and the potential used 
  // Do necessary scaling of the pre-computed variables
  if(bond_functional=="Harmonic"||bond_functional=="Quartic"){
    prms["K"] = 0.5*K;
    prms["r0"]= r0;
    status = is_K * is_r0 ;
  }
  else if(bond_functional=="Morse"){
    prms["D"] = D;
    prms["r0"]= r0;
    prms["alpha"]=alpha;
    status = is_K * is_r0 * is_alpha;
  }
  return status;



}



}// namespace libforcefield
}// namespace liblibra


