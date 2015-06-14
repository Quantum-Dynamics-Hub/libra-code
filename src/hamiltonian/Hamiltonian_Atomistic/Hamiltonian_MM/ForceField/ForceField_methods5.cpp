#include "../ForceField.h"

int ForceField::get_oop_parameters(string ff_type1, string ff_type2, /*Inputs*/
                                   string ff_type3, string ff_type4,
                                   map<string,double>& prms          /*Outputs*/
                                  ){
/******************************************************************
  Function returns the status of parameters assigned:
  1 - if the parameters were successfully obtained
  0 - otherwise 
******************************************************************/
  double Vphi,Vphi1,Vphi2,Vphi3,phi0;
  int n,opt;
  int is_Vphi,is_Vphi1,is_Vphi2,is_Vphi3,is_phi0,is_n,is_opt;
  is_Vphi = is_Vphi1 = is_Vphi2 = is_Vphi3 = is_phi0 = is_n = is_opt = 0;
  int status = 0;
  //-------------- Start with looking the whole record ---------------------
  // Find index of Bond_Record corresponding to the bond formed by atom types
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
/*
  //------------- Apply rules to calculate missing parameters --------------
  if(!is_r0){  bond_r0_rule(ff_type1,ff_type2,bond_order,r0,is_r0); } 
  if(!is_K) {  bond_K_rule(ff_type1,ff_type2,bond_order,K,is_K); }
  if(!is_D) {  bond_D_rule(ff_type1,ff_type2,bond_order,D,is_D); }
  if(!is_alpha){ bond_alpha_rule(ff_type1,ff_type2,bond_order,alpha,is_alpha); }

  //------------ Additional relations --------------------------------------
  if(!is_alpha){ if(is_D && is_K){ alpha = sqrt(0.5*K/D); is_alpha = 1;}  }
  if(!is_D)    { if(is_K && is_alpha){ D = 0.5*K/(alpha*alpha); is_D = 1; }  }
  if(!is_K)    { if(is_D && is_alpha){ K = 2.0*D*alpha*alpha; is_K = 1; } }


  // Assign parameters according to the force field and the potential used 
  // Do necessary scaling of the pre-computed variables
  if(bond_functional=="Harmonic"||bond_functional=="Quartic"){
    prms["K"] = K;
    prms["r0"]= r0;
    status = is_K * is_r0 ;
  }
  else if(bond_functional=="Morse"){
    prms["D"] = D;
    prms["r0"]= r0;
    prms["alpha"]=alpha;
    status = is_K * is_r0 * is_alpha;
  }
*/
  return status;
}

