#include "ForceField.h"

namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_mm{
namespace libforcefield{



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



  // Convert to atomic units
  //Vphi *= (1.0/hartree);
  //Vphi1 *= (1.0/hartree);
  //Vphi2 *= (1.0/hartree);
  //Vphi3 *= (1.0/hartree);

/*
  // Assign parameters according to the force field and the potential used 
  // Do necessary scaling of the pre-computed variables
  if(oop_functional=="General0"||oop_functional=="General1"||
     oop_functional=="General2"||oop_functional=="General3"){
    prms["Vphi"] = 0.5*Vphi;
    prms["phi0"]= phi0;
    prms["n"]=n;
    status = is_Vphi * is_phi0 * is_n ;
  }
  else if(oop_functional=="Fourier0"||oop_functional=="Fourier1"){
    prms["Vphi1"] = Vphi1;
    prms["Vphi2"] = Vphi2;
    prms["Vphi3"] = Vphi3;
    status = is_Vphi1 * is_Vphi2 * is_Vphi3;
  }
*/

  return status;

}


}// namespace libforcefield
}// namespace libhamiltonian_mm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian

