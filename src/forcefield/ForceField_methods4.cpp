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



void ForceField::dihedral_rule(std::string ff_type1,std::string ff_type2,
                               std::string ff_type3,std::string ff_type4,
                               double bond_order12,double bond_order23,double bond_order34,
                               double& Vphi,int& is_Vphi,
                               double& phi0,int& is_phi0,
                               int& n,int& is_n){
  Vphi = 0.0;
  is_n = is_phi0 = is_Vphi = 1;
  
  if(ForceField_Name=="UFF"){
    map<std::string,double> Vsp3;
    Vsp3["C_3"]   = 2.119;      Vsp3["S_3+2"] = 0.484;    Vsp3["Sb3"] = 1.1;
    Vsp3["N_3"]   = 0.450;      Vsp3["S_3+4"] = 0.484;    Vsp3["Te3"] = 0.3;
    Vsp3["O_3"]   = 0.018;      Vsp3["S_3+6"] = 0.484;    Vsp3["Pb3"] = 0.1;
    Vsp3["Si3"]   = 1.225;      Vsp3["Ge3"]   = 0.701;    Vsp3["Bi3"] = 1.0;
    Vsp3["P_3+3"] = 2.400;      Vsp3["As3"]   = 1.5;      Vsp3["Po3"] = 0.3;
    Vsp3["P_3+5"] = 2.400;      Vsp3["Se3"]   = 0.335;
    Vsp3["P_3+q"] = 2.400;      Vsp3["Sn3"]   = 0.199;

    map<std::string,double> Usp2;
    // I am not sure what they mean in the article - "from first through sixth periods"
    // I assume they mean from secont through sixth periods.
    Usp2["B_2"] = Usp2["C_2"] = Usp2["N_2"] = Usp2["O_2"] = 2.0;
    Usp2["S_2"] = 1.25;

    char hyb1,hyb2,hyb3,hyb4;
    std::string elt2,elt3;
    map<std::string,double>::iterator it1,it2;    
    if(ff_type1.size()>=3){ hyb1          = ff_type1.c_str()[2]; }
    else                  { hyb1          = ' ';}

    if(ff_type2.size()>=3){ hyb2          = ff_type2.c_str()[2]; }
    else                  { hyb2          = ' ';}

    if(ff_type3.size()>=3){ hyb3          = ff_type3.c_str()[2]; }
    else                  { hyb3          = ' ';}

    if(ff_type4.size()>=3){ hyb4          = ff_type4.c_str()[2]; }
    else                  { hyb4          = ' ';}

    elt2 = ff_type2.c_str()[0] + ff_type2.c_str()[1];
    elt3 = ff_type3.c_str()[0] + ff_type3.c_str()[1];

    map<std::string,double> G6sp3;
    G6sp3["O_"] = 2.0;
    G6sp3["S_"] = G6sp3["Se"] = G6sp3["Te"] = G6sp3["Po"] = 6.8;

    // Case 1: (J,K=X_3)
    if((hyb2=='3')&&(hyb3=='3')){
      // General case
      it1 = Vsp3.find(ff_type2);
      it2 = Vsp3.find(ff_type3);
      if((it1!=Vsp3.end())&&(it2!=Vsp3.end())){  Vphi = (sqrt((it1->second) * (it2->second))/9.0);  }
      else{  Vphi = (2.0/9.0);   }
      n  = 3;
      phi0 = 180.0;
      // Special case - group 6 sp3 atoms
      if( (elt2=="O_"||elt2=="S_"||elt2=="Se"||elt2=="Te"||elt2=="Po" )  &&
          (elt3=="O_"||elt3=="S_"||elt3=="Se"||elt3=="Te"||elt3=="Po" )
         ){
          Vphi = sqrt(G6sp3[elt2]*G6sp3[elt3]);
          phi0 = 90.0;
          n = 2;
      }
    }
    // Special case : (I = X_2,X_R, J = X_2,X_R, K = X_3 or L = X_2,X_R, K = X_2,X_R, J = X_3)
    else if( (((hyb1=='2')||(hyb1=='R'))&&((hyb2=='2')||(hyb2=='R'))&&(hyb3=='3')) ||
             (((hyb4=='2')||(hyb4=='R'))&&((hyb3=='2')||(hyb3=='R'))&&(hyb2=='3'))
           ){
      Vphi = (2.0/4.0); n = 3;  phi0 = 180.0;
    }

    // Case 2: (J = X_2, X_R; K = X_3 and vice versa)
    else if(  (((hyb2=='2')||(hyb2=='R'))&&(hyb3=='3'))||
              (((hyb3=='2')||(hyb3=='R'))&&(hyb2=='3'))
            ){
      // General case
      n  = 6;  Vphi  = (1.0/6.0);   phi0 = 0.0;  
      // Special case
      if(((elt2=="O_"||elt2=="S_"||elt2=="Se"||elt2=="Te"||elt2=="Po" )&& hyb2=='3') ||
         ((elt3=="O_"||elt3=="S_"||elt3=="Se"||elt3=="Te"||elt3=="Po" )&& hyb3=='3')
        ){
        it1 = Usp2.find(ff_type2);
        it2 = Usp2.find(ff_type3);
        if((it1!=Usp2.end())&&(it2!=Usp2.end())){  Vphi = 0.25*5.0*sqrt((it1->second) * (it2->second))*(1.0+4.18*log(bond_order23));  }
        else{ is_Vphi = 0; }
        phi0 = 90.0; 
        n = 2;
      }

    }

    // Case 3: (J,K=X_2)
    else if((hyb2=='2'||hyb2=='R')&&(hyb3=='2'||hyb3=='R')){
      // There are only a few types of atoms which correspond to such case
      it1 = Usp2.find(ff_type2);
      it2 = Usp2.find(ff_type3);
      if((it1!=Usp2.end())&&(it2!=Usp2.end())){  Vphi = 0.25*5.0*sqrt((it1->second) * (it2->second))*(1.0+4.18*log(bond_order23));  }
      else{ Vphi = (45.0/4.0);  }
      n = 2;
      phi0 = 180.0;
    }
    else{ is_Vphi = 0; is_phi0 = 0; is_n = 0;}

    // Convert to radians
    phi0 = deg_to_rad*phi0 ;
       

  }// UFF
  else if(ForceField_Name=="DREIDING"){

    char hyb1,hyb2,hyb3,hyb4;
    std::string elt1,elt2,elt3,elt4;
    if(ff_type1.size()>=3){ hyb1          = ff_type1.c_str()[2]; }
    else                  { hyb1          = ' ';}

    if(ff_type2.size()>=3){ hyb2          = ff_type2.c_str()[2]; }
    else                  { hyb2          = ' ';}

    if(ff_type3.size()>=3){ hyb3          = ff_type3.c_str()[2]; }
    else                  { hyb3          = ' ';}

    if(ff_type4.size()>=3){ hyb4          = ff_type4.c_str()[2]; }
    else                  { hyb4          = ' ';}

    elt1 = ff_type1.c_str()[0] + ff_type1.c_str()[1];
    elt2 = ff_type2.c_str()[0] + ff_type2.c_str()[1];
    elt3 = ff_type3.c_str()[0] + ff_type3.c_str()[1];
    elt4 = ff_type4.c_str()[0] + ff_type4.c_str()[1];

    // Case (a): (J,K=X_3)
    if((hyb2=='3')&&(hyb3=='3')){
      // General case
      Vphi = (2.0/9.0); n  = 3;  phi0 = 180.0;

      // Special case - group 6 sp3 atoms
      // Case (h)
      if( (elt2=="O_"||elt2=="S_"||elt2=="Se"||elt2=="Te"||elt2=="Po" )  &&
          (elt3=="O_"||elt3=="S_"||elt3=="Se"||elt3=="Te"||elt3=="Po" )
         ){
          Vphi = (2.0/4.0);
          phi0 = 90.0;
          n = 2;
      }
    }
    // Case (j): (I = X_2,X_R, J = X_2,X_R, K = X_3 or L = X_2,X_R, K = X_2,X_R, J = X_3)
    else if( (((hyb1=='2')||(hyb1=='R'))&&((hyb2=='2')||(hyb2=='R'))&&(hyb3=='3')) ||
             (((hyb4=='2')||(hyb4=='R'))&&((hyb3=='2')||(hyb3=='R'))&&(hyb2=='3'))
           ){
      Vphi = (2.0/4.0); n = 3;  phi0 = 180.0;
    }

    // Case (b): (J = X_2, X_R; K = X_3 and vice versa)
    else if(  (((hyb2=='2')||(hyb2=='R'))&&(hyb3=='3'))||
              (((hyb3=='2')||(hyb3=='R'))&&(hyb2=='3'))
            ){
      // General case
      Vphi  = (1.0/6.0); n= 6;  phi0 = 0.0;

      // Special case
      // Case (i)
      if(((elt2=="O_"||elt2=="S_"||elt2=="Se"||elt2=="Te"||elt2=="Po" )&& hyb2=='3') ||
         ((elt3=="O_"||elt3=="S_"||elt3=="Se"||elt3=="Te"||elt3=="Po" )&& hyb3=='3')
        ){
        Vphi = (2.0/4.0); n = 2; phi0 = 180.0;
      }
    }
    // Case (e): (J = X_2,X_R, K = X_2,X_R and bond_order23 == 1)
    else if((hyb2=='2'||hyb2=='R')&&(hyb3=='2'||hyb3=='R') && (bond_order23==1.0))
    {   Vphi = (5.0/4.0); n = 2;  phi0 = 180.0;    }
    
    // Case (f) omitted

    // Case (c): (J,K=X_2)
    else if((hyb2=='2')&&(hyb3=='2')){   Vphi = (45.0/4.0); n = 2;  phi0 = 180.0;    }

    // Case (d): (J,K=X_R)
    else if((hyb2=='R')&&(hyb3=='R')){   Vphi = (25.0/4.0); n = 2;  phi0 = 180.0;    }

    // Case (g)
    else{ is_Vphi = 0; is_phi0 = 0; is_n = 0;}

    // Convert to radians
    phi0 = deg_to_rad*phi0 ;



  }
  else{  is_Vphi = 0; is_phi0 = 0; is_n = 0;  }

}

int ForceField::get_dihedral_parameters(string ff_type1, string ff_type2, /*Inputs*/
                                        string ff_type3, string ff_type4,
                                        double bond_order12,double bond_order23, 
                                        double bond_order34,
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
      phi0 = deg_to_rad*Dihedral_Records[ff_dihedral_indx].Dihedral_phase; is_phi0 = 1;
    }
    // Potntial periodicity/multiplicity
    if(Dihedral_Records[ff_dihedral_indx].is_Dihedral_mult){
      n = Dihedral_Records[ff_dihedral_indx].Dihedral_mult; is_n = 1;
    }

  }

  //------------- Apply rules to calculate missing parameters --------------
  if(!(is_Vphi && is_n && is_phi0)){  
    dihedral_rule(ff_type1,ff_type2,ff_type3,ff_type4,bond_order12,bond_order23,bond_order34,Vphi,is_Vphi,phi0,is_phi0,n,is_n); 
  }
//  cout<<"dihedral_functional = "<<dihedral_functional<<endl;
//  cout<<"is_Vphi = "<<is_Vphi<<" is_phi0 = "<<is_phi0<<" is_n = "<<is_n<<endl;

  //------------ Additional relations --------------------------------------

  // Convert to atomic units
  Vphi *= (1.0/hartree);
  Vphi1 *= (1.0/hartree);
  Vphi2 *= (1.0/hartree);
  Vphi3 *= (1.0/hartree);



  // Assign parameters according to the force field and the potential used 
  // Do necessary scaling of the pre-computed variables
  if(dihedral_functional=="General0"||dihedral_functional=="General1"||
     dihedral_functional=="General2"||dihedral_functional=="General3"){
    prms["Vphi"] = 0.5*Vphi;
    prms["phi0"]= phi0;
    prms["n"]=n;
    status = is_Vphi * is_phi0 * is_n ;
  }
  else if(dihedral_functional=="Fourier0"||dihedral_functional=="Fourier1"){
    prms["Vphi1"] = Vphi1;
    prms["Vphi2"] = Vphi2;
    prms["Vphi3"] = Vphi3;
    status = is_Vphi1 * is_Vphi2 * is_Vphi3;
  }
  return status;
}


}// namespace libforcefield
}// namespace liblibra

