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



string ForceField::gaff_type(string elt,int geometry,string func_grp,int min_ring,int& coordination){

  //------- Transform element symbol to common format  --------------
  string type = "none"; // defaul value
  int sz = elt.size();
  coordination = 2;

  if(sz==1){   type  = toupper(elt.c_str()[0]);  type += "_";    }
  else if(sz>=2){  type  = toupper(elt.c_str()[0]); type += tolower(elt.c_str()[1]);  }

  //-----------------Type determination ------------------------------
  if(type=="H_"){ type = "hc";  // default: aliphatic hydrogen
                  if(func_grp=="Sulfhydryl"||func_grp=="Sulfonamide_H_on_S"){  type = "hs"; }
                  else if(func_grp=="Hydroxyl"||func_grp=="Carboxyl"||
                          func_grp=="Hydroperoxy"||func_grp=="Hemiketal"||
                          func_grp=="Sulfino"||func_grp=="Sulfo"||
                          func_grp=="Phosphono"||func_grp=="Phosphate"||
                          func_grp=="Phosphodiester"||func_grp=="Hemiacetal_H_on_O"){ type = "ho"; }
                  else if(func_grp=="Phosphino"){ type = "hp"; }
                  else if(func_grp=="Primary_amine"||func_grp=="Secondary_amine"||
                          func_grp=="Ammonium"||func_grp=="Sulfonamide_H_on_N"||
                          func_grp=="Primary_aldimine_H_on_N"||func_grp=="Imide_H_on_N"||
                          func_grp=="Carboxamide_H_on_N"){ type = "hn"; }
                  else if(func_grp=="Aldehyde"||func_grp=="Carboxamide_H_on_C"||
                          func_grp=="Imide_H_on_C"){ type = "ha"; }
                     }
  else if(type=="C_"){ if(func_grp=="Thione"||func_grp=="Carbonyl"||
                          func_grp=="Aldehyde"||func_grp=="Haloformyl"||
                          func_grp=="Carbonate_ester"||func_grp=="Ester"||
                          func_grp=="Carboxamide"||func_grp=="Carboxyl"||
                          func_grp=="Carboxylate"){ type = "c"; coordination = 3; }
                       else{
                         if(geometry==2){type = "c1"; coordination = 1;} // linear, sp carbon
                         else if(geometry==3){type = "c2"; coordination = 3;} // default: sp2
                         else if((geometry>=4)||(geometry==0)){type = "c3";}   // tetrahedral
                       }
                     } // type = "C_"
  else if(type=="N_"){ if(func_grp=="Carboxamide"||func_grp=="Sulfonamide"||
                          func_grp=="Imide"){ type = "n"; coordination = 3; }
                       else{
                         if(geometry==1){type = "n1"; coordination = 1;} // linear, sp nitrogen
                         else if(geometry==2){type = "n2"; coordination = 3; } // trigonal-planar, sp2 with 2 connected atoms
                         else if(geometry==3){type = "n3";  // default: Sp3 N with three connected atoms
                           if(func_grp=="Nitro"){ type = "no"; }
                         }
	  	       else if((geometry>=4)||(geometry==0)){ type = "n4";  } // Sp3 N with four connected atoms
                       }
                     } // type = "N_"
  else if(type=="O_"){ if(geometry==1){type = "o"; coordination = 1;} // Oxygen with one connected atom
                       else if(geometry==2){ type = "os"; // default: Ether and ester oxygen
                         if(func_grp=="Hydroxyl"){ type = "oh";} // Oxygen in hydroxyl group
                       }
                       else{ std::cout<<"Wrong connectivity of the oxygen atom \n"; exit(2); }
                     }// type = "O_"
  else if(type=="S_"){ if(geometry==1){ type = "s2"; }
                       else if(geometry==2){ type = "ss"; // default: Sp3 S in thio-ester and thio-ether
                         if(func_grp=="Sulfhydryl"){ type = "sh"; } // Sp3 S connected with hydrogen
                       }
                       else if(geometry==3){ type = "s4"; } // S with three connected atoms
                       else if(geometry==4){ type = "s6"; } // S with four connected atoms
                       else{ std::cout<<"Wrong connectivity of the sulfur atom \n"; exit(2); }
                     }// type = "S_"
  else if(type=="P_"){ if(geometry==2){ type = "p2"; }// Phosphate with two connected atoms
                       else if(geometry==3){ type = "p3";} // default: Phosphate with three connected atoms, such as PH3
                       else if(geometry==4){ type = "p5"; } // Phosphate with four connected atoms, such as O=P(OH)3
                       else{ std::cout<<"Wrong connectivity of the phosphorus atom \n"; exit(2); }
                     }// type = "P_"
  else if(type=="F_"){type = "f";}
  else if(type=="Cl"){type = "cl";}
  else if(type=="Br"){type = "br";}
  else if(type=="I_"){type = "i";}

  return type;
}

}// namespace libforcefield
}// namespace liblibra
