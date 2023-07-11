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


string ForceField::tripos_type(string elt,int geometry,string func_grp,int min_ring,int& coordination){

  //------- Transform element symbol to common format  --------------
  string type = "";
  int sz = elt.size();
  coordination = 2;

  if(sz==1){   type  = toupper(elt.c_str()[0]);  type += "_";    }
  else if(sz>=2){  type  = toupper(elt.c_str()[0]); type += tolower(elt.c_str()[1]);  }

  //-----------------Type determination ------------------------------
  if(type=="H_")     { type = "H";   }
  else if(type=="C_"){ if(geometry==2){type = "C.1"; coordination = 1;} // linear
                       else if(geometry==3){type = "C.2"; coordination = 3; } // trigonal-planar
                       else if((geometry>=4)||(geometry==0)){type = "C.3";}   // tetrahedral
                       else{ std::cout<<"Wrong connectivity of the carbon atom \n"; exit(2);  }
                     }
  else if(type=="N_"){ if(func_grp=="Carboxamide"||func_grp=="Sulfonamide"||func_grp=="Imide"){
                         type = "N.am"; coordination = 3;
                       }
                       else{
                         if(geometry==1){type = "N.1"; coordination = 1;} // linear
                         else if(geometry==2){type = "N.2"; coordination = 3; } // trigonal-planar
                         else if((geometry==3)||(geometry==0)){type += "3";} // tetragonal
                         else if(geometry==4){ type = "N.4"; }
                         else{ std::cout<<"Wrong connectivity of the nitrogen atom \n"; exit(2); }
                       }
                     }// N.3, N.ar, N.pl3 are not included yet
  else if(type=="O_"){ if(geometry==1){type = "O.2"; coordination = 1;}
                       else if((geometry==2)||(geometry==4)||(geometry==0)){type = "O.3";}
                       else{ std::cout<<"Wrong connectivity of the oxygen atom \n"; exit(2);   }
                     }
  else if(type=="S_"){ if(func_grp=="Sulfinyl"||func_grp=="Sulfino"){ type = "S.o"; }
                       else if(func_grp=="Sulfonyl"||func_grp=="Sulfo"||func_grp=="Sulfonamide"){ type = "S.o2"; }
                       else{
                         if(geometry==1){type = "S.2"; coordination = 1;}
                         else if((geometry==2)||(geometry==4)||(geometry==0)){type = "S.3";}
                         else{ std::cout<<"Wrong connectivity of the oxygen atom \n"; exit(2);   }
                       }
                     }

  else if(type=="Li"){type = "Li";}
  else if(type=="Na"){type = "Na";}
  else if(type=="K_"){type = "K";}
  else if(type=="F_"){type = "F";}
  else if(type=="Cl"){type = "Cl";}
  else if(type=="Br"){type = "Br";}
  else if(type=="I_"){type = "I";}
  else if(type=="Ca"){type = "Ca";}
  else if(type=="Al"){type = "Al";}
  else if(type=="Si"){type = "Si";}
  else if(type=="P_"){type = "P.3";}

  return type;
}


}// namespace libforcefield
}// namespace liblibra

