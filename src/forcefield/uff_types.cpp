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


string ForceField::uff_type(string elt,int geometry,string func_grp,int min_ring,int& coordination){

  //------- Transform element symbol to common format  --------------
  string type = "";
  int sz = elt.size();
  coordination = 2;

  if(sz==1){   type  = toupper(elt.c_str()[0]);  type += "_";    }
  else if(sz>=2){  type  = toupper(elt.c_str()[0]); type += tolower(elt.c_str()[1]);  }

  //-----------------Type determination ------------------------------
  if(type=="H_")     { if(geometry==2){type += "b"; }      }
  else if(type=="B_"){ if(geometry==3){type += "2"; coordination = 3;}
                       else {          type += "3";}
                     }
  else if(type=="C_"){ if(func_grp=="Carboxamide"||func_grp=="Carboxilate"||
                          func_grp=="Primary_aldimine"||func_grp=="Secondary_aldimine"){ type = "C_R";}
                       else{
                         if(geometry==2){type += "1"; coordination = 1;} // linear
                         else if(geometry==3){type += "2"; coordination = 3; } // trigonal-planar
                         else if((geometry>=4)||(geometry==0)){type += "3";}   // tetrahedral
                         else{ std::cout<<"Wrong connectivity of the carbon atom \n"; exit(2);  }
                       }
                     }
  else if(type=="N_"){ if(func_grp=="Carboxamide"||func_grp=="Primary_aldimine"||
                          func_grp=="Secondary_aldimine"||func_grp=="Nitro"){ type = "N_R"; }
                       else{
                         if(geometry==1){type += "1"; coordination = 1;} // linear
                         else if(geometry==2){type += "2"; coordination = 3; } // trigonal-planar
                         else if((geometry==3)||(geometry==4)||(geometry==0)){type += "3";} // tetragonal
                         else{ std::cout<<"Wrong connectivity of the nitrogen atom \n"; exit(2); }
                       }
                     }
  else if(type=="O_"){ if(func_grp=="Carboxamide"||func_grp=="Carboxilate"||func_grp=="Nitro"){ type = "O_R"; }
                       else{
                         if(geometry==1){type += "2"; coordination = 1;}
                         else if((geometry==2)||(geometry==4)||(geometry==0)){type += "3";}
                         else{ std::cout<<"Wrong connectivity of the oxygen atom \n"; exit(2);   }
                       }
                     }
  else if(type=="P_"){ if(geometry<=3){type += "3+3";}
                       else if(geometry==4){
                         if(func_grp=="Phosphono"||func_grp=="Phosphate"||func_grp=="Phosphodiester"){ type += "3+5"; }
                         else{ type +="3_q";}
                       }
                       else if(geometry==5){ type += "3+5"; } 
                     }
  else if(type=="S_"){ if(geometry==1){type += "2";coordination = 1;} // linear
                       else if((geometry==2)||(geometry==0)){type += "3+2";}
                       else if(geometry==3){ type += "3+4"; }
                       else if(geometry==4){ type += "3+6"; }
                       else{ std::cout<<"Wrong connectivity of the sulfur atom \n"; exit(2);   }
                     }

  // Relatively rare elements which have non-unique atom types
  else if(type=="Ti"){}
  else if(type=="Fe"){}
  else if(type=="Mo"){}
  else if(type=="W_"){}
  else if(type=="Re"){}

  // Most of these atom type do not require other information (type is unique) 
  //---------------------- 1 ----------------------------
  else if((type=="Ag"))     {type += "1+1"; coordination = 1;}
  else if((type=="Hg"))     {type += "1+2"; coordination = 1;}


  //---------------------- 3 -----------------------------
  else if((type=="Cd"))  {type += "3+2"; }
  else if((type=="Al")||(type=="Si")||(type=="Ge")||(type=="Sn")||(type=="Pb")) {type += "3";}
  else if((type=="Cu")) {type += "3+1";}
  else if((type=="Be")||(type=="Mg")||(type=="Zn")||(type=="Se")||(type=="Te")||(type=="Po")) {type += "3+2";}
  else if((type=="Sc")||(type=="Ga")||(type=="As")||(type=="Y_")||
          (type=="In")||(type=="Sb")||(type=="La")||(type=="Tl")||(type=="Bi")) {type += "3+3";}
  else if((type=="Zr")||(type=="Hf")) {type += "3+4";}
  else if((type=="V_")||(type=="Nb")||(type=="Ta")) {type += "3+5";}

  //---------------------- 4 --------------------
  else if((type=="Ni")||(type=="Pd")||(type=="Pt")) {type += "4+2"; coordination = 4;}
  else if((type=="Au")) {type += "4+3"; coordination = 4;}
  else if((type=="He")||(type=="Ne")||(type=="Ar")||(type=="Kr")||(type=="Xe")||(type=="Rn")) {type += "4+4"; coordination = 4;}

  //--------------------- 6 ----------------------
  else if((type=="Ca")||(type=="Mn")||(type=="Sr")||(type=="Ru")||(type=="Ba")||(type=="Ra")) {type += "6+2"; coordination = 4;}
  else if((type=="Cr")||(type=="Co")||(type=="Rh")||(type=="Ce")||(type=="Pr")||
          (type=="Nd")||(type=="Pm")||(type=="Sm")||(type=="Eu")||(type=="Gd")||
          (type=="Tb")||(type=="Dy")||(type=="Ho")||(type=="Er")||(type=="Tm")||
          (type=="Yb")||(type=="Lu")||(type=="Ir")||(type=="Ac")||(type=="Cm")||
          (type=="Bk")||(type=="Cf")||(type=="Es")||(type=="Fm")||(type=="Md")||
          (type=="No")||(type=="Lw")) {type += "6+3"; coordination = 4;}
  else if((type=="Tc")) {type += "6+5"; coordination = 4;}
  else if((type=="Th")||(type=="Pa")||(type=="U_")||(type=="Np")||(type=="Pu")||(type=="Am")){type += "6+4"; coordination = 4;}
  else if((type=="Os")) {type += "6+6"; coordination = 4;}

  //-------------------------------------------------------------------------

  return type;
}



}// namespace libforcefield
}// namespace liblibra


