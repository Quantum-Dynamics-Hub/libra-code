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

using namespace libio;

namespace libforcefield{


string ForceField::dreiding_type(string elt,int geometry,string func_grp,int min_ring,int& coordination){

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
  else if(type=="C_"){ if(geometry==2){type += "1"; coordination = 1;} // linear
                       else if(geometry==3){type += "2"; coordination = 3; } // trigonal-planar
                       else if((geometry>=4)||(geometry==0)){type += "3";}   // tetrahedral
                       else{ std::cout<<"Wrong connectivity of the carbon atom \n"; exit(2);  }
                     }
  else if(type=="N_"){ if(geometry==1){type += "1"; coordination = 1;} // linear
                       else if(geometry==2){type += "2"; coordination = 3; } // trigonal-planar
                       else if((geometry==3)||(geometry==4)||(geometry==0)){type += "3";} // tetragonal
                       else{ std::cout<<"Wrong connectivity of the nitrogen atom \n"; exit(2); }
                     }
  else if(type=="O_"){ if(geometry==1){type += "2"; coordination = 1;}
                       else if((geometry==2)||(geometry==4)||(geometry==0)){type += "3";}
                       else{ std::cout<<"Wrong connectivity of the oxygen atom \n"; exit(2);   }
                     }

  else if(type=="Al"){type += "3";}
  else if(type=="Si"){type += "3";}
  else if(type=="P_"){type += "3";}
  else if(type=="S_"){type += "3";}
  else if(type=="Ga"){type += "3";}
  else if(type=="Ge"){type += "3";}
  else if(type=="As"){type += "3";}
  else if(type=="Se"){type += "3";}
  else if(type=="In"){type += "3";}
  else if(type=="Sn"){type += "3";}
  else if(type=="Sb"){type += "3";}
  else if(type=="Te"){type += "3";}


  return type;
}


}// namespace libforcefield
}// namespace liblibra


