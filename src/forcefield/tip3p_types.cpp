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


string ForceField::tip3p_type(string elt,int geometry,string func_grp,int min_ring,int& coordination){

  //------- Transform element symbol to common format  --------------
  string type = "";
  int sz = elt.size();
  coordination = 2;

  if(sz==1){   type  = toupper(elt.c_str()[0]);  type += "_";    }
  else if(sz>=2){  type  = toupper(elt.c_str()[0]); type += tolower(elt.c_str()[1]);  }

  //-----------------Type determination ------------------------------
  if(type=="H_")     { type = "H";   }
  else if(type=="O_"){ type = "O";   } 
  else{ type = "";}
  return type;
}


}// namespace libforcefield
}// namespace liblibra

