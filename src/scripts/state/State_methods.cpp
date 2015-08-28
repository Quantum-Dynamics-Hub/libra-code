/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "State.h"

namespace libscripts{
namespace libstate{


void State::update(){
/********************************************
 This function updated state valiables using 
 corresponding parameters and/or data from its
 sub-objects (system,thermostat,etc.)
*********************************************/

 if(syst!=NULL){
   curr_V = syst->volume(); is_curr_V = 1;
 }
}

void State::cool(){
  //------------- System -------------------
  if(syst!=NULL){    syst->cool();  }
  if(thermostat!=NULL){ thermostat->cool();   }
  if(barostat!=NULL){ barostat->cool(); }
  //------------ State variables -----------
  E_kin = 0.0;
  E_pot = 0.0;
  H0 = E_pot; is_H0 = 0;
}



}// namespace libstate
}// namespace libscripts

