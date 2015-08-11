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

