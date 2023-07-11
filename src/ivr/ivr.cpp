/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file ivr.cpp
  \brief This file implements some of the auxiliary classes for 
         storing parameters and dynamical date used in the
         IVR procedures

*/

#include "ivr.h"


/// liblibra namespace
namespace liblibra{

/// libivr namespace
namespace libivr{


ivr_params::ivr_params(int ndof_){

  Ndof = ndof_; 

  qIn = new MATRIX(Ndof, 1);
  pIn = new MATRIX(Ndof, 1);

  Width0 = new MATRIX(Ndof, Ndof);  
  WidthT = new MATRIX(Ndof, Ndof);
  invWidth0 = new MATRIX(Ndof, Ndof);
  invWidthT = new MATRIX(Ndof, Ndof);
                       
  TuningQ = new MATRIX(Ndof, Ndof);
  TuningP = new MATRIX(Ndof, Ndof);
  invTuningQ = new MATRIX(Ndof, Ndof);
  invTuningP = new MATRIX(Ndof, Ndof);

}

ivr_params::ivr_params(const ivr_params& prms){

  Ndof = prms.Ndof; 

  qIn = new MATRIX(Ndof, 1);   *qIn = *prms.qIn;
  pIn = new MATRIX(Ndof, 1);   *pIn = *prms.pIn;

  Width0 = new MATRIX(Ndof, Ndof);     *Width0 = *prms.Width0;
  WidthT = new MATRIX(Ndof, Ndof);     *WidthT = *prms.WidthT;
  invWidth0 = new MATRIX(Ndof, Ndof);  *invWidth0 = *prms.invWidth0;
  invWidthT = new MATRIX(Ndof, Ndof);  *invWidthT = *prms.invWidthT;
                       
  TuningQ = new MATRIX(Ndof, Ndof);    *TuningQ = *prms.TuningQ;
  TuningP = new MATRIX(Ndof, Ndof);    *TuningP = *prms.TuningP;
  invTuningQ = new MATRIX(Ndof, Ndof); *invTuningQ = *prms.invTuningQ;
  invTuningP = new MATRIX(Ndof, Ndof); *invTuningP = *prms.invTuningP;

}

ivr_params::~ivr_params(){

  delete qIn;    qIn = NULL;
  delete pIn;    pIn = NULL;

  delete Width0; Width0 = NULL; 
  delete WidthT; WidthT = NULL; 
  delete invWidth0; invWidth0 = NULL; 
  delete invWidthT; invWidthT = NULL; 

  delete TuningQ; TuningQ = NULL; 
  delete TuningP; TuningP = NULL; 
  delete invTuningQ; invTuningQ = NULL; 
  delete invTuningP; invTuningP = NULL; 

}

void ivr_params::set_qIn(MATRIX& qIn_){  *qIn = qIn_; }
void ivr_params::set_pIn(MATRIX& pIn_){  *pIn = pIn_; }

void ivr_params::set_Width0(double w){

  for(int i=0; i<Ndof; i++){  
    Width0->set(i,i, w); 
    invWidth0->set(i,i, 1.0/w); 
  }
}

void ivr_params::set_WidthT(double w){

  for(int i=0; i<Ndof; i++){  
    WidthT->set(i,i, w); 
    invWidthT->set(i,i, 1.0/w); 
  }
}

void ivr_params::set_TuningQ(double cq){

  for(int i=0; i<Ndof; i++){  
    TuningQ->set(i,i, cq); 
    invTuningQ->set(i,i, 1.0/cq); 
  }
}

void ivr_params::set_TuningP(double cp){

  for(int i=0; i<Ndof; i++){  
    TuningP->set(i,i, cp); 
    invTuningP->set(i,i, 1.0/cp); 
  }
}



}/// namespace libivr
}/// namespace liblibra


