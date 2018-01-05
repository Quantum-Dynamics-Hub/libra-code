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
  \file ivr_prefactors.cpp
  \brief These are the C++ re-implementation of the Fortran codes from the Ananth group.
         The original codes can be found here:
         https://github.com/AnanthGroup/SC-IVR-Code-Package    

  According to original documentation:

! This file contains subroutines that 
! compute prefactors for any multidimensional  SC-IVR
!                                                                         
! MQC-IVR Prefactors:                                                     
! - Forward-Backward: General B    
! - Forward-Forward : General B    
! DHK-IVR Prefactor                                  



*/

#include "ivr.h"


/// liblibra namespace
namespace liblibra{

/// libivr namespace
namespace libivr{




complex<double> MQC_prefactor_FB_G
(vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, 
 CMATRIX& Width0,  CMATRIX& invWidth0,  CMATRIX& WidthT,  CMATRIX& invWidthT,
 CMATRIX& TuningQ, CMATRIX& invTuningQ, CMATRIX& TuningP, CMATRIX& invTuningP
){

/**
  \brief Multi-dimensional  FB-MQC-IVR Prefactor with general B
  
  \param[i]n Mfwd - Monodromy matrix array (forward) defined as Mfwd = (MqqF,MqpF,MpqF,MppF)
             each element of the array is a Ndof x Ndof matrix
  \param[i]n Mbck - Monodromy matrix array (backward) defined as Mbck = (MqqB,MqpB,MpqB,MppB)
             each element of the array is a Ndof x Ndof matrix
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of CS at t = 0
  \param[in] invWidth0 - a Ndof x Ndof matrix  inverse of Width0
  \param[in] WidthT - a Ndof x Ndof matrix with the width parameters of CS at t = T
  \param[in] invWidthT - a Ndof x Ndof matrix  inverse of WidthT
  \param[in] TuningQ - a Ndof x Ndof matrix , cq,  see theory
  \param[in] invTuningQ - a Ndof x Ndof matrix , inverse of TuningQ
  \param[in] TuningP - a Ndof x Ndof matrix , cp, see theory
  \param[in] invTuningP - a Ndof x Ndof matrix , inverse of TuningP

      

*/

  int Ndof = Mfwd[0].n_rows;
  complex<double> one(0.0, 1.0);


  CMATRIX Id(Ndof, Ndof), G(Ndof, Ndof), invG(Ndof, Ndof);
  G = (TuningQ + WidthT) * TuningP + TuningQ * (TuningP + invWidthT);

  for(int i=0;i<Ndof;i++){
    Id.set(i,i, 1.0, 0.0);
    invG.set(i,i, 1.0/G.get(i,i));
  }


  CMATRIX A(Ndof, Ndof), B(Ndof, Ndof), C(Ndof, Ndof), D(Ndof, Ndof);
  A = Mbck[3] - one * Width0 * Mbck[1];
  B = Mfwd[3] * Width0 + one * Mfwd[2];
  C = Width0 * Mbck[0] + one * Mbck[2];
  D = Mfwd[0] - one * Mfwd[1] * Width0;


  CMATRIX G1(Ndof,Ndof), Gp(Ndof,Ndof), Gq(Ndof,Ndof);  
  G1   = invG + Id;
  Gq = (0.5 * Width0 + TuningQ) * invG;
  Gp = (0.5 * invWidthT + TuningP) * invG; 


  CMATRIX coeff(Ndof, Ndof);
  coeff = 0.5* invWidth0 * G * ( A * (0.5 * G1 * B + Gq * D) + C * (0.5 * G1 * D + Gp * B));


  return det(coeff);

}




complex<double> MQC_prefactor_FF_G
(vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, 
 CMATRIX& Width0,  CMATRIX& invWidth0,  CMATRIX& WidthT,  CMATRIX& invWidthT,
 CMATRIX& TuningQ, CMATRIX& invTuningQ, CMATRIX& TuningP, CMATRIX& invTuningP
){

/**
  \brief Multi-dimensional  FF-MQC-IVR Prefactor with general B
  
  \param[i]n Mfwd - Monodromy matrix array (forward) defined as Mfwd = (MqqF,MqpF,MpqF,MppF)
             each element of the array is a Ndof x Ndof matrix
  \param[i]n Mbck - Monodromy matrix array (backward) defined as Mbck = (MqqB,MqpB,MpqB,MppB)
             each element of the array is a Ndof x Ndof matrix
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of CS at t = 0
  \param[in] invWidth0 - a Ndof x Ndof matrix  inverse of Width0
  \param[in] WidthT - a Ndof x Ndof matrix with the width parameters of CS at t = T
  \param[in] invWidthT - a Ndof x Ndof matrix  inverse of WidthT
  \param[in] TuningQ - a Ndof x Ndof matrix , cq,  see theory
  \param[in] invTuningQ - a Ndof x Ndof matrix , inverse of TuningQ
  \param[in] TuningP - a Ndof x Ndof matrix , cp, see theory
  \param[in] invTuningP - a Ndof x Ndof matrix , inverse of TuningP

      

*/

  int Ndof = Mfwd[0].n_rows;
  complex<double> one(0.0, 1.0);

  CMATRIX A(Ndof, Ndof), B(Ndof, Ndof), C(Ndof, Ndof), D(Ndof, Ndof);

  A = Mfwd[3] - one * WidthT * Mfwd[1];
  B = Mbck[3] * WidthT + one * Mbck[2];
  C = WidthT * Mfwd[0] + one * Mfwd[2];
  D = Mbck[0] - one * Mbck[1] * WidthT;

  CMATRIX G(Ndof, Ndof), invG(Ndof, Ndof), Id(Ndof, Ndof);
  G = TuningP * (TuningQ + Width0) + TuningQ * (TuningP + invWidth0);

  for(int i=0;i<Ndof;i++){
    Id.set(i,i, 1.0, 0.0);
    invG.set(i,i, 1.0/G.get(i,i));
  }

  CMATRIX Gq(Ndof, Ndof), Gp(Ndof, Ndof);
  Gq = (0.5 * Width0+TuningQ) * invG;
  Gp = (0.5 * invWidth0+TuningP) * invG; 


  CMATRIX coeff(Ndof, Ndof);
  coeff = 0.5 * invWidthT * G * (0.5 * A * (invG + Id) * B + C * Gp * B + 0.5 * C * (invG + Id) * D + A * Gq * D);


  return det(coeff);

}





vector<complex<double> > DHK_prefactor
(vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, 
 CMATRIX& Width0,  CMATRIX& invWidth0,  CMATRIX& WidthT,  CMATRIX& invWidthT,
 CMATRIX& TuningQ, CMATRIX& invTuningQ, CMATRIX& TuningP, CMATRIX& invTuningP
){

/**
  \brief Multi-dimensional  Herman-Kluk Prefactors for DHK TCF
  
  \param[i]n Mfwd - Monodromy matrix array (forward) defined as Mfwd = (MqqF,MqpF,MpqF,MppF)
             each element of the array is a Ndof x Ndof matrix
  \param[i]n Mbck - Monodromy matrix array (backward) defined as Mbck = (MqqB,MqpB,MpqB,MppB)
             each element of the array is a Ndof x Ndof matrix
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of CS at t = 0
  \param[in] invWidth0 - a Ndof x Ndof matrix  inverse of Width0
  \param[in] WidthT - a Ndof x Ndof matrix with the width parameters of CS at t = T
  \param[in] invWidthT - a Ndof x Ndof matrix  inverse of WidthT
  \param[in] TuningQ - a Ndof x Ndof matrix , cq,  see theory
  \param[in] invTuningQ - a Ndof x Ndof matrix , inverse of TuningQ
  \param[in] TuningP - a Ndof x Ndof matrix , cp, see theory
  \param[in] invTuningP - a Ndof x Ndof matrix , inverse of TuningP

      

*/

  int Ndof = Mfwd[0].n_rows;
  complex<double> one(0.0, 1.0);

  CMATRIX sqrt_Width0(Ndof, Ndof), sqrt_invWidth0(Ndof, Ndof); 
  CMATRIX sqrt_WidthT(Ndof, Ndof), sqrt_invWidthT(Ndof, Ndof);

  for(int i=0;i<Ndof;i++){
    sqrt_Width0.set(i,i, std::sqrt(Width0.get(i,i)));
    sqrt_invWidth0.set(i,i, std::sqrt(invWidth0.get(i,i)));
    sqrt_WidthT.set(i,i, std::sqrt(WidthT.get(i,i)));
    sqrt_invWidthT.set(i,i, std::sqrt(invWidthT.get(i,i)));
  }


  vector<complex<double> > res(2, complex<double>(0.0, 0.0));

  CMATRIX A(Ndof, Ndof), B(Ndof, Ndof), C(Ndof, Ndof), D(Ndof, Ndof);

  A = sqrt_Width0 * Mfwd[0] * sqrt_invWidthT;
  B = sqrt_invWidth0 * Mfwd[3] * sqrt_WidthT;
  C = sqrt_Width0 * Mfwd[1] * sqrt_WidthT;
  D = sqrt_invWidth0 * Mfwd[2] * sqrt_invWidthT;

  CMATRIX tmp(Ndof, Ndof);
  tmp = 0.5*(A + B + one*(D - C));
  res[0] =  det(tmp);

  A = sqrt_WidthT * Mbck[0] * sqrt_invWidth0;
  B = sqrt_invWidthT * Mbck[3] * sqrt_Width0;
  C = sqrt_WidthT * Mbck[1] * sqrt_Width0;
  D = sqrt_invWidthT * Mbck[2] * sqrt_invWidth0;

  tmp = 0.5*(A + B + one*(D - C));
  res[1] =  det(tmp);


  return res;

}


}/// namespace libivr
}/// liblibra
