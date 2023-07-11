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

#include "Overlaps.h"
#include "Integral_Nuclear_Attraction.h"

/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{


// Auxiliary functions

void Aux_Function4(int n1,int n2,double PA,double PB,double PC,double gamma,
                   double* G,double* dGdA, double* dGdB, double* dGdC,double* f, double* dfda, double* dfdb, int n_aux){
// This is G function used in multiple summs for NAI calculations
// dGdA = dG / dPA
// dGdB = dG / dPB
// dGdC = dG / dPC - derivatives


   // dfda = df / dPA;  dfda = df / dPB
  binomial_expansion(n1,n2,PA,PB,f,dfda,dfdb,1);  // 1 = is_derivs

  int n = n1+n2;
  double p = PC; // But keep in mind: that Taketa uses p = CP = -PC, so define PC outside this function accordingly

  zero_array(G,n_aux);
  zero_array(dGdA,n_aux);
  zero_array(dGdB,n_aux);
  zero_array(dGdC,n_aux);

  
  if(n<=4){ // Explicit formula
    if(n==0){ 
      G[0] = 1.0; dGdA[0] = 0.0; dGdB[0] = 0.0; dGdC[0] = 0.0;
    }
    else if(n==1){
      G[0] = f[0];     dGdA[0] = dfda[0];  dGdB[0] = dfdb[0]; dGdC[0] = 0.0;
      G[1] = -p;       dGdA[1] = 0.0;      dGdB[1] = 0.0;     dGdC[1] = -1.0;
    }
    else if(n==2){
      double scl = (0.5/gamma);      
      G[0] = f[0] + scl*f[2];     dGdA[0] = dfda[0] + scl*dfda[2];    dGdB[0] = dfdb[0] + scl*dfdb[2];      dGdC[0] = 0.0;
      G[1] = -f[1]*p - scl*f[2];  dGdA[1] = -dfda[1]*p - scl*dfda[2]; dGdB[1] = -dfdb[1]*p - scl*dfdb[2];   dGdC[1] = -f[1];
      G[2] = p*p;                 dGdA[2] = 0.0;                      dGdB[2] = 0.0;                        dGdC[2] = 2.0*p;
    }
    else if(n==3){
      double scl = (0.5/gamma);
      double scl1 = 3.0*p*scl;
      double p2 = p*p;

      G[0] = f[0] + scl*f[2];           dGdA[0] = dfda[0] + scl*dfda[2];    dGdB[0] = dfdb[0] + scl*dfdb[2];     dGdC[0] = 0.0;
      G[1] = -f[1]*p - scl*f[2] - scl1; dGdA[1] = -dfda[1]*p - scl*dfda[2]; dGdB[1] = -dfdb[1]*p - scl*dfdb[2];  dGdC[1] = -f[1] - 3.0*scl;
      G[2] = f[2]*p2 + scl1;            dGdA[2] = dfda[2]*p2;               dGdB[2] = dfdb[2]*p2;                dGdC[2] = 2.0*p*f[2] + 3.0*scl;
      G[3] = -p2*p;                     dGdA[3] = 0.0;                      dGdB[3] = 0.0;                       dGdC[3] = -3.0*p2;
    }
    else if(n==4){
      double scl = (0.5/gamma);
      double scl1 = 3.0*p*scl;
      double scl2 = 3.0*scl*scl;
      double p2 = p*p;
      double p3 = p2*p;

      G[0] = f[0] + scl*f[2] + scl2;                    dGdA[0] = dfda[0] + scl*dfda[2];                   dGdB[0] = dfdb[0] + scl*dfdb[2];                   dGdC[0] = 0.0;
      G[1] = -f[1]*p - scl*f[2] - scl1*f[3] - 2.0*scl2; dGdA[1] = -dfda[1]*p - scl*dfda[2] - scl1*dfda[3]; dGdB[1] = -dfdb[1]*p - scl*dfdb[2] - scl1*dfdb[3]; dGdC[1] = -f[1] - f[3]*3.0*scl;
      G[2] = f[2]*p2 + scl1*f[3] + scl2 + 6.0*scl*p2;   dGdA[2] = dfda[2]*p2 + scl1*dfda[3];               dGdB[2] = dfdb[2]*p2 + scl1*dfdb[3];               dGdC[2] = 2.0*p*f[2] + 3.0*scl*f[3] + 12.0*scl*p;
      G[3] = -f[3]*p3 - 6.0*scl*p2;                     dGdA[3] = -dfda[3]*p3;                             dGdB[3] = -dfdb[3]*p3;                             dGdC[3] = -3.0*f[3]*p2 - 12.0*scl*p;
      G[4] = p2*p2;                                     dGdA[4] = 0.0;                                     dGdB[4] = 0.0;                                     dGdC[4] = 4.0*p3;
    }
  
  }
  else{
    for(int i=0;i<=n;i++){
  
      for(int r=0;r<=(i/2);r++){ // use integer division! so the upper limit is [i/2]

        for(int u=0;u<=((i-2*r)/2);u++){ // use integer division! so the upper limit is [(i-2r)/2]

          int I = (i-2*r-u);
          double den; den = FACTORIAL(r)*FACTORIAL(u)*FACTORIAL(i-2*r-2*u);
          double pref =  pow(-1.0,(i+u))*(FACTORIAL(i)/den)*FAST_POW((0.25/gamma),(r+u));
          double pw; 

          pw = FAST_POW(p,(i-2*r-2*u));
          G[I] += pref*pw*f[i];

          dGdA[I] += pref*pw*dfda[i];
          dGdB[I] += pref*pw*dfdb[i];


          if((i-2*r-2*u)==0){ pw = 0.0; } else{ pw = (i-2*r-2*u)*FAST_POW(p,(i-2*r-2*u-1)); }
          dGdC[I] += pref*pw*f[i];
            
        }// for u
      }// for r
    }// for i
  }// for n1+n2 > 4
}




double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                   VECTOR& Rc,int is_normalize, 
                                   int is_derivs, VECTOR& DA,VECTOR& DB, VECTOR& DC,
                                   vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv                                   
                                  ){
/***********************************************************************************

  < G(nxa,nya,nza,alp_a, Ra) |  1/|r-Rc|  | G(nxb,nyb,nnb,alp_b, Rb) >

// This is an accelerated version of NAI - using pre-allocated memory - do not create and destroy arrays during run-time
// Ra, Rb - locations of primitive Gaussians
// Rc - coordinates of the nucleus with which the distribution G(a)*G(b) interacts
// DA,DB,DC - derivatives, w.r.t the positions of centers A, B and C 
***********************************************************************************/

  VECTOR R_AB,P,PA,PB,PC;
  double gamma = alp_a + alp_b;
  R_AB = Ra - Rb;
  P = (alp_a*Ra + alp_b*Rb)/gamma;
  PA = P - Ra;
  PB = P - Rb;
  double sgn = -1.0; // -1 - by Taketa
  PC = sgn*(P - Rc);
  
  // Jacobian
  double dPA_dA  = (alp_a/gamma - 1.0);
  double dPA_dB  = (alp_b/gamma);
  double dPA_dC  = 0.0;

  double dPB_dA  = (alp_a/gamma);
  double dPB_dB  = (alp_b/gamma - 1.0);
  double dPB_dC  = 0.0;

  double dPC_dA  = sgn*(alp_a/gamma);
  double dPC_dB  = sgn*(alp_b/gamma);
  double dPC_dC  =  -1.0*sgn;


  // Aliases
  double *GI, *GJ, *GK;
  double *dGI_dPA, *dGJ_dPA, *dGK_dPA; // derivatives w.r.t. PA
  double *dGI_dPB, *dGJ_dPB, *dGK_dPB; // derivatives w.r.t. PB
  double *dGI_dPC, *dGJ_dPC, *dGK_dPC; // derivatives w.r.t. PC

  // any order
  GI = aux[0]; GJ = aux[1]; GK = aux[2];
  dGI_dPA = aux[3]; dGJ_dPA = aux[4]; dGK_dPA = aux[5];
  dGI_dPB = aux[6]; dGJ_dPB = aux[7]; dGK_dPB = aux[8];
  dGI_dPC = aux[9]; dGJ_dPC = aux[10]; dGK_dPC = aux[11];


  // probably need to use -PC, becuase in Taketa they use CP, not PC!!!!!!!!!!!!!!!!!  FIX THIS  !!!!
  Aux_Function4(nxa,nxb,PA.x,PB.x,PC.x,gamma,GI,dGI_dPA,dGI_dPB,dGI_dPC,aux[12],aux[13],aux[14],n_aux); 
  Aux_Function4(nya,nyb,PA.y,PB.y,PC.y,gamma,GJ,dGJ_dPA,dGJ_dPB,dGJ_dPC,aux[12],aux[13],aux[14],n_aux);
  Aux_Function4(nza,nzb,PA.z,PB.z,PC.z,gamma,GK,dGK_dPA,dGK_dPB,dGK_dPC,aux[12],aux[13],aux[14],n_aux);
                                           
  
  int max_exp = (nxa + nxb + nya + nyb + nza + nzb);
  double *C_nu;  
  VECTOR* dCdA_nu; // derivatives of C_nu with respect to A vector
  VECTOR* dCdB_nu; // derivatives of C_nu with respect to B vector
  VECTOR* dCdC_nu; // derivatives of C_nu with respect to C vector
  
  C_nu = aux[12];
  dCdA_nu = auxv[0];
  dCdB_nu = auxv[1];
  dCdC_nu = auxv[2];

  /***********************************
   If no external memory would be used
  C_nu = new double[max_exp + 1];
  dCdA_nu = new VECTOR[max_exp + 1];
  dCdB_nu = new VECTOR[max_exp + 1];
  dCdC_nu = new VECTOR[max_exp + 1];
  ***********************************/

  /// Initialize to zero
  for(int nu=0;nu<=max_exp; nu++){ C_nu[nu] = 0.0; dCdA_nu[nu] = 0.0; dCdB_nu[nu] = 0.0; dCdC_nu[nu] = 0.0; }


  for(int I=0;I<=(nxa + nxb); I++){
    for(int J=0;J<=(nya + nyb); J++){
      for(int K=0;K<=(nza + nzb); K++){

        C_nu[I+J+K] += GI[I]*GJ[I]*GK[K];

        if(is_derivs){
          // Derivatives with respect to A coordinates
          dCdA_nu[I+J+K].x += (dGI_dPA[I]*dPA_dA + dGI_dPB[I]*dPB_dA + dGI_dPC[I]*dPC_dA )*GJ[J]*GK[K];
          dCdA_nu[I+J+K].y += (dGJ_dPA[J]*dPA_dA + dGJ_dPB[J]*dPB_dA + dGJ_dPC[J]*dPC_dA )*GI[I]*GK[K];
          dCdA_nu[I+J+K].z += (dGK_dPA[K]*dPA_dA + dGK_dPB[K]*dPB_dA + dGK_dPC[K]*dPC_dA )*GI[I]*GJ[J];

          // Derivatives with respect to B coordinates
          dCdB_nu[I+J+K].x += (dGI_dPA[I]*dPA_dB + dGI_dPB[I]*dPB_dB + dGI_dPC[I]*dPC_dB )*GJ[J]*GK[K];
          dCdB_nu[I+J+K].y += (dGJ_dPA[J]*dPA_dB + dGJ_dPB[J]*dPB_dB + dGJ_dPC[J]*dPC_dB )*GI[I]*GK[K];
          dCdB_nu[I+J+K].z += (dGK_dPA[K]*dPA_dB + dGK_dPB[K]*dPB_dB + dGK_dPC[K]*dPC_dB )*GI[I]*GJ[J];

          // Derivatives with respect to C coordinates
          dCdC_nu[I+J+K].x += (dGI_dPA[I]*dPA_dC + dGI_dPB[I]*dPB_dC + dGI_dPC[I]*dPC_dC )*GJ[J]*GK[K];
          dCdC_nu[I+J+K].y += (dGJ_dPA[J]*dPA_dC + dGJ_dPB[J]*dPB_dC + dGJ_dPC[J]*dPC_dC )*GI[I]*GK[K];
          dCdC_nu[I+J+K].z += (dGK_dPA[K]*dPA_dC + dGK_dPB[K]*dPB_dC + dGK_dPC[K]*dPC_dC )*GI[I]*GJ[J];

        }// is_derivs

      }// for K
    }// for J
  }// for I

  /// Precompute inclomplete Gamma functions:
  double* F_nu;  F_nu = aux[13];
  ///  F_nu = new double[max_exp+2]; // +2 -to accomodate 1 extra nu value - for derivatives

  for(int nu=0;nu<=(max_exp+1); nu++){  F_nu[nu] = Fn(nu,gamma*PC.length2()); }


  // Now compute NAI and its derivative
  double NAI = 0.0;
  VECTOR dNAI_dA, dNAI_dB, dNAI_dC; 

  for(int nu=0;nu<=max_exp; nu++){

    NAI += C_nu[nu] * F_nu[nu];

    if(is_derivs){

      dNAI_dA += dCdA_nu[nu]*F_nu[nu] + C_nu[nu]*(-2.0*gamma*F_nu[nu+1]*PC*dPC_dA);
      dNAI_dB += dCdB_nu[nu]*F_nu[nu] + C_nu[nu]*(-2.0*gamma*F_nu[nu+1]*PC*dPC_dB);
      dNAI_dC += dCdC_nu[nu]*F_nu[nu] + C_nu[nu]*(-2.0*gamma*F_nu[nu+1]*PC*dPC_dC);

    }// is_derivs
  }// for nu


  // Include other terms
  double pref = (2.0*M_PI/gamma) * exp(-alp_a*alp_b*R_AB.length2()/gamma);

  if(is_derivs){

    VECTOR extra = pref*(-alp_a*alp_b/gamma)*2.0*NAI * R_AB;
    dNAI_dA *= pref; dNAI_dA += extra;
    dNAI_dB *= pref; dNAI_dB -= extra;
    dNAI_dC *= pref;

  }// is_derivs

  NAI = pref * NAI;// this must be done after derivatives, because we use NAI as temporary variable

  double N = 1.0;
  if(is_normalize){
    N = (gaussian_normalization_factor(nxa,alp_a) * gaussian_normalization_factor(nxb,alp_b));
    N *= (gaussian_normalization_factor(nya,alp_a) * gaussian_normalization_factor(nyb,alp_b));
    N *= (gaussian_normalization_factor(nza,alp_a) * gaussian_normalization_factor(nzb,alp_b));
  }


  DA = dNAI_dA * N;
  DB = dNAI_dB * N;
  DC = dNAI_dC * N; 
  
  return NAI * N;

}// nuclear_attraction_integral - most general version


double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                   VECTOR& Rc,int is_normalize, 
                                   int is_derivs, VECTOR& DA, VECTOR& DB, VECTOR& DC
                                  ){

  // Allocate working memory
  int i;
  int n_aux = 20;
  int n_auxv = 10;
  vector<double*> auxd(20);
  for(i=0;i<20;i++){ auxd[i] = new double[n_aux]; }
  vector<VECTOR*> auxv(5);
  for(i=0;i<5;i++){ auxv[i] = new VECTOR[n_auxv]; }

  // Do computations
  double res = nuclear_attraction_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb,
                                           Rc, is_normalize, is_derivs, DA, DB, DC, auxd, n_aux, auxv, n_auxv);

  // Clean working memory
  for(i=0;i<20;i++){ delete [] auxd[i]; }  
  auxd.clear();
  for(i=0;i<5;i++){ delete [] auxv[i]; }  
  auxv.clear();

  return res;

}// nuclear_attraction_integral - without external memory


boost::python::list nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                     int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                     VECTOR& Rc, int is_normalize, int is_derivs
                                    ){
  VECTOR dIdA, dIdB, dIdC;

  double I = nuclear_attraction_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb,
                                           Rc, is_normalize, is_derivs, dIdA, dIdB, dIdC );

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dIdA);
    res.append(dIdB);
    res.append(dIdC);
  }

  return res;
 
}

double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        VECTOR& Rc, int is_normalize
                       ){

  VECTOR dIdA, dIdB, dIdC;
  double res = nuclear_attraction_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb,
                                           Rc, is_normalize, 0, dIdA, dIdB, dIdC );
  return res;
}

double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb, VECTOR& Rc
                       ){
  double res = nuclear_attraction_integral(nxa,nya,nza,alp_a,Ra,nxb,nyb,nzb,alp_b,Rb, Rc, 1);
  return res;
}




}// namespace libmolint
}// namespace liblinalg

