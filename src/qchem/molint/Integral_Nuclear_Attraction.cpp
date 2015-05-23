#include "Overlaps.h"
#include "Integral_Nuclear_Attraction.h"
#include "Aux_Functions.h"

namespace libqchem{
namespace libmolint{


double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                   VECTOR& Rc,int c_indx,VECTOR& DA,VECTOR& DB, VECTOR& DC,
                                   vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv,
                                   int is_normalize, int is_derivs
                                  ){

// This is an accelerated version of NAI - using pre-allocated memory - do not create and destroy arrays during run-time
// Ra, Rb - locations of primitive Gaussians
// Rc - coordinates of the nucleus with which the distribution G(a)*G(b) interacts
// DA,DB,DC - derivatives, w.r.t the positions of centers A, B and C 

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
    N = (gaussian_norm(nxa,alp_a) * gaussian_norm(nxb,alp_b));
    N *= (gaussian_norm(nya,alp_a) * gaussian_norm(nyb,alp_b));
    N *= (gaussian_norm(nza,alp_a) * gaussian_norm(nzb,alp_b));
  }


  DA = dNAI_dA * N;
  DB = dNAI_dB * N;
  DC = dNAI_dC * N; 
  
  return NAI * N;

}// nuclear_attraction_integral

}// namespace libmolint
}// namespace libqchem
