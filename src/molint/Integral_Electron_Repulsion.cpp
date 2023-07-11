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
#include "Integral_Electron_Repulsion.h"


/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;
using namespace liblinalg;

namespace libmolint{


//void Aux_Function5

void eri_aux1
( int l1,int l2,double a,double b,double gamma,
  double* H,double* dHda, double* dHdb,
  double* f, double* dfda, double* dfdb, int n_aux
){

  int maxL = (l1+l2); // this is correct
  int L,i,r,l;

  binomial_expansion(l1,l2,a,b,f,dfda,dfdb, 1); // 1 - is_derivs

  zero_array(H,n_aux);
  zero_array(dHda,n_aux);
  zero_array(dHdb,n_aux);

  
  for(i=0;i<=maxL;i++){
    for(r=0;r<=(i/2);r++){   //!!! Integer division
      L = i - 2*r;    
      if(L>=0){
  
        double prefact = ( (FACTORIAL(i)/(1.0*FACTORIAL(r)*FACTORIAL(L))) ) * FAST_POW((0.25/gamma),(i-r));
        H[L] += prefact*f[i];  
        dHda[L] += prefact*dfda[i];
        dHdb[L] += prefact*dfdb[i];         
 
      }// L>0
    }// for r
  }// for i

}



//void Aux_Function6
void eri_aux2
( int l1,int l2,int l3,int l4,double PA,double PB,double QC,double QD,double p,double gamma1,double gamma2,
  double* G, double* dGdA, double* dGdB, double* dGdC, double* dGdD, double* dGdp, 
  double* HL,double* dHLda, double* dHLdb,
  double* HM,double* dHMda, double* dHMdb,
  double* f, double* dfda, double* dfdb,
  int n_aux
){

// This is C(G) function used in multiple summs for ERI calculations
// dGdA = dG / dPA
// dGdB = dG / dPB
// dGdC = dG / dQC 
// dGdD = dG / dQD
// dGdp = dG / dp


  zero_array(G,n_aux);
  zero_array(dGdA,n_aux);
  zero_array(dGdB,n_aux);
  zero_array(dGdC,n_aux);
  zero_array(dGdD,n_aux);
  zero_array(dGdp,n_aux);


// This is C function used in multiple summs for ERI calculations
  int maxL = l1 + l2;
  int maxM = l3 + l4;

  eri_aux1(l1,l2,PA,PB,gamma1,HL,dHLda,dHLdb,f,dfda,dfdb,n_aux); // HL of size l1+l2
  eri_aux1(l3,l4,QC,QD,gamma2,HM,dHMda,dHMdb,f,dfda,dfdb,n_aux); // HM of size l3+l4


  double d = 0.25*((1.0/gamma1) + (1.0/gamma2));

  int L,M,u;
  for(L=0;L<=maxL;L++){
    for(M=0;M<=maxM;M++){
      int maxU = ((L+M)/2);  //INTEGER_DIVISION!!! "largest integer less of equal (L+M)/2"
      for(u=0;u<=maxU;u++){
        int I = L + M - u;
        if(I>=0){
          double prefac = FAST_POW(-1.0,(M+u)) * FACTORIAL(L+M)/ (1.0*FACTORIAL(u) * FACTORIAL(L+M-2*u) * FAST_POW(d,(L+M-u)));
          double fp = FAST_POW(p,(L+M-2*u));

         
          G[I] += prefac * HL[L] * HM[M] * fp;
          dGdA[I] += prefac * dHLda[L] * HM[M] * fp;
          dGdB[I] += prefac * dHLdb[L] * HM[M] * fp;
          dGdC[I] += prefac * HL[L] * dHMda[M] * fp;
          dGdD[I] += prefac * HL[L] * dHMdb[M] * fp;

          if((L+M-2*u)>0){ // to avoid singularity when p = 0
            dGdp[I] += (L+M-2*u) * prefac * HL[L] * HM[M] * FAST_POW(p,(L+M-2*u-1));
          }

        }
      }// for u
    }// for M
  }// for L

}// eri_aux2





double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize, 
   int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD,
    vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
){

/// This is equivalent to chemists' notation: (ab|cd) 
///                 or
///                    to physicists notation: <ac|bd>


    VECTOR R_AB,R_CD,P,Q,PA,PB,QC,QD,PQ;
    R_AB = Ra - Rb;
    R_CD = Rc - Rd;

    int maxI = nxa + nxb + nxc + nxd;
    int maxJ = nya + nyb + nyc + nyd;
    int maxK = nza + nzb + nzc + nzd;

    double gamma1,gamma2;
    gamma1 = alp_a + alp_b;
    gamma2 = alp_c + alp_d;


    P = (alp_a*Ra + alp_b*Rb)/gamma1;
    Q = (alp_c*Rc + alp_d*Rd)/gamma2;

    PA = P - Ra;
    PB = P - Rb;
    QC = Q - Rc;
    QD = Q - Rd;

    PQ = Q - P;


    // Jacobian
    double dPA_dA  = (alp_a/gamma1 - 1.0);
    double dPA_dB  = (alp_b/gamma1);
    double dPA_dC  = 0.0;
    double dPA_dD  = 0.0;

    double dPB_dA  = (alp_a/gamma1);
    double dPB_dB  = (alp_b/gamma1 - 1.0);
    double dPB_dC  = 0.0;
    double dPB_dD  = 0.0;

    double dQC_dA  = 0.0;
    double dQC_dB  = 0.0;
    double dQC_dC  = (alp_c/gamma2 - 1.0);
    double dQC_dD  = (alp_d/gamma2);

    double dQD_dA  = 0.0;
    double dQD_dB  = 0.0;
    double dQD_dC  = (alp_c/gamma2);
    double dQD_dD  = (alp_d/gamma2 - 1.0);

    double dPQ_dA  = -(alp_a/gamma1);
    double dPQ_dB  = -(alp_b/gamma1);
    double dPQ_dC  = (alp_c/gamma2);
    double dPQ_dD  = (alp_d/gamma2);



    // Aliases
    double *GI, *GJ, *GK;
    double *dGI_dPA, *dGJ_dPA, *dGK_dPA; // derivatives w.r.t. PA
    double *dGI_dPB, *dGJ_dPB, *dGK_dPB; // derivatives w.r.t. PB
    double *dGI_dQC, *dGJ_dQC, *dGK_dQC; // derivatives w.r.t. QC
    double *dGI_dQD, *dGJ_dQD, *dGK_dQD; // derivatives w.r.t. QD
    double *dGI_dPQ, *dGJ_dPQ, *dGK_dPQ; // derivatives w.r.t. PQ

    // any order
         GI = aux[0];       GJ = aux[1];       GK = aux[2];
    dGI_dPA = aux[3];  dGJ_dPA = aux[4];  dGK_dPA = aux[5];
    dGI_dPB = aux[6];  dGJ_dPB = aux[7];  dGK_dPB = aux[8];
    dGI_dQC = aux[9];  dGJ_dQC = aux[10]; dGK_dQC = aux[11];
    dGI_dQD = aux[12]; dGJ_dQD = aux[13]; dGK_dQD = aux[14];
    dGI_dPQ = aux[15]; dGJ_dPQ = aux[16]; dGK_dPQ = aux[17];


//    // Testing derivatives:
//    double dx = 0.0001;
//    Aux_Function6(GA.x_exp,GB.x_exp,GC.x_exp,GD.x_exp, PA.x+dx,PB.x,QC.x,QD.x,PQ.x, gamma1,gamma2,  GI, dGI_dPA, dGI_dPB, dGI_dQC, dGI_dQD, dGI_dPQ, aux[18],aux[19],aux[20], aux[21],aux[22],aux[23], aux[24],aux[25],aux[26], n_aux);
//    Aux_Function6(GA.x_exp,GB.x_exp,GC.x_exp,GD.x_exp, PA.x-dx,PB.x,QC.x,QD.x,PQ.x, gamma1,gamma2,  GJ, dGJ_dPA, dGI_dPB, dGI_dQC, dGI_dQD, dGI_dPQ, aux[18],aux[19],aux[20], aux[21],aux[22],aux[23], aux[24],aux[25],aux[26], n_aux);
//
//    for(int i=0;i<=maxI; i++){
//      cout<<i<<" numer()="<<(GJ[i]-GI[i])/(2.0*dx)<<" anal= "<<dGI_dPA[i]<<endl;
//    }


    eri_aux2(nxa, nxb, nxc, nxd, PA.x,PB.x,QC.x,QD.x,PQ.x, gamma1,gamma2,  GI, dGI_dPA, dGI_dPB, dGI_dQC, dGI_dQD, dGI_dPQ, aux[18],aux[19],aux[20], aux[21],aux[22],aux[23], aux[24],aux[25],aux[26], n_aux);
    eri_aux2(nya, nyb, nyc, nyd, PA.y,PB.y,QC.y,QD.y,PQ.y, gamma1,gamma2,  GJ, dGJ_dPA, dGJ_dPB, dGJ_dQC, dGJ_dQD, dGJ_dPQ, aux[18],aux[19],aux[20], aux[21],aux[22],aux[23], aux[24],aux[25],aux[26], n_aux);
    eri_aux2(nza, nzb, nzc, nzd, PA.z,PB.z,QC.z,QD.z,PQ.z, gamma1,gamma2,  GK, dGK_dPA, dGK_dPB, dGK_dQC, dGK_dQD, dGK_dPQ, aux[18],aux[19],aux[20], aux[21],aux[22],aux[23], aux[24],aux[25],aux[26], n_aux);



    double *C_nu;  C_nu = aux[27];

    VECTOR* dCdA_nu; // derivatives of C_nu with respect to A vector
    VECTOR* dCdB_nu; // derivatives of C_nu with respect to B vector
    VECTOR* dCdC_nu; // derivatives of C_nu with respect to C vector
    VECTOR* dCdD_nu; // derivatives of C_nu with respect to D vector
    dCdA_nu = auxv[0];
    dCdB_nu = auxv[1];
    dCdC_nu = auxv[2];
    dCdD_nu = auxv[3];


    int max_exp = maxI + maxJ + maxK;
    for(int nu=0;nu<=max_exp; nu++){ C_nu[nu] = 0.0; dCdA_nu[nu] = 0.0; dCdB_nu[nu] = 0.0; dCdC_nu[nu] = 0.0; dCdD_nu[nu] = 0.0; }


    for(int I=0;I<=maxI; I++){
      for(int J=0;J<=maxJ; J++){
        for(int K=0;K<=maxK; K++){
          double tmp = GI[I]*GJ[J]*GK[K];
 
          C_nu[I+J+K] += tmp;

          if(is_derivs){
            // Derivatives with respect to A coordinates
            dCdA_nu[I+J+K].x += ( dGI_dPA[I]*dPA_dA + dGI_dPB[I]*dPB_dA + dGI_dQC[I]*dQC_dA + dGI_dQD[I]*dQD_dA + dGI_dPQ[I]*dPQ_dA)*GJ[J]*GK[K];
            dCdA_nu[I+J+K].y += GI[I]*( dGJ_dPA[J]*dPA_dA + dGJ_dPB[J]*dPB_dA + dGJ_dQC[J]*dQC_dA + dGJ_dQD[J]*dQD_dA + dGJ_dPQ[J]*dPQ_dA)*GK[K];
            dCdA_nu[I+J+K].z += GI[I]*GJ[J]*( dGK_dPA[K]*dPA_dA + dGK_dPB[K]*dPB_dA + dGK_dQC[K]*dQC_dA + dGK_dQD[K]*dQD_dA + dGK_dPQ[K]*dPQ_dA);

            // Derivatives with respect to B coordinates
            dCdB_nu[I+J+K].x += ( dGI_dPA[I]*dPA_dB + dGI_dPB[I]*dPB_dB + dGI_dQC[I]*dQC_dB + dGI_dQD[I]*dQD_dB + dGI_dPQ[I]*dPQ_dB)*GJ[J]*GK[K];
            dCdB_nu[I+J+K].y += GI[I]*( dGJ_dPA[J]*dPA_dB + dGJ_dPB[J]*dPB_dB + dGJ_dQC[J]*dQC_dB + dGJ_dQD[J]*dQD_dB + dGJ_dPQ[J]*dPQ_dB)*GK[K];
            dCdB_nu[I+J+K].z += GI[I]*GJ[J]*( dGK_dPA[K]*dPA_dB + dGK_dPB[K]*dPB_dB + dGK_dQC[K]*dQC_dB + dGK_dQD[K]*dQD_dB + dGK_dPQ[K]*dPQ_dB);

            // Derivatives with respect to C coordinates
            dCdC_nu[I+J+K].x += ( dGI_dPA[I]*dPA_dC + dGI_dPB[I]*dPB_dC + dGI_dQC[I]*dQC_dC + dGI_dQD[I]*dQD_dC + dGI_dPQ[I]*dPQ_dC)*GJ[J]*GK[K];
            dCdC_nu[I+J+K].y += GI[I]*( dGJ_dPA[J]*dPA_dC + dGJ_dPB[J]*dPB_dC + dGJ_dQC[J]*dQC_dC + dGJ_dQD[J]*dQD_dC + dGJ_dPQ[J]*dPQ_dC)*GK[K];
            dCdC_nu[I+J+K].z += GI[I]*GJ[J]*( dGK_dPA[K]*dPA_dC + dGK_dPB[K]*dPB_dC + dGK_dQC[K]*dQC_dC + dGK_dQD[K]*dQD_dC + dGK_dPQ[K]*dPQ_dC);

            // Derivatives with respect to D coordinates
            dCdD_nu[I+J+K].x += ( dGI_dPA[I]*dPA_dD + dGI_dPB[I]*dPB_dD + dGI_dQC[I]*dQC_dD + dGI_dQD[I]*dQD_dD + dGI_dPQ[I]*dPQ_dD)*GJ[J]*GK[K];
            dCdD_nu[I+J+K].y += GI[I]*( dGJ_dPA[J]*dPA_dD + dGJ_dPB[J]*dPB_dD + dGJ_dQC[J]*dQC_dD + dGJ_dQD[J]*dQD_dD + dGJ_dPQ[J]*dPQ_dD)*GK[K];
            dCdD_nu[I+J+K].z += GI[I]*GJ[J]*( dGK_dPA[K]*dPA_dD + dGK_dPB[K]*dPB_dD + dGK_dQC[K]*dQC_dD + dGK_dQD[K]*dQD_dD + dGK_dPQ[K]*dPQ_dD);

          }// if is_derivs

        }// for K
      }// for J
    }// for I



    // Precompute inclomplete Gamma functions:
    double* F_nu;  F_nu = aux[28];
    double d4 = ((1.0/gamma1) + (1.0/gamma2));
    for(int nu=0;nu<=(maxI+maxJ+maxK+1); nu++){
        F_nu[nu] = Fn(nu,PQ.length2()/d4); // Aux_Function2(nu,PQ.length2()/d4);
    }// for nu



    // Now compute ERI and its derivative
    double ERI = 0.0;
    VECTOR dERI_dA, dERI_dB, dERI_dC, dERI_dD;
    dERI_dA = 0.0;
    dERI_dB = 0.0;
    dERI_dC = 0.0;
    dERI_dD = 0.0;

    for(int nu=0;nu<=(maxI+maxJ+maxK); nu++){

      ERI += C_nu[nu] * F_nu[nu];

      if(is_derivs){
        dERI_dA += ( dCdA_nu[nu]*F_nu[nu] + C_nu[nu]*(-(2.0/d4)*F_nu[nu+1]*PQ*dPQ_dA) );
        dERI_dB += ( dCdB_nu[nu]*F_nu[nu] + C_nu[nu]*(-(2.0/d4)*F_nu[nu+1]*PQ*dPQ_dB) );
        dERI_dC += ( dCdC_nu[nu]*F_nu[nu] + C_nu[nu]*(-(2.0/d4)*F_nu[nu+1]*PQ*dPQ_dC) );
        dERI_dD += ( dCdD_nu[nu]*F_nu[nu] + C_nu[nu]*(-(2.0/d4)*F_nu[nu+1]*PQ*dPQ_dD) );
      }

    }// for nu


    double N = 1.0;
    if(is_normalize){

      N = (gaussian_normalization_factor(nxa,alp_a) * gaussian_normalization_factor(nxb,alp_b));
      N *= (gaussian_normalization_factor(nya,alp_a) * gaussian_normalization_factor(nyb,alp_b));
      N *= (gaussian_normalization_factor(nza,alp_a) * gaussian_normalization_factor(nzb,alp_b));

      N *= (gaussian_normalization_factor(nxa,alp_c) * gaussian_normalization_factor(nxb,alp_d));
      N *= (gaussian_normalization_factor(nya,alp_c) * gaussian_normalization_factor(nyb,alp_d));
      N *= (gaussian_normalization_factor(nza,alp_c) * gaussian_normalization_factor(nzb,alp_d));

    }
    


    double pref0 = N*(2.0*M_PI*M_PI/(gamma1*gamma2))*sqrt(M_PI/(gamma1+gamma2));
    double pref_AB = exp(-alp_a*alp_b*R_AB.length2()/gamma1);
    double pref_CD = exp(-alp_c*alp_d*R_CD.length2()/gamma2);


    ERI = pref0 * pref_AB * pref_CD * ERI;

    DA = 0.0; DB = 0.0; DC = 0.0; DD = 0.0;
    if(is_derivs){
      DA = pref0 * pref_AB * pref_CD * ( dERI_dA + ERI*(-2.0*alp_a*alp_b/gamma1)*R_AB );   
      DB = pref0 * pref_AB * pref_CD * ( dERI_dB + ERI*( 2.0*alp_a*alp_b/gamma1)*R_AB );
      DC = pref0 * pref_AB * pref_CD * ( dERI_dC + ERI*(-2.0*alp_c*alp_d/gamma2)*R_CD );
      DD = pref0 * pref_AB * pref_CD * ( dERI_dD + ERI*( 2.0*alp_c*alp_d/gamma2)*R_CD );
    }
    

    return ERI;

}// eri - the very general version


double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize, 
   int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD
){

  // Allocate working memory
  int i;
  int n_aux = 40;
  int n_auxv = 40;
  vector<double*> auxd(30);
  for(i=0;i<30;i++){ auxd[i] = new double[n_aux]; }
  vector<VECTOR*> auxv(5);
  for(i=0;i<5;i++){ auxv[i] = new VECTOR[n_auxv]; }

  // Do computations
  double res = electron_repulsion_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb,
                                           nxc,nyc,nzc,alp_c,Rc, nxd,nyd,nzd,alp_d,Rd,
                                           is_normalize, is_derivs, DA, DB, DC, DD,
                                           auxd, n_aux, auxv, n_auxv);
  // Clean working memory
  for(i=0;i<30;i++){ delete [] auxd[i]; }  
  auxd.clear();
  for(i=0;i<5;i++){ delete [] auxv[i]; }  
  auxv.clear();

  return res;

}// eri - without externam memory allocation



boost::python::list electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize, int is_derivs
){

  VECTOR DA,DB,DC,DD;

  double I = electron_repulsion_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb,
                                         nxc,nyc,nzc,alp_c,Rc, nxd,nyd,nzd,alp_d,Rd,
                                         is_normalize, is_derivs, DA, DB, DC, DD     );
  boost::python::list res;
  res.append(I);
 
  if(is_derivs){
    res.append(DA);    res.append(DB);    res.append(DC);    res.append(DD);
  }

  return res;
 
}// eri - python-exported version


double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize
){

  VECTOR DA, DB, DC, DD;

  // Do computations
  double res = electron_repulsion_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb,
                                           nxc,nyc,nzc,alp_c,Rc, nxd,nyd,nzd,alp_d,Rd,
                                           is_normalize, 0, DA, DB, DC, DD
                                          );
  return res;

}// eri - no derivatives

double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd
){

  double res = electron_repulsion_integral(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb,
                                           nxc,nyc,nzc,alp_c,Rc, nxd,nyd,nzd,alp_d,Rd,  1   );
  return res;

}// eri - no derivatives, normalization



}// namespace libqchem
}// namespace liblinalg


