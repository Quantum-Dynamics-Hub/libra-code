#include "Overlaps.h"
#include "Integral_Electron_Repulsion.h"
#include "Aux_Functions.h"


namespace libqchem{
namespace libmolint{


/*
double ELECTRON_REPULSION_INTEGRAL(PrimitiveG& GA,PrimitiveG& GB,PrimitiveG& GC,PrimitiveG& GD,VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD,
       vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv){
// This is equivalent to chemists' notation: (ab|cd) 
//                 or
//                    to physicists notation: <ac|bd>

//    cout<<"In ERI calculations:\n";

//   cout<<"A="<<*GA.G_center<<" B="<<*GB.G_center<<" C= "<<*GC.G_center<<" D= "<<*GD.G_center<<endl;


    VECTOR R_AB,R_CD,P,Q,PA,PB,QC,QD,PQ;
    R_AB = *GA.G_center - *GB.G_center;
    R_CD = *GC.G_center - *GD.G_center;

    int maxI = GA.x_exp + GB.x_exp + GC.x_exp + GD.x_exp;
    int maxJ = GA.y_exp + GB.y_exp + GC.y_exp + GD.y_exp;
    int maxK = GA.z_exp + GB.z_exp + GC.z_exp + GD.z_exp;

    double gamma1,gamma2;
    gamma1 = GA.G_alpha + GB.G_alpha;
    gamma2 = GC.G_alpha + GD.G_alpha;
//    cout<<"gamma1 = "<<gamma1<<" gamma2 = "<<gamma2<<endl;

    P = (GA.G_alpha*(*GA.G_center) + GB.G_alpha*(*GB.G_center))/gamma1;
    Q = (GC.G_alpha*(*GC.G_center) + GD.G_alpha*(*GD.G_center))/gamma2;

    PA = P - (*GA.G_center);
    PB = P - (*GB.G_center);
    QC = Q - (*GC.G_center);
    QD = Q - (*GD.G_center);

    PQ = Q - P;


    // Jacobian
    double dPA_dA  = (GA.G_alpha/gamma1 - 1.0);
    double dPA_dB  = (GB.G_alpha/gamma1);
    double dPA_dC  = 0.0;
    double dPA_dD  = 0.0;

    double dPB_dA  = (GA.G_alpha/gamma1);
    double dPB_dB  = (GB.G_alpha/gamma1 - 1.0);
    double dPB_dC  = 0.0;
    double dPB_dD  = 0.0;

    double dQC_dA  = 0.0;
    double dQC_dB  = 0.0;
    double dQC_dC  = (GC.G_alpha/gamma2 - 1.0);
    double dQC_dD  = (GD.G_alpha/gamma2);

    double dQD_dA  = 0.0;
    double dQD_dB  = 0.0;
    double dQD_dC  = (GC.G_alpha/gamma2);
    double dQD_dD  = (GD.G_alpha/gamma2 - 1.0);

    double dPQ_dA  = -(GA.G_alpha/gamma1);
    double dPQ_dB  = -(GB.G_alpha/gamma1);
    double dPQ_dC  = (GC.G_alpha/gamma2);
    double dPQ_dD  = (GD.G_alpha/gamma2);






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




    Aux_Function6(GA.x_exp,GB.x_exp,GC.x_exp,GD.x_exp, PA.x,PB.x,QC.x,QD.x,PQ.x, gamma1,gamma2,  GI, dGI_dPA, dGI_dPB, dGI_dQC, dGI_dQD, dGI_dPQ, aux[18],aux[19],aux[20], aux[21],aux[22],aux[23], aux[24],aux[25],aux[26], n_aux);
    Aux_Function6(GA.y_exp,GB.y_exp,GC.y_exp,GD.y_exp, PA.y,PB.y,QC.y,QD.y,PQ.y, gamma1,gamma2,  GJ, dGJ_dPA, dGJ_dPB, dGJ_dQC, dGJ_dQD, dGJ_dPQ, aux[18],aux[19],aux[20], aux[21],aux[22],aux[23], aux[24],aux[25],aux[26], n_aux);
    Aux_Function6(GA.z_exp,GB.z_exp,GC.z_exp,GD.z_exp, PA.z,PB.z,QC.z,QD.z,PQ.z, gamma1,gamma2,  GK, dGK_dPA, dGK_dPB, dGK_dQC, dGK_dQD, dGK_dPQ, aux[18],aux[19],aux[20], aux[21],aux[22],aux[23], aux[24],aux[25],aux[26], n_aux);





//    cout<<"maxI = "<<maxI<<" maxJ = "<<maxJ<<" maxK = "<<maxK<<endl;



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


        }// for K
      }// for J
    }// for I



    // Precompute inclomplete Gamma functions:
    double* F_nu;  F_nu = aux[28];
    double d4 = ((1.0/gamma1) + (1.0/gamma2));
    for(int nu=0;nu<=(maxI+maxJ+maxK+1); nu++){
        F_nu[nu] = Fn(nu,PQ.length2()/d4) // Aux_Function2(nu,PQ.length2()/d4);
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

//      cout<<"nu="<<nu<<" C[nu]= "<<C_nu[nu]<<" F[nu]= "<<F_nu[nu]<<endl;

      dERI_dA += ( dCdA_nu[nu]*F_nu[nu] + C_nu[nu]*(-(2.0/d4)*F_nu[nu+1]*PQ*dPQ_dA) );
      dERI_dB += ( dCdB_nu[nu]*F_nu[nu] + C_nu[nu]*(-(2.0/d4)*F_nu[nu+1]*PQ*dPQ_dB) );
      dERI_dC += ( dCdC_nu[nu]*F_nu[nu] + C_nu[nu]*(-(2.0/d4)*F_nu[nu+1]*PQ*dPQ_dC) );
      dERI_dD += ( dCdD_nu[nu]*F_nu[nu] + C_nu[nu]*(-(2.0/d4)*F_nu[nu+1]*PQ*dPQ_dD) );

    }// for nu



    double pref0 = (2.0*M_PI*M_PI/(gamma1*gamma2))*sqrt(M_PI/(gamma1+gamma2));
    double pref_AB = exp(-GA.G_alpha*GB.G_alpha*R_AB.length2()/gamma1);
    double pref_CD = exp(-GC.G_alpha*GD.G_alpha*R_CD.length2()/gamma2);


    ERI = pref0 * pref_AB * pref_CD * ERI;

    DA = pref0 * pref_AB * pref_CD * ( dERI_dA + ERI*(-2.0*GA.G_alpha*GB.G_alpha/gamma1)*R_AB );   
    DB = pref0 * pref_AB * pref_CD * ( dERI_dB + ERI*( 2.0*GA.G_alpha*GB.G_alpha/gamma1)*R_AB );
    DC = pref0 * pref_AB * pref_CD * ( dERI_dC + ERI*(-2.0*GC.G_alpha*GD.G_alpha/gamma2)*R_CD );
    DD = pref0 * pref_AB * pref_CD * ( dERI_dD + ERI*( 2.0*GC.G_alpha*GD.G_alpha/gamma2)*R_CD );

//    cout<<"ERI = "<<ERI<<endl;
//    cout<<"pref0= "<<pref0<<endl;
//    cout<<"pref_AB= "<<pref_AB<<endl;
//    cout<<"pref_CD= "<<pref_CD<<endl;
//    cout<<"dERI_dA= "<<dERI_dA<<endl;
//    cout<<"dERI_dB= "<<dERI_dB<<endl;
//    cout<<"dERI_dC= "<<dERI_dC<<endl;
//    cout<<"dERI_dD= "<<dERI_dD<<endl;


    return ERI;

}// eri


*/

}// namespace libqchem
}// namespace libmolint


