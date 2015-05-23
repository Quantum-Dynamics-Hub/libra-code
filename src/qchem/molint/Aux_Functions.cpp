#include "Aux_Functions.h"


namespace libqchem{
namespace libmolint{


//===========================================================================================

double Aux_Integral2(int n1,int n2,double x1,double x2,double alp,double& dI_dx1,double& dI_dx2,vector<double*>& aux,int n_aux){
// This version will not call outside function - all is here
// dI_dx1, dI_dx2 - are the derivatives of the integral w.r.t. variables x1 and x2, respectively.

  double* f;  f = aux[0];
  double* dfdx1; dfdx1 = aux[1];
  double* dfdx2; dfdx2 = aux[2];
  double* PW1; PW1 = aux[3];
  double* dPW1; dPW1 = aux[4];
  double* PW2; PW2 = aux[5];
  double* dPW2; dPW2 = aux[6];

  for(int k=0;k<=(n1+n2);k++){ f[k] = 0.0; dfdx1[k] = 0.0; dfdx2[k] = 0.0; }

  // Pre-compute PW1 and PW2 arrays independently:
  double pw,dpw,bin;
  
  for(int i=n1;i>=0;i--){

    int ni1 = n1 - i;

    if(ni1==0){ pw = 1.0; dpw = 0.0;    }
    else{  dpw = pw + x1*dpw;  pw *= x1;  } // order matters here!!!
   
    bin = BINOM(i,n1);

    PW1[i] = bin * pw;  
    dPW1[i] = bin * dpw;

  }

  for(int i=n2;i>=0;i--){

    int ni2 = n2 - i;

    if(ni2==0){ pw = 1.0; dpw = 0.0;    }
    else{  dpw = pw + x2*dpw;  pw *= x2;  } // order matters here!!!
   
    bin = BINOM(i,n2);

    PW2[i] = bin * pw;  
    dPW2[i] = bin * dpw;

  }

  // Now do the double loop
  for(int i=0;i<=n1;i++){
    for(int j=0;j<=n2;j++){

      f[i+j] += PW1[i]*PW2[j];
      dfdx1[i+j] += dPW1[i]*PW2[j];
      dfdx2[i+j] += PW1[i]*dPW2[j];

    }// for j
  }// for i

  

  double I = 0.0;
  dI_dx1 = 0.0;
  dI_dx2 = 0.0;


  double rat = (0.5/alp);
  double w = sqrt(M_PI/alp);
  for(int k=0; k<=(n1+n2); k+=2){   // this implies that k%2 = 0 - always

    // The following few lines inline this one:
    if(k>0){   w *= (k-1) * rat;  }

    I += f[k] * w;
    dI_dx1 += dfdx1[k] * w;
    dI_dx2 += dfdx2[k] * w;
  }

  return I;
}




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


void Aux_Function5(int l1,int l2,double a,double b,double gamma,
                   double* H,double* dHda, double* dHdb,
                   double* f, double* dfda, double* dfdb, int n_aux){

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



void Aux_Function6(int l1,int l2,int l3,int l4,double PA,double PB,double QC,double QD,double p,double gamma1,double gamma2,
 double* G, double* dGdA, double* dGdB, double* dGdC, double* dGdD, double* dGdp, 
 double* HL,double* dHLda, double* dHLdb,
 double* HM,double* dHMda, double* dHMdb,
 double* f, double* dfda, double* dfdb,
 int n_aux ){

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

  Aux_Function5(l1,l2,PA,PB,gamma1,HL,dHLda,dHLdb,f,dfda,dfdb,n_aux); // HL of size l1+l2
  Aux_Function5(l3,l4,QC,QD,gamma2,HM,dHMda,dHMdb,f,dfda,dfdb,n_aux); // HM of size l3+l4


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

}




}// namespace libmolint
}// namespace libqchem


//=============================================================================================

