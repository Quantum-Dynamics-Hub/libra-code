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

#include "Integral_Approx1.h"
#include "../math_specialfunctions/libspecialfunctions.h"


/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;

namespace libmolint{


double A_plus_1(int n,double alpha){

  double res,d;
  double a_n,n_fact;

  d      = 1.0;
  n_fact = 1.0;
  a_n    = 1.0;

  res = d;
	
  for(double v=1.0;v<=n;v++){

      d      = d*(alpha/v);
      n_fact = n_fact*v;
      a_n    = a_n*alpha;

      res += d;
  }
      a_n    = a_n*alpha;

      res = res*exp(-alpha)*(n_fact/a_n);

  return res;
}


double A_minus_1(int n,double alpha){

  double res,d;
  double a_n,n_fact;

  d      = 1.0;
  n_fact = 1.0;
  a_n    = 1.0;

  res = d;
	
  for(float v=1;v<=n;v++){

      d      = d*(-alpha/v);
      n_fact = n_fact*v;
      a_n    = a_n*alpha;

      res += d;
  }
      a_n    = a_n*alpha;

      res = res*exp(alpha)*(n_fact/a_n);

  return res;
}

double B(int n,double alpha){

  double res;
  res = A_minus_1(n,alpha)-A_plus_1(n,alpha);

  return res;
}

double D(int m,int n,int p){

  double res=0.0;
  int k_min,k_max;

  k_min = ((p-m)>0)?(p-m):0;
  k_max = (p<n)?p:n;

  for(int k=k_min;k<=k_max;k++){

      res += pow(-1.0,k)*BINOM((p-k),m)*BINOM(k,n);
  }

  return res;

}

double K2aa(int m,int n,double alpha,double betha,double R){
//  Coulomb integral on Slater functions
//  [alpha == betha]

  double K1,K2;
  double r,t,f,sum;
  double rel;

  r = alpha*R;
  r = pow(r,2*m);
  t = pow(2.0,2*m);
  f = FACTORIAL(2*m);

  if(R!=0.0){

      K1  = (1.0/R);
      sum = (2.0*alpha*R-2.0*m)*(A_plus_1((2*m-1),(2.0*alpha*R)) - exp(-2.0*alpha*R));
      K1 += sum*(t*r/(f*R));

      rel=1.0;
      K2 = 0.0;

      for(float v=0;v<=(2*n-1);v++){

	   sum =0.0;
	   for(int p=0;(2.0*p<=(2.0*m+v-1) );p++){

	      sum += (1.0/(2.0*p+1.0))*D((2*m-1),v,2*p)*A_plus_1((2*m-1+v-2*p),2.0*R*alpha);
	   }

	   if(v==0){rel = 1.0;}
	   else{rel = rel*(betha*R/v);} 

	   K2 += (2*n-v)*rel*sum;

      }// for v

      double norm;
      norm = (alpha/((double)n));
      norm = norm*(r/f);

      K2 = K1-K2*norm;

  }// if R!=0
	       
  else if(R==0.0){

      K1 = (alpha/(double)m);
      K2 = 0.0;

      for(float v=0;v<=(2*n-1);v++){
							
	   if(v==0){rel = FACTORIAL(2*m-1);}
	   else{rel = rel*(0.5*(2*m-1+v)/v);} 

	   K2 += (2*n-v)*rel;
      }// for v

      K2 = K1 - K2*(alpha/((double)n))*(1.0/(t*f));

  }// if R==0.0
	
  return K2;

}

double K2ab(int m,int n,double alpha,double betha,double R){
//  Coulomb integral on Slater functions
//  [alpha != betha]

  double K1,K2;
  double r,t,f,sum;

  r = alpha*R;
  r = pow(r,2*m);
  t = pow(2.0,2*m);
  f = FACTORIAL(2*m);

  K1  = (1.0/R);
  sum = (2.0*alpha*R-2.0*m)*(A_plus_1((2*m-1),(2.0*alpha*R)) - exp(-2.0*alpha*R));
  K1 += sum*(t*r/(f*R));

  double rel=1.0;
  K2 = 0.0;

  for(float v=0;v<=2*n;v++){

      sum =0.0;
      for(int p=0;p<=(2*m-1+v);p++){

	   sum += D((2*m-1),v,p)*B(p,R*(alpha-betha))*A_plus_1((2*m-1+v-p),R*(alpha+betha));
      }

      if(v==0){rel = 1.0;}
      else{rel = rel*(betha*R/v);} 

      K2 += (2*n-v)*rel*sum;

  }// for v

  double norm;
  norm = 0.5*(alpha/((double)n));
  norm = norm*(r/f);

  K2 = K1-K2*norm;

  return K2;
}


double Jab(int m,int n,double alpha,double betha,double R){
// This is the same as Kab
//  Coulomb integral on Slater functions
//  [alpha != betha]

  double K1,K2;
  double r,t,f,sum;
  double rel;

  if(alpha!=betha){//----------------------------------------------------------

      r = alpha*R;
      r = pow(r,2*m);
      t = pow(2.0,2*m);
      f = FACTORIAL(2*m);

      K1  = (1.0/R);
      sum = (2.0*alpha*R-2.0*m)*(A_plus_1((2*m-1),(2.0*alpha*R)) - exp(-2.0*alpha*R));
      K1 += sum*(t*r/(f*R));

      rel=1.0;
      K2 = 0.0;

      for(float v=0;v<=2*n;v++){ //< because if v==2n => 2n-v = 0, zero addition

	   sum =0.0;
	   for(int p=0;p<=(2*m-1+v);p++){

	      sum += D((2*m-1),v,p)*B(p,R*(alpha-betha))*A_plus_1((2*m-1+v-p),R*(alpha+betha));
	   }

	   if(v==0){rel = 1.0;}
 	   else{rel = rel*(betha*R/v);} 

	   K2 += (2*n-v)*rel*sum;

      }// for v

      double norm;
      norm = 0.5*(alpha/((double)n));
      norm = norm*(r/f);

      K2 = K1-K2*norm;

  } //if(alpha!=betha)
  else{//---------------------------------------------------------------------

      betha = alpha = 0.5*(alpha+betha);

      r = alpha*R;
      r = pow(r,2*m);
      t = pow(2.0,2*m);
      f = FACTORIAL(2*m);

      K1  = (1.0/R);
      sum = (2.0*alpha*R-2.0*m)*(A_plus_1((2*m-1),(2.0*alpha*R)) - exp(-2.0*alpha*R));
      K1 += sum*(t*r/(f*R));
	
      rel=1.0;
      K2 = 0.0;

      for(float v=0;v<=(2*n-1);v++){

	   sum =0.0;
	   for(float p=0;(p<=(m+0.5*(v-1)) );p++){

	      sum += (1.0/(2.0*p+1.0))*D((2*m-1),v,2*p)*A_plus_1((2*m-1+v-2*p),2.0*R*alpha);
	   }

	   if(v==0){rel = 1.0;}
	   else{rel = rel*(alpha*R/v);} 

	   K2 += (2*n-v)*rel*sum;

      }// for v

      double norm;
      norm = (alpha/((double)n));
      norm = norm*(r/f);

      K2 = K1-K2*norm;

  }// if alpha==betha

  double convert=1.0;

// Keep result in atomic units (Hartee)
//  convert = 27.211383; // Hartree -> eV
 

  return K2*convert;


}

double Jab(int m,int n,double alpha,double betha){

// Idempotential for some particular atom

  double K1,K2;
  double r,t,f,sum;
  double rel;

  if(alpha==betha){

      t = pow(2.0,2*m);
      f = FACTORIAL(2*m);
      K1 = (alpha/((double)m));

      rel=1.0;
      K2 = 0.0;

      for(float v=0;v<=(2*n-1);v++){

	   if(v==0){rel = (f/(2.0*m));}
	   else{rel = rel*((2.0*m+v-1)/(2.0*v));} 

 	   K2 += (2*n-v)*rel;

      }

      K2 = K2*alpha/(t*n*f);
      K2 = K1 - K2;

  }// if alpha==betha

  double convert=1.0;

// Keep result in atomic units (Hartree)
//  convert = 27.211383; // Hartree -> eV

  return K2*convert;

}


double dJab_dqa(int m,int n,double alpha,double betha,double R){

// This derivative is valid only for H atom where alpha = alpha0 + q;

  double K2,sum,sum1;
  double rel=1.0;
  K2 = 0.0;

  for(float i=0;i<=2*n;i++){

      sum =0.0;
      for(int j=0;j<=(1+i);j++){

	   sum1 = (-i*alpha*alpha+(i-2.0*j+1.0)*alpha*betha-3.0*betha*betha)/(alpha*alpha-betha*betha);
	   sum1 = sum1*B(j,R*(alpha-betha))*A_plus_1(1+i-j,R*(alpha+betha));

	   sum1 += (alpha/(alpha-betha))*A_plus_1(1+i-j,R*(alpha+betha))*exp(-R*(alpha-betha));

	   sum1 += pow(-1.0,j)*exp(R*(alpha-betha));

	   sum1 += (alpha/(alpha+betha))*B(j,R*(alpha-betha))*exp(-R*(alpha+betha));
	
	   sum += D(1,j,j)*sum1;
	
      }// for j

      if(i==0){rel = 1.0;}
      else{rel = rel*(betha*R/i);} 

      K2 += (2*n-i)*rel*sum;

  }

  K2 = K2*(alpha*alpha*R*R)/(4.0*(double)n);
  K2 = (2.0*alpha*R)*(2.0*alpha*R)*A_plus_1(1,(2.0*alpha*R)) - K2;


  return K2;
}



double Coulomb_Integral(double R,int n_i, double Jii, double ksi_i, std::string type_i, double q_i,
                                 int n_j, double Jjj, double ksi_j, std::string type_j, double q_j,
                                 double epsilon, int mode){
/**
  Literature references:

[1] - A. K. Rappe, W. A. Goddard III "Charge Equilibration for Molecular Dynamics
      Simulations" J. Phys. Chem. 1991 V. 95, P. 3358-3363

[2] - A. Oda, S. Hirono "Geometry-dependent charge calculations using charge
      equilibration method with empirical two-center Coulombic terms" J. Mol. Struct.
      (Theochem) 2003 V. 634, P. 159-170
     Here they found that Ohno-Klopman and DasGupta-Huzinaga equationa are most
     efficient and that Ohno-Klopman rule with parameters of Bakowies and Thiel is
     the best combination for calculation of dipole moments

[3] - J. N. Louwen, E. T. C. Vogt "Semi-empirical atomic charges for use in
      computational chemistry of molecular sieves" J. Mol. Catal. A. 1998, V. 134,
      P. 63 - 77

*/

  // Jij is in Hartrees (a.u. of energy)
  // R should be in Atomic units
  // epsilon - is a dielectric constant

  double Jij=0.0;

  switch(mode){

  case 0: {//--------- Calculate (charge-dependent idempotential)----------------
           // This method obeys to those formula presented in original QEq approach
           // by Rappe and Goddard [1]
           Jij = Jab(n_i,n_j,ksi_i,ksi_j,R);
           if( (type_i=="H") && (type_j=="H") ){
           Jij = Jij * (1.0 + q_i/ksi_i);  // Make only H idempotentials charge dependent,[ what is q_i != q_j ?? ]
           }
          }break;

  case 1: {//-------- Calculate by simple formula ---------------------------------
           // This is a simplest way to calculate Jab integral proposed in original
           // QEq paper [1]
           Jij = 1.0/(epsilon*R);
          }break;

  case 2: {//-------- Nishimoto-Mataga equation (from [2]) -----------------------
           Jij = 1.0/(R + 2.0/(Jii + Jjj));
          }break;

  case 3: {//------- Nishimoto-Mataga-Weiss equation (from [2]) ------------------
           double fr = 1.2;
           Jij = fr/(R +2.0*fr/(Jii + Jjj));
          }break;

  case 4: {//-------------- Ohno equation (from [2]) ----------------------------
           double a = 2.0/(Jii + Jjj);
           Jij = 1.0/sqrt(R*R + a*a);
          }break;

  case 5: {//------------- Ohno-Klopman equation (from [2]) ---------------------
           double a = 0.5*((1.0/Jii) + (1.0/Jjj));
                      Jij = 1.0/sqrt(R*R + a*a);
          }break;

  case 6: {//----------- DasGupta - Huzinaga equation (from [2]) ------------------
           double ki,kj; ki = kj = 0.4; // Klondike parameter
           double a = 0.5*(Jii*exp(ki*R) + Jjj*exp(kj*R));
           Jij = 1.0/(R + (1.0/a));
          }break;
  case 7: {//-------------------- Louwen-Vogt [3] --------------------------------
           double a = 2.0/(Jii + Jjj);
           Jij = 1.0/std::pow((R*R*R + a*a*a),(1.0/3.0));
          }break;

  default: {
  std::cout<<"Warning: To calculate two-center Coulomb integral there are options:"<<std::endl;
  std::cout<<"0, 1, ... , 7. Exiting"<<std::endl;
  exit(0);
           }
  }

  return Jij;

}// Coulomb_Integral


}// namespace libmolint
}// namespace liblibra

