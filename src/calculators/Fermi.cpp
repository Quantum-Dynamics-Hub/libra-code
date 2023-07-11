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
/**
  \file Fermi.cpp
  \brief The file implements functions for Fermi energy/population calculations
    
*/

#include "Fermi.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libcalculators namespace
namespace libcalculators{



double fermi_population(double e,double ef,double degen, double kT){
/** 
  \brief Fermi populations.

  Compute Fermi-based occupation of energy level with energy e
  ef - is the assumed Fermi energy
  double kT = 0.00095;  kT in a.u. for T = 300 K

  \param[in] e  Energy level
  \param[in] ef Fermi energy
  \param[in] degen Degeneracy of the energy levels
  \param[in] kT  Broadening factor for Fermi distribution

*/


  double argg = ((e-ef)/kT);
  double pop;

  if(argg>50.0){ pop = 0.0;  }
  else if(argg<-50.0){   pop = degen;   }
  else{ pop = degen/(1.0 + std::exp(argg));  }  

  return pop;

}// fermi_population


double fermi_integral(std::vector<double>& bnds, double ef, double degen, double kT){
/**
  \brief Auxiliary function to compute Fermi integral

  \param[in] bnds 
  Compute integral(actually the sum):
    sum [degen / (1 + exp(e-ef))]   
     i
  where N - is a number of electrons (valence)
  double kT = 0.00095; // kT in a.u. for T = 300 K

  \param[in] bnds Input band (energy level index and the energy)
  \param[in] ef Fermi energy
  \param[in] degen Degeneracy of the energy levels
  \param[in] kT  Broadening factor for Fermi distribution

*/


  int Norb = bnds.size();

  double sum = 0.0;
  for(int i=0;i<Norb;i++){

    double argg = ((bnds[i] - ef)/kT);

    if(argg>50.0){   }
    else if(argg<-50.0){   sum += degen;   }
    else{  sum += degen/(1.0 + std::exp(argg));   }  

  }// for i

  return sum;

}// double fermi_integral


double fermi_integral(boost::python::list bnds,double ef,double degen, double kT){
/**
  \brief Auxiliary function to compute Fermi integral (Python-friendly)

  \param[in] bnds 
  Compute integral(actually the sum):
    sum [degen / (1 + exp(e-ef))]   
     i
  where N - is a number of electrons (valence)
  double kT = 0.00095; // kT in a.u. for T = 300 K

  \param[in] bnds Input band (energy level index and the energy)
  \param[in] ef Fermi energy
  \param[in] degen Degeneracy of the energy levels
  \param[in] kT  Broadening factor for Fermi distribution

*/


  int sz = boost::python::len(bnds);
  vector<double> int_bnds(sz,0.0);
  for(int i=0;i<sz;i++){ 
    int_bnds[i] = boost::python::extract<double>(bnds[i]);
  }

  return fermi_integral(int_bnds, ef, degen, kT);

}


double fermi_energy(std::vector<double>& bnds,double Nel,double degen, double kT, double etol){
/**
  \brief Fermi energy solver

  Computes Fermi energy by solving equation
  fermi_integral( ...  ef ... )  = Nel
  Using bisection method

  \param[in] bnds Input band (energy level index and the energy)
  \param[in] Nel The number of electrons to distribute on energy levels
  \param[in] degen Degeneracy of the energy levels
  \param[in] kT  Broadening factor for Fermi distribution
  \param[in] etol Tolerance level (stop when 0.5*|e_f(old) - e_f(new)|<tol)

*/


  double ef_l,ef_m,ef_r,i_l,i_m,i_r;
  double err = 2.0*etol;
  int Norb = bnds.size();

  ef_l = bnds[0] - 10.0;
  ef_r = bnds[Norb-1] + 10.0;

  i_l = fermi_integral(bnds,ef_l,degen,kT) - Nel;
  i_r = fermi_integral(bnds,ef_r,degen,kT) - Nel;

  do{
    ef_m = 0.5*(ef_l + ef_r);
    i_m = fermi_integral(bnds,ef_m,degen,kT) - Nel;

    if(0){ // For debug
      cout<<"ef_l= "<<ef_l<<" ef_m= "<<ef_m<<" ef_r= "<<ef_r<<endl;
      cout<<"i_l= "<<i_l<<" i_m= "<<i_m<<" i_r= "<<i_r<<endl;
    }

    int var;
    if(i_m*i_r<=0.0 && i_m*i_l>=0.0){ var = 1; }
    else if(i_m*i_r>=0.0 && i_m*i_l<=0.0){ var = 2; }
    else{ cout<<"Error in fermi_energy\n"; exit(0); }

    switch(var){
      case 1: {i_l = i_m; ef_l = ef_m; break; }
      case 2: {i_r = i_m; ef_r = ef_m; break; }
      default: break;
    }

    err = 0.5*(ef_r-ef_l); 
//    err = 0.5*(i_r - i_l);
    
  }while(err>etol);
 
  return 0.5*(ef_l+ef_r);

}// double fermi_energy


double fermi_energy(boost::python::list bnds,double Nel,double degen, double kT, double etol){
/**
  \brief Fermi energy solver (Python-friendly version)

  Computes Fermi energy by solving equation
  fermi_integral( ...  ef ... )  = Nel
  Using bisection method

  \param[in] bnds Input band (energy level index and the energy)
  \param[in] Nel The number of electrons to distribute on energy levels
  \param[in] degen Degeneracy of the energy levels
  \param[in] kT  Broadening factor for Fermi distribution
  \param[in] etol Tolerance level (stop when 0.5*|e_f(old) - e_f(new)|<tol)

*/


  int sz = boost::python::len(bnds);
  vector<double> int_bnds(sz,0.0);
  for(int i=0;i<sz;i++){ 
    int_bnds[i] = boost::python::extract<double>(bnds[i]);
  }

  return fermi_energy(int_bnds, Nel, degen, kT, etol);

}





///========= The following functions are needed for 
/// FOE: Fermi Operator Expansion - a linear-scaling electronic 
/// structure methodology

double p_up(double e, double e_up, double de){
  double res = 0.0;
  if(e<e_up){ res = 0.0; }
  else{ 
    double argg = (e-e_up)/de; 
    if(argg<=100){ res = exp(argg); }
    else{ res = exp(100.0); }
  }
  return res;
}

double p_dn(double e, double e_dn, double de){
  double res = 0.0;
  if(e>e_dn){ res = 0.0; }
  else{ 
    double argg = -(e-e_dn)/de; 
    if(argg<=100){ res = exp(argg); }
    else{ res = exp(100.0); }
  }
  return res;
}

double p_ef(double e, double ef, double de){
  double res = 0.0;
  double argg = (e - ef)/de;
  if(argg<-100){ res = 1.0; }
  else if(argg>100.0){ res = 0.0; }
  else{ res = 1.0/(1.0+exp(argg)); }
  return res;
}


void Chebyshev_coeff(vector<double>& C, double (*f)(double x, double y, double z), double ef, double de, int N){
/** According to Numerical Recipes

// Computes expansion coefficients of scalar Fermi function
// C  - must be initialized and contains N elements: C[0]... C[N-1]
// double (*f)(double x, double y, double z) -is a scalar function which we are approximating
// ef - trial parameter (Fermi energy, upper or lower eigenvalue)
// de - energy width parameter (smaller it is, the higher the accuracy, but is slower calculation)
// np - degree of polynomial expansion

*/

  float np = N;

  for(int j=0;j<N;j++){  C[j] = 0.0;   }// for k

  for(int k=0;k<N;k++){                 // over all points

    double x_k = cos((k+0.5)*M_PI/np);   // k-th zero of polynomial T_N

    // Function (Fermi) at the x_k point
    double f_k = f(x_k, ef, de) * (2.0/np);

    double t_curr = 0.0;
    double t_prev1, t_prev2;
      
    for(int j=0;j<N;j++){

      if(j==0){  t_curr = 1.0; } // T_0
      else if(j==1){  t_prev1 = t_curr; t_curr = x_k; } // T_1
      else{                                             // T_j
        t_prev2 = t_prev1; t_prev1 = t_curr;
        t_curr = 2.0*x_k * t_prev1 - t_prev2; 
      }

      C[j] += f_k * t_curr;

    }// for j all coefficients
  }// for k all points


//  cout<<"Chebyshev coefficients\n";
//  for(int k=0;k<=np;k++){  cout<<"k= "<<k<<"  C["<<k<<"]= "<<C[k]<<endl;   }// for k

}

double Chebyshev_fit(MATRIX& H, MATRIX& P, double (*f)(double _x, double _y, double _z), double ef, double de, int np){
/**
// Compute density matrix using Chebyshev polynomials
// Return trace of the density matrix, for it should be equal to the number of electrons
*/

  int N = H.n_cols; // assume square matrix

  // Compute coefficients
  vector<double> c_k(np+1,0.0);   Chebyshev_coeff(c_k, f, ef,de,np);


  // Recursive definition of Chebyshev matrices
  P = 0.0;
  MATRIX T_k(N,N);     
  MATRIX T_prev1(N,N); 
  MATRIX T_prev2(N,N); 
  
  for(int k=0;k<=np;k++){

    if(k==0){  T_k.Init_Unit_Matrix(1.0);}
    else if(k==1){  T_prev1 = T_k; T_k = H; }
    else{
      T_prev2 = T_prev1;  T_prev1 = T_k;
      T_k = 2.0*H*T_prev1 - T_prev2; 
    }


    P += c_k[k] * T_k; 
    if(k==0){  P -= 0.5*c_k[0]*T_k; }    // T_k is simply unity matrix in this case

  }// for k

  return P.tr();    

}





}// namespace libcalculators
}// liblibra


