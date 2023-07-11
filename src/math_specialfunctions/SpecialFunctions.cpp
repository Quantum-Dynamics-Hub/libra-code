/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "SpecialFunctions.h"
#include "../math_meigen/mEigen.h"


//================== Functions ==========================
// Taylor expansions for some formula are from:
// http://web.mit.edu/kenta/www/three/taylor.html

/// liblibra namespace
namespace liblibra{


/// libspecialfunctions namespace
namespace libspecialfunctions{


double FAST_POW(double x,int n){

    double res,res1;    
    if(n==0)      { res = 1.0; }
    else if(n==1) { res = x; }                                    // 0 mult
    else if(n==2) { res = x*x;}                                   // 1 mult
    else if(n==3) { res = x*x*x; }                                // 2 mult
    else if(n==4) { res = x*x; res = res*res; }                   // 2 mult
    else if(n==5) { res = x*x; res = res*res*x; }                 // 3 mult
    else if(n==6) { res = x*x*x; res = res*res; }                 // 3 mult
    else if(n==7) { res = x*x*x; res = res*res*x; }               // 4 mult
    else if(n==8) { res = x*x; res = res*res; res = res*res;}     // 3 mult
    else if(n==9) { res = x*x*x; res = res*res*res; }             // 4 mult
    else if(n==10){ res1= res = x*x; 
                    res = res*res;
                    res = res*res;
                    res = res*res1;
                  }                                               // 4 mult
    else if(n==11){ res1= res = x*x;
                    res = res*res;
                    res = res*res;
                    res = res*res1*x;
                  }                                               // 5 mult
    else if(n==12){ res = x*x*x;
                    res = res*res;
                    res = res*res;
                  }                                               // 4 mult
    else{   res = pow(x,n);   }

    return res;

}


double sinh_(double x){
/**  This is sinh(x)/x function  */

    double tol = 1e-12;
    double res,x2,x4,x6,x8;
    double c1 = (1.0/6.0);
    double c2 = (1.0/120.0);
    double c3 = (1.0/5040.0);
    double c4 = (1.0/362880);

    if(fabs(x)>tol){  res = (sinh(x)/x); }
    else{ // Use Maclaurin series
          x2 = x*x;
          x4 = x2*x2;
          x6 = x2*x4;
          x8 = x2*x6;

          res = 1.0 + c1*x2 + c2*x4 + c3*x6 + c4*x8;
    }

    return res;

}

double sin_(double x){
/**  This is sin(x)/x function */

    double tol = 1e-12;
    double res,x2,x4,x6,x8,x10;
    double c1 = (1.0/6.0);
    double c2 = (1.0/120.0);
    double c3 = (1.0/5040.0);
    double c4 = (1.0/362880.0);
    double c5 = (1.0/39916800.0);
    double c6 = (1.0/6227020800.0);

    if(fabs(x)>tol){  res = (sin(x)/x); }
    else{ // Use Maclaurin series
          x2 = x*x;
          x4 = x2*x2;
          x6 = x2*x4;
          x8 = x2*x6;
          x10 = x2*x6;

          res = 1.0 - c1*x2 + c2*x4 - c3*x6 + c4*x8 - c5*x10;
    }

    return res;

}

double ERF(double x){
/** Error function approximation:

Source:
http://homepages.physik.uni-muenchen.de/~Winitzki/erf-approx.pdf

Article:
S. Winitzki, Uniform approximations for transcendental functions, in Proc. ICCSA-
2003, LNCS, 2667/2003, p. 962

/Theory/erf-approx.pdf

*/

   int signx = 1;
   double x2,ax2;

   if(x<=0) { x = -x; signx = -1; }

   x2 = x * x;
   ax2 = 0.147*x2;

   double res = sqrt(1.0 - exp(-x2*(1.273239545 + ax2)/(1.0 + ax2)) );

   res = res*signx;

   return res;

}

double ERFC(double x){
/** Error function approximation:

Source:
http://homepages.physik.uni-muenchen.de/~Winitzki/erf-approx.pdf

Article:
S. Winitzki, Uniform approximations for transcendental functions, in Proc. ICCSA-
2003, LNCS, 2667/2003, p. 962

/Theory/erf-approx.pdf

*/

   int signx = 1;
   double x2,ax2;

   if(x<=0) { x = -x; signx = -1; }

   x2 = x * x;
   ax2 = 0.147*x2;

   double res = sqrt(1.0 - exp(-x2*(1.273239545 + ax2)/(1.0 + ax2)) );

   res = 1.0 - res*signx;

   return res;

}// ERFC


double gamma_lower(double s,double x){
/** This computes lower incomplete gamma function divided by power:
 gamma_tilda(s,x) = gamma(s,x)/ x^s, 
 gamma(s,x) is defined in the link below:
 http://en.wikipedia.org/wiki/Incomplete_gamma_function
 We use series expansion, as in the link
*/

  double res = 0.0;
  double t1,term;

  if(x<0){ x = 0.0; }

  term = exp(-x)/s;
  res = term;

  int iter = 1;
  while((term>0.0) && (iter<1000)){
    term = term * (x/(s + iter));
    res += term;
    iter++;
  }

    
  return res;
}


double Fn(int n,double t){
/** This computes the incomplete gamma function given by an integral
            1
  F_n(t) =  | u^2n * exp(-t*u^2) du 
            0
 
  which is equal to: gamma(n+1/2,t)/ (2* t^{n+1/2}) = 0.5*gamma_lower(n+1/2,t)
  where gamma(s,x) - is lower incomplete gamma-function:
  http://en.wikipedia.org/wiki/Incomplete_gamma_function
*/

  double res = 0.5*gamma_lower((n+0.5),t);

  return res; 

}


double gaussian_int(int n, double alp){
/****************************************************************************
 This function computes the elementary integral

   +inf
  int   {   x^n * exp(-alp*x^2) dx }   = 0 , if n is odd,  !=0, if n is even:
   -inf

 if n is even, we use:  n = 2k

   +inf
  int   {   x^2k * exp(-alp*x^2) dx }   = ( (2k-1)!! /( (2a)^k ) ) * sqrt(pi/a)
   -inf

  http://en.wikipedia.org/wiki/Gaussian_integral
****************************************************************************/

  double res = 0.0;

  if(n%2==0){  // even    
    if(n==0){  res = sqrt(M_PI/alp); }
    else{ 
      res = (DFACTORIAL(n-1)/pow(2.0*alp, n/2))*sqrt(M_PI/alp);
    }
  }
  else{        // odd
    res = 0.0;
  }

  return res;
}

double gaussian_norm2(int n,double alp){
/** This is scalar product: <G(n,alp)|G(n,alp)>  */

  return (gaussian_int(2.0*n,2.0*alp));
}

double gaussian_norm1(int n,double alp){
/** This is scalar product: <G(n,alp)|G(n,alp)>  */

  return sqrt(gaussian_int(2.0*n,2.0*alp));
}

double gaussian_normalization_factor(int n,double alp){
/** This is scalar product: <G(n,alp)|G(n,alp)>   */

  return (1.0/sqrt(gaussian_int(2.0*n,2.0*alp)));
}





double FACTORIAL(int n) {
/** Factorial function :
  n! = n * (n-1) * (n-2) * ... * 1
  0! = 1  
*/

  double res = 1.0;
  if(n<0){ cout<<"Factorial of negative number is not defined!\n"; }
  else if(n==1||n==0){ res = 1.0; }
  else if(n==2){ res = 2.0; }
  else if(n==3){ res = 6.0; }
  else if(n==4){ res = 24.0; }
  else if(n==5){ res = 120.0; }
  else if(n==6){ res = 720.0; }
  else if(n==7){ res = 5040.0; }
  else if(n==8){ res = 40320.0; }
  else if(n==9){ res = 362880.0; }
  else if(n==10){ res = 3628800.0; }
  else{ res = n * FACTORIAL(n-1);  }

  return res;
}

double DFACTORIAL(int n) {
/** Double factorial :
  n!! = n * (n-2)!!, that is 1*3*5*... 

*/
  double res; 

  if(n<=1){ res = 1.0;}
  else if(n==3){ res = 3.0; }
  else if(n==5){ res = 15.0; }
  else if(n==7){ res = 105.0; }
  else if(n==9){ res = 945.0; }
  else if(n==11){ res = 10395.0; }
  else if(n==13){ res = 135135.0; }

  else{  res = n*DFACTORIAL(n-2); }

  return res;
}


double BINOM(int i,int n) {
/** Binomial coefficient : 
  C^i_n =  n! /( i! * (n-i)!) ,  0<=i<=n
*/

  if(i<0) { 
    cout<<"Error in BINOM: i must be non-negative\n"; 
    exit(0);
  }
  if(n<0) { 
    cout<<"Error in BINOM: n must be non-negative\n"; 
    exit(0);
  }
  if(i>n) { 
    cout<<"Error in BINOM: i must be not larger than n\n"; 
    exit(0);
  }


  double res;

  if(n==0){ 
    if(i==0){ res = 1.0; }
  }
  else if(n==1){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 1.0; }
  }
  else if(n==2){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 2.0; }
    else if(i==2){ res = 1.0; }
  }
  else if(n==3){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 3.0; }
    else if(i==2){ res = 3.0; }
    else if(i==3){ res = 1.0; }
  }
  else if(n==4){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 4.0; }
    else if(i==2){ res = 6.0; }
    else if(i==3){ res = 4.0; }
    else if(i==4){ res = 1.0; }
  }
  else if(n==5){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 5.0; }
    else if(i==2){ res = 10.0; }
    else if(i==3){ res = 10.0; }
    else if(i==4){ res = 5.0; }
    else if(i==5){ res = 1.0; }
  }
  else if(n==6){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 6.0; }
    else if(i==2){ res = 15.0; }
    else if(i==3){ res = 20.0; }
    else if(i==4){ res = 15.0; }
    else if(i==5){ res = 6.0; }
    else if(i==6){ res = 1.0; }
  }

  else{

    double a = FACTORIAL(n);
    double b = FACTORIAL(i);
    double c = FACTORIAL(n-i);

    res = (a/(b*c));
  }

  return res;
}


void zero_array(double* X,int n){
  for(int i=0;i<n;i++){ X[i] = 0.0; }
}

boost::python::list binomial_expansion(int n1,int n2,double x1,double x2,int is_derivs){

  int n = n1+n2+1;
  double* f; f = new double[n];
  double* dfdx1; dfdx1 = new double[n];
  double* dfdx2; dfdx2 = new double[n];

  binomial_expansion(n1,n2,x1,x2,f,dfdx1,dfdx2,is_derivs);

  boost::python::list lf; 
  boost::python::list ldfdx1; 
  boost::python::list ldfdx2; 

  for(int i=0;i<n;i++){
    lf.append(f[i]);

    if(is_derivs){
      ldfdx1.append(dfdx1[i]);
      ldfdx2.append(dfdx2[i]);
    }
  }

  boost::python::list res;   
  res.append(lf);
  if(is_derivs){
    res.append(ldfdx1);
    res.append(ldfdx2);
  }

  delete [] f;
  delete [] dfdx1;
  delete [] dfdx2;

  return res;

}

void binomial_expansion(int n1,int n2,double x1,double x2,double* f,double* dfdx1, double* dfdx2, int is_derivs){
/********************************************************************************
 This function calculates coefficients f_i in the series expansion

   (x+x1)^n1 * (x+x2)^n2 = summ x^i * f_i (x1,x2,n1,n2)
                             i
   f_i coefficients of the expansion

   dfdx1 - derivatives of f w.r.t. x1:  df / dx1
   dfdx2 - derivatives of f w.r.t. x2:  df / dx2

   is_derivs - defines if we want to compute derivatives 
*********************************************************************************/
  // Maximal expansion power
  int n = n1 + n2; 

  // Zero arrays
  zero_array(f,n+1);
  zero_array(dfdx1,n+1);
  zero_array(dfdx2,n+1);

  // Low-order expansions are set explicitly
  if(n<=4){

    //>>>>>>>>>>>>>>>>>>>>>>>> n == 0  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if(n1==0 && n2==0){  
      f[0] = 1; 

      if(is_derivs){
      dfdx1[0] = 0.0;
      dfdx2[0] = 0.0;
      }
    }/// n1==0 && n2==0

    //>>>>>>>>>>>>>>>>>>>>>>>> n == 1  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    else if(n1==1 && n2==0){ 
      f[0] = x1;       f[1] = 1.0;    

      if(is_derivs){
      dfdx1[0] = 1.0;  dfdx1[1] = 0.0;
      dfdx2[0] = 0.0;  dfdx2[1] = 0.0;
      }    
    }/// n1==1 && n2==0

    else if(n1==0 && n2==1){ 
      f[0] = x2;       f[1] = 1.0;    

      if(is_derivs){
      dfdx1[0] = 0.0;  dfdx1[1] = 0.0;
      dfdx2[0] = 1.0;  dfdx2[1] = 0.0; 
      }
    }/// n1==0 && n2==1

    //>>>>>>>>>>>>>>>>>>>>>>>> n == 2  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    else if(n1==1 && n2==1){ 
      f[0] = x1*x2;    f[1] = x1+x2;     f[2] = 1.0;    

      if(is_derivs){
      dfdx1[0] = x2;   dfdx1[1] = 1.0;   dfdx1[2] = 0.0;
      dfdx2[0] = x1;   dfdx2[1] = 1.0;   dfdx2[2] = 0.0; 
      }
    }/// n1==1 && n2==1

    else if(n1==2 && n2==0){
      // Optimized version: 2 mults
      f[0] = x1*x1;       f[1] = 2.0*x1;     f[2] = 1.0;    

      if(is_derivs){
      dfdx1[0] = f[1];    dfdx1[1] = 2.0;    dfdx1[2] = 0.0;
      dfdx2[0] = 0.0;     dfdx2[1] = 0.0;    dfdx2[2] = 0.0;
      }         
/**  Fully expanded version - retain for clarity:  3  mults
      f[0] = x1*x1;    dfdx1[0] = 2.0*x1;  dfdx2[0] = 0.0;
      f[1] = 2.0*x1;   dfdx1[1] = 2.0;     dfdx2[1] = 0.0;
      f[2] = 1.0;      dfdx1[2] = 0.0;     dfdx2[2] = 0.0;
**/
    }/// n1==2 && n2==0

    else if(n1==0 && n2==2){
      // Optimized version: 2 mults
      f[0] = x2*x2;       f[1] = 2.0*x2;     f[2] = 1.0;    

      if(is_derivs){
      dfdx1[0] = 0.0;     dfdx1[1] = 0.0;    dfdx1[2] = 0.0;
      dfdx2[0] = f[1];    dfdx2[1] = 2.0;    dfdx2[2] = 0.0;
      }      
/**  Fully expanded version - retain for clarity:  3  mults
      f[0] = x2*x2;    dfdx1[0] = 0.0;  dfdx2[0] = 2.0*x2;
      f[1] = 2.0*x2;   dfdx1[1] = 0.0;  dfdx2[1] = 2.0;
      f[2] = 1.0;      dfdx1[2] = 0.0;  dfdx2[2] = 0.0;
**/
    }/// n1==0 && n2==2

    //>>>>>>>>>>>>>>>>>>>>>>>> n == 3  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    else if(n1==3 && n2==0){ 
      // Optimized version: 5 mults
      double x11 = x1*x1;
      f[0] = x11*x1;       f[1] = 3.0*x11;     f[2] = 3.0*x1;   f[3] = 1.0;    

      if(is_derivs){
      dfdx1[0] = f[1];     dfdx1[1] = 6.0*x1;  dfdx1[2] = 3.0;  dfdx1[3] = 0.0;
      dfdx2[0] = 0.0;      dfdx2[1] = 0.0;     dfdx2[2] = 0.0;  dfdx2[3] = 0.0;
      }
/**  Fully expanded version - retain for clarity:  8  mults
      f[0] = x1*x1*x1;  dfdx1[0] = 3.0*x1*x1;  dfdx2[0] = 0.0;
      f[1] = 3.0*x1*x1; dfdx1[1] = 6.0*x1;     dfdx2[1] = 0.0;
      f[2] = 3.0*x1;    dfdx1[2] = 3.0;        dfdx2[2] = 0.0;
      f[3] = 1.0;       dfdx1[3] = 0.0;        dfdx2[3] = 0.0;
**/
    }/// n1==3 && n2==0

    else if(n1==0 && n2==3){ 
      // Optimized version: 5 mults
      double x22 = x2*x2;
      f[0] = x22*x2;       f[1] = 3.0*x22;      f[2] = 3.0*x2;    f[3] = 1.0;    

      if(is_derivs){
      dfdx1[0] = 0.0;      dfdx1[1] = 0.0;      dfdx1[2] = 0.0;   dfdx1[3] = 0.0;
      dfdx2[0] = f[1];     dfdx2[1] = 6.0*x2;   dfdx2[2] = 3.0;   dfdx2[3] = 0.0;
      }
/**  Fully expanded version - retain for clarity:  8  mults
      f[0] = x2*x2*x2;  dfdx1[0] = 0.0;        dfdx2[0] = 3.0*x2*x2;
      f[1] = 3.0*x2*x2; dfdx1[1] = 0.0;        dfdx2[1] = 6.0*x2;
      f[2] = 3.0*x2;    dfdx1[2] = 0.0;        dfdx2[2] = 3.0;
      f[3] = 1.0;       dfdx1[3] = 0.0;        dfdx2[3] = 0.0;
**/
    }/// n1==0 && n2==3

    else if(n1==2 && n2==1){ 
      // Optimized version: 7 mults
      double x12 = 2.0*x1*x2;
      double x11 = x1*x1;
      f[0] = x11*x2;  f[1] = x12 + x11;         f[2] = 2.0*x1+x2;  f[3] = 1.0; 

      if(is_derivs){ 
      dfdx1[0] = x12; dfdx1[1] = 2.0*(x2 + x1); dfdx1[2] = 2.0;    dfdx1[3] = 0.0; 
      dfdx2[0] = x11; dfdx2[1] = 2.0*x1;        dfdx2[2] = 1.0;    dfdx2[3] = 0.0; 
      }
/**  Fully expanded version - retain for clarity:  12  mults
      f[0] = x1*x1*x2;        dfdx1[0] = 2.0*x1*x2;     dfdx2[0] = x1*x1;
      f[1] = 2.0*x1*x2+x1*x1; dfdx1[1] = 2.0*x2+2.0*x1; dfdx2[1] = 2.0*x1;   
      f[2] = 2.0*x1+x2;       dfdx1[2] = 2.0;           dfdx2[2] = 1.0;      
      f[3] = 1.0;             dfdx1[3] = 0.0;           dfdx2[3] = 0.0;      
**/
    }/// n1==2 && n1==1

    else if(n1==1 && n2==2){ 
      // Optimized version: 7 mults
      double x12 = 2.0*x1*x2;
      double x22 = x2*x2;
      f[0] = x22*x1;   f[1] = x12 + x22;         f[2] = 2.0*x2 + x1;  f[3] = 1.0;     

      if(is_derivs){ 
      dfdx1[0] = x22;  dfdx1[1] = 2.0*x2;        dfdx1[2] = 1.0;      dfdx1[3] = 0.0; 
      dfdx2[0] = x12;  dfdx2[1] = 2.0*(x1 + x2); dfdx2[2] = 2.0;      dfdx2[3] = 0.0; 
      }

/**  Fully expanded version - retain for clarity:  12  mults
      f[0] = x2*x2*x1;        dfdx1[0] = x2*x2;      dfdx2[0] = 2.0*x2*x1;
      f[1] = 2.0*x2*x1+x2*x2; dfdx1[1] = 2.0*x2;     dfdx2[1] = 2.0*x1+2.0*x2;   
      f[2] = 2.0*x2+x1;       dfdx1[2] = 1.0;        dfdx2[2] = 2.0;      
      f[3] = 1.0;             dfdx1[3] = 0.0;        dfdx2[3] = 0.0;      
**/
    }/// n1==1 && n1==2

    //>>>>>>>>>>>>>>>>>>>>>>>> n == 4  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    else if(n1==4 && n2==0){ 
      // Optimized version: 8 mults
      double x11 = x1*x1;
      double x111 = 4.0*x11*x1;
      f[0] = x11*x11;  f[1] = x111;         f[2] = 6.0*x11;     f[3] = 4.0*x1;  f[4] = 1.0;     

      if(is_derivs){ 
      dfdx1[0] = x111; dfdx1[1] = 12.0*x11; dfdx1[2] = 12.0*x1; dfdx1[3] = 4.0; dfdx1[4] = 0.0; 
      dfdx2[0] = 0.0;  dfdx2[1] = 0.0;      dfdx2[2] = 0.0;     dfdx2[3] = 0.0; dfdx2[4] = 0.0; 
      }
/**  Fully expanded version - retain for clarity:  15 mults
      f[0] = x1*x1*x1*x1;   dfdx1[0] = 4.0*x1*x1*x1;   dfdx2[0] = 0.0;
      f[1] = 4.0*x1*x1*x1;  dfdx1[1] = 12.0*x1*x1;     dfdx2[1] = 0.0;
      f[2] = 6.0*x1*x1;     dfdx1[2] = 12.0*x1;        dfdx2[2] = 0.0;
      f[3] = 4.0*x1;        dfdx1[3] = 4.0;            dfdx2[3] = 0.0;
      f[4] = 1.0;           dfdx1[4] = 0.0;            dfdx2[4] = 0.0;
**/
    }/// n1==4 && n2==0

    else if(n1==0 && n2==4){ 
      // Optimized version: 8 mults
      double x22 = x2*x2;
      double x222 = 4.0*x22*x2;
      f[0] = x22*x22;  f[1] = x222;         f[2] = 6.0*x22;     f[3] = 4.0*x2;  f[4] = 1.0;     

      if(is_derivs){ 
      dfdx2[0] = x222; dfdx2[1] = 12.0*x22; dfdx2[2] = 12.0*x2; dfdx2[3] = 4.0; dfdx2[4] = 0.0; 
      dfdx1[0] = 0.0;  dfdx1[1] = 0.0;      dfdx1[2] = 0.0;     dfdx1[3] = 0.0; dfdx1[4] = 0.0; 
      }
/**  Fully expanded version - retain for clarity:  15 mults
      f[0] = x2*x2*x2*x2;   dfdx2[0] = 4.0*x2*x2*x2;   dfdx1[0] = 0.0;
      f[1] = 4.0*x2*x2*x2;  dfdx2[1] = 12.0*x2*x2;     dfdx1[1] = 0.0;
      f[2] = 6.0*x2*x2;     dfdx2[2] = 12.0*x2;        dfdx1[2] = 0.0;
      f[3] = 4.0*x2;        dfdx2[3] = 4.0;            dfdx1[3] = 0.0;
      f[4] = 1.0;           dfdx2[4] = 0.0;            dfdx1[4] = 0.0;
**/
    }/// n1==0 && n2==4

    else if(n1==3 && n2==1){ 
      // Optimized version: 15 mults
      double x11 = x1*x1;
      double x22 = x2*x2;
      double x12 = x1*x2;
      f[0] = x11*x12;           f[1] = x11*x1 + 3.0*x11*x2;   f[2] = 3.0*(x12 + x11);     f[3] = 3.0*x1+x2; f[4] = 1.0;     

      if(is_derivs){ 
      dfdx2[0] = x11*x1;        dfdx2[1] = 3.0*x11;            dfdx2[2] = 3.0*x1;          dfdx2[3] = 1.0;   dfdx2[4] = 0.0; 
      dfdx1[0] = f[1]-dfdx2[0]; dfdx1[1] = dfdx2[1] + 6.0*x12; dfdx1[2] = 3.0*x2 + 6.0*x1; dfdx1[3] = 3.0;   dfdx1[4] = 0.0; 
      }
/**  Fully expanded version - retain for clarity:  27 mults
      f[0] = x1*x1*x1*x2;             dfdx1[0] = 3.0*x1*x1*x2;            dfdx2[0] = x1*x1*x1;
      f[1] = x1*x1*x1+3.0*x1*x1*x2;   dfdx1[1] = 3.0*x1*x1+6.0*x1*x2;     dfdx2[1] = 3.0*x1*x1;
      f[2] = 3.0*x1*x2+3.0*x1*x1;     dfdx1[2] = 3.0*x2+6.0*x1;           dfdx2[2] = 3.0*x1;
      f[3] = 3.0*x1+x2;               dfdx1[3] = 3.0;                     dfdx2[3] = 1.0;
      f[4] = 1.0;                     dfdx1[4] = 0.0;                     dfdx2[4] = 0.0;
**/
    }/// n1==3 && n2==1

    else if(n1==1 && n2==3){ 
      // Optimized version: 15 mults
      double x11 = x1*x1;
      double x22 = x2*x2;
      double x12 = x1*x2;
      f[0] = x22*x12;           f[1] = 3.0*x22*x1 + x22*x2;   f[2] = 3.0*(x12 + x22);     f[3] = 3.0*x2+x1; f[4] = 1.0;     

      if(is_derivs){ 
      dfdx1[0] = x22*x2;        dfdx1[1] = 3.0*x22;            dfdx1[2] = 3.0*x2;          dfdx1[3] = 1.0;   dfdx1[4] = 0.0; 
      dfdx2[0] = f[1]-dfdx1[0]; dfdx2[1] = dfdx1[1] + 6.0*x12; dfdx2[2] = 3.0*x1 + 6.0*x2; dfdx2[3] = 3.0;   dfdx2[4] = 0.0; 
      }
/**  Fully expanded version - retain for clarity: 27 mults
      f[0] = x2*x2*x2*x1;             dfdx2[0] = 3.0*x2*x2*x1;            dfdx1[0] = x2*x2*x2;
      f[1] = x2*x2*x2+3.0*x2*x2*x1;   dfdx2[1] = 3.0*x2*x2+6.0*x2*x1;     dfdx1[1] = 3.0*x2*x2;
      f[2] = 3.0*x2*x1+3.0*x2*x2;     dfdx2[2] = 3.0*x1+6.0*x2;           dfdx1[2] = 3.0*x2;
      f[3] = 3.0*x2+x1;               dfdx2[3] = 3.0;                     dfdx1[3] = 1.0;
      f[4] = 1.0;                     dfdx2[4] = 0.0;                     dfdx1[4] = 0.0;
**/
    }/// n1==1 && n2==3

    else if(n1==2 && n2==2){ 
      // Optimized version: 18 mults
      double x11 = x1*x1;
      double x22 = x2*x2;
      double x12 = x1*x2;
      f[0] = x11*x22;           f[1] = 2.0*(x1*x22 + x11*x2);   f[2] = x11 + x22 + 4.0*x12; f[3] = 2.0*(x1 + x2); f[4] = 1.0;     

      if(is_derivs){ 
      dfdx1[0] = 2.0*x1*x22;    dfdx1[1] = 2.0*x22 + 4.0*x12; dfdx1[2] = 2.0*x1 + 4.0*x2; dfdx1[3] = 2.0;       dfdx1[4] = 0.0; 
      dfdx2[0] = f[1]-dfdx1[0]; dfdx2[1] = 4.0*x12 + 2.0*x11; dfdx2[2] = 2.0*x2 + 4.0*x1; dfdx2[3] = 2.0;       dfdx2[4] = 0.0; 
      }
/** Fully expanded verion - retain for clarity: 33 mults
      f[0] = x1*x1*x2*x2;                 dfdx1[0] = 2.0*x1*x2*x2;            dfdx2[0] = 2.0*x1*x1*x2;
      f[1] = 2.0*x1*x2*x2+2.0*x2*x1*x1;   dfdx1[1] = 2.0*x2*x2+4.0*x2*x1;     dfdx2[1] = 4.0*x1*x2+2.0*x1*x1;
      f[2] = x1*x1+x2*x2+4.0*x1*x2;       dfdx1[2] = 2.0*x1+4.0*x2;           dfdx2[2] = 2.0*x2+4.0*x1;
      f[3] = 2.0*x1+2.0*x2;               dfdx1[3] = 2.0;                     dfdx2[3] = 2.0;
      f[4] = 1.0;                         dfdx1[4] = 0.0;                     dfdx2[4] = 0.0;
**/
    }/// n1==2 && n2==2


  }else{ // general expansion

    for(int i=0;i<=n1;i++){
      for(int j=0;j<=n2;j++){
        double b = BINOM(i,n1)*BINOM(j,n2);
        double p1 = FAST_POW(x1,n1-i);
        double p2 = FAST_POW(x2,n2-j);
        //
        f[i+j] += b*p1*p2;

        if(is_derivs){ 

          double pw1,pw2;
          if(n1-i==0){ pw1 = 0.0; } else{ pw1 = (n1-i)*FAST_POW(x1,n1-i-1); }
          if(n2-j==0){ pw2 = 0.0; } else{ pw2 = (n2-j)*FAST_POW(x2,n2-j-1); }

          dfdx1[i+j] += b*pw1*p2;
          dfdx2[i+j] += b*p1*pw2;

        }// is_derivs
      }// for j
    }// for i
  }// general expansion

}// binomial_expansion







void LEGENDRE(int n,double x,double a,double b,double& p,double& q){
/** Legendre polynomials and its derivatives
 Algorithm is taken from: Flowers, B. H. "An Introduction to Numerical Methods in C++", Clarendon Press, Oxford, 1995

 <> - subscript
 recursive relation: (n+1)P<n+1>(x) = (2n+1)xP<n>(x) - nP<n-1>(x)
 normaliztion:        = 2/(2n+1)
 recursive derivatives: P'<n+1>(x) = xP'<n>(x) + (n+1)P<n>(x)    
 weighting function: w(x) = 1
        
 defined for an interval [a,b]

*/
    if(x<a){ 
      cout<<"Error: Legendre polynomial is not defined for x = "<<x<<" < a = "<<a<<endl; 
      exit(0);
    }
    if(x>b){ 
      cout<<"Error: Legendre polynomial is not defined for x = "<<x<<" > b = "<<b<<endl; 
      exit(0);
    }

    x = 0.5*((b+a)+(b-a)*x);

    double r = 0, s = 0, t = 0;

    // r - value of polynomial one degree  less
    // t - value of polynomial two degrees less
    // s - differential coefficient corresponding to r

    p = 1; q = 0;

    for(float m=0;m<n;m++){
        r = p; s = q;
        p = (2.0*m+1)*x*r - m*t;
        p /= m+1;
        q = x*s + (m+1)*r;
        t = r;
    }

}

void CHEBYSHEV(int n,double x,double a,double b,double& p,double& q){
/** Chebyshev polynomials and its derivatives    

 <> - subscript

 recursive relation: T<n+1>(x) = 2xT<n>(x) - T<n-1>(x)
 normaliztion:        = PI (n==0);  PI/2 (n>0)
 recursive derivatives: T'<n+1>(x) = [(n+1)x/n]T'<n>(x) + (n+1)T<n>(x)   
 weighting function: w(x) = 1/sqrt(1-x^2)
        
 defined for an interval [a,b]

*/
    if(x<a){ 
      cout<<"Error: Chebychev polynomial is not defined for x = "<<x<<" < a = "<<a<<endl; 
      exit(0);
    }
    if(x>b){ 
      cout<<"Error: Chebychev polynomial is not defined for x = "<<x<<" > b = "<<b<<endl; 
      exit(0);
    }


    x = 0.5*((b+a)+(b-a)*x);

    double r = 0, s = 0, t = x;

    // r - value of polynomial one degree  less
    // t - value of polynomial two degrees less
    // s - differential coefficient corresponding to r

    p = 1; q = 0;

    for(float m=0;m<n;m++){

        r = p; s = q;
        p = 2.0*x*r - t;
        if(m==0){  q = (m+1)*r;  }
        else{ q = ((m+1)/m)*x*s + (m+1)*r;  }
        t = r;
    }

}

void LAGUERRE(int n,double x,double& p,double& q){
/** Laguerre polynomials and its derivatives     

 <> - subscript

 recursive relation: (n+1)L<n+1>(x) = (2n+1-x)L<n>(x) - nL<n-1>(x)
 normaliztion:        = 1
 recursive derivatives: L'<n+1>(x) = L'<n>(x) - L<n>(x)  
 weighting function: w(x) = exp(-x)
        
 defined for an interval [0,+inf)

*/

    if(x<0.0){ 
      cout<<"Error: Laguerre polynomial is not defined for x = "<<x<<" < 0.0"<<endl; 
      exit(0);
    }


    double r = 0, s = 0, t = 0;

    // r - value of polynomial one degree  less
    // t - value of polynomial two degrees less
    // s - differential coefficient corresponding to r

    p = 1; q = 0;

    for(float m=0;m<n;m++){

      r = p; s = q;
      p = (2.0*m+1-x)*r - m*t;
      p /= m+1;
      q = s - r;
      t = r;
    }

}


void HERMITE(int n,double x,double& p,double& q){

        // Hermite polynomials and its derivatives      
        /*

        <> - subscript

                recursive relation: H<n+1>(x) = 2xH<n>(x) - 2nH<n-1>(x)
                normaliztion:        = 2^n*n!*sqrt(PI)
        recursive derivatives: H'<n+1>(x) = 2(n+1)H<n>(x)       
                weighting function: w(x) = exp(-x^2)
        
            defined for an interval (-inf,+inf)

        */

        double r = 0, s = 0, t = 0;

        // r - value of polynomial one degree  less
        // t - value of polynomial two degrees less
        // s - differential coefficient corresponding to r

        p = 1; q = 0;

        for(float m=0;m<n;m++){

                r = p; s = q;
                p = 2.0*x*r - 2.0*m*t;
                q = 2.0*(m+1)*r;
                t = r;
        }

}

double Ellipe(double phi0,double k,int N)  // Elliptic integral of first kind (Jacobi elliptic function)
{

        double cp,sp,sp2;
        double k2 = k*k;
        double phi=phi0;
        double S,M,M1,ADD;

        cp = cos(phi);
        sp = sin(phi);
        sp2= sp*sp;

        S  = 0.5*(phi-cp*sp);
        M  = 0.5*k2;
        M1 = cp*sp*sp2;

        double sum=phi;

        ADD = M*S;

        for(int n=1;n<N;n++){


                sum+=ADD;

                S = (1.0/(2.0*n+2.0))*(-M1+(2.0*n+1.0)*S);

                M1 = M1*sp2;

                M *= ((2.0*(n+1.0)-1.0)*k2/(2.0*(n+1.0)));

                ADD = M*S;


        }

        //sum = phi/(div*a1);

        return sum;
}
void Ellipe2(double u,double k,double& sn,double& cn,double& dn){

        double k2,k4,k6;
        double u2,u4,u6,u8;
        double f2,f3,f4,f5,f6,f7,f8;

        f2 = 0.5;
        f3 = f2/3.0;
        f4 = f3/4.0;
        f5 = f4/5.0;
        f6 = f5/6.0;
        f7 = f6/7.0;
        f8 = f7/8.0;

        u2 = u*u;
        u4 = u2*u2;
        u6 = u2*u4;
        u8 = u2*u6;

        k2 = k*k;
        k4 = k2*k2;
        k6 = k2*k4;


        sn = u * (1.0 - (1.0+k2)*u2*f3 + (1.0+14.0*k2+k4)*u4*f5 - (1.0+135.0*k2+135.0*k4+k6)*u6*f7);

    cn = sqrt(1.0-sn*sn);

        dn = sqrt(1.0-k*k*sn*sn);

}
void Jacobi_Elliptic(double u,double m,double tol,double& am, double& sn, double& cn, double& dn){

        // [A. S.] = M. Abramowitz, I.A. Stegun, Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables, Dover, New
    //York, 1965

        double K,Kprime;

        double k = sqrt(m);

        K      = Km(m,tol);
        Kprime = Km(1-m,tol);

        double q  = exp(-M_PI*(Kprime/K)); // nome
        double q2 = q*q;

        double v = (M_PI*u/(2.0*K));

        double q_sqrt = sqrt(q);
        double q_n, q_2n;
        double term,term1,term2,term3;
        double mt1, mt2, mt3;

        am = v;

        q_n  = q;
        q_2n = q2;

        double n = 1;
        while(true){

                // [A. S. 17.3.24]
                term = 2.0*q_n*sin(2.0*n*v)/(n*(1+q_2n));

                am += term;

                n++;
                q_n *= q;
                q_2n *= q2;

                if(fabs(term)<tol){ break; }

        }

        sn = sin(am);          // [A.S. 17.2.2]
        cn = cos(am);          // [A.S. 17.2.3]
        dn = sqrt(1-m*sn*sn);  // [A.S. 17.2.4]

}
double Km(double m,double tol){

        // Return function:
        // K(m) = F(PI/2|m) = complete elliptic integral of first kind;

        double a0,b0,c0;
        double a1,b1,c1;
        double k   = sqrt(m); // m = k*k;


        a1 = a0 = 1;
        b1 = b0 = sqrt(1-m);
        c1 = c0 = sqrt(m);

        double n = 1;

        while(true){

                a1 = 0.5*(a0+b0);
                b1 = sqrt(a0*b0);
                c1 = 0.5*(a0-b0);

                if(fabs(c1*n)<tol) break;

                a0 = a1;
                b0 = b1;
                c0 = c1;

                n++;

        }

//        cout<<"m = "<<m<<endl;
//        cout<<"# of iterations = "<<n<<endl;

        return (0.5*M_PI/a1);

}

void Ellint(double m,double sinphi,double tol,double& Km_,double& value){
// This function calculates incomplete elliptic integral of first kind (returned in value)
// It also returns the complete elliptic integral of first kind (in Km_)
// Arguments: sinphi - sine of angle phi
// m - parameter, tol - accuracy 

//------ Precompute some quantities and the complete integral -------------
  if (m < 0.0) {
    fprintf(stderr, "ERROR: m < 0.0 in Ellint (m = %f): FIX IT!\n", m);
    exit(1);
  }

  int MAXN = 25;
  double* a; a = new double[MAXN];
  double* b; b = new double[MAXN];
  double* d; d = new double[MAXN];

  int N = 0;
  a[N] = 1.0;
  b[N] = sqrt(1.0-m);
  d[N] = b[N]/a[N];
  double err;

  while(N<MAXN){
    N++;
    a[N] = 0.5*(a[N-1] + b[N-1]);
    b[N] = sqrt(a[N-1]*b[N-1]);
    d[N] = b[N]/a[N];
    err = (a[N-1]-b[N-1]);
    
    if(fabs(err*(double)N)<tol) break;
  }

  Km_      = (0.5*M_PI/a[N]);
//  cout<<"m = "<<m<<endl;
//  cout<<"# of iterations = "<<N<<endl;

  double Kscale = (1.0/((1<<N)*a[N]));

//---------- Now calculate incomplete integral -----------
  if (sinphi <= -1.0){   value =  -Km_; }
  else if (sinphi >= 1.0){   value = Km_; }
  else{
    double t = (sinphi/sqrt(1.0-sinphi*sinphi));
    int s = 0;
    for (int i = 0; i < N; i++) {
      double c = 1.0 - d[i]*t*t;
      s *= 2;
      if (c < 0.0){ s += (t>0.0?1:-1);}
      t = (1.0 + d[i])*t/c;
    }
    value = Kscale*(atan(t)+s*M_PI);
  }// else

  delete [] a;
  delete [] b;
  delete [] d;

}


int randperm(int size,int of_size,vector<int>& result){
// Makes a random permutation of 'of_size' numbers and places 'size'
// of them into the 'result' vector

  int* all_range;

  all_range = new int[of_size];

  // Initialize
  for(int i=0;i<of_size;i++){
     all_range[i] = i;
  }
  // Make permutations
  int indx;
  int lowest=0, highest=of_size-1;
  int range=(highest-lowest)+1;
  int tmp;

  for(int i=0;i<size;i++){

    indx = rand();
    indx = i+lowest+int((range-i)*(indx/(RAND_MAX + 1.0)));

    // Swap all_range[i] and all_range[indx]
    tmp = all_range[i];
    all_range[i] = all_range[indx];
    all_range[indx] = tmp;

  }
  // Now get 'size' first elements
  if(result.size()>0) {result.clear();}
  for(int i=0;i<size;i++){
     result.push_back(all_range[i]);
  }

  delete [] all_range;

  return 0;
}



MATRIX3x3 exp_(MATRIX3x3& x,double dt){
  // This function calculates exp(x*dt)
  // C.T() * x * C = D

  MATRIX x3x3(3,3);
  x3x3.set(0,0, x.xx); x3x3.set(0,1, x.xy); x3x3.set(0,2, x.xz);
  x3x3.set(1,0, x.yx); x3x3.set(1,1, x.yy); x3x3.set(1,2, x.yz);
  x3x3.set(2,0, x.zx); x3x3.set(2,1, x.zy); x3x3.set(2,2, x.zz);

  MATRIX res(3,3);
  MATRIX C(3,3);
  MATRIX D(3,3);

  libmeigen::solve_eigen(x3x3, D, C, 0);


  D.M[0] = exp(dt*D.M[0]);
  D.M[4] = exp(dt*D.M[4]);
  D.M[8] = exp(dt*D.M[8]);

  res = C * D * C.T();
  MATRIX3x3 res3x3;
  res3x3.xx = res.M[0];  res3x3.xy = res.M[1];  res3x3.xz = res.M[2];
  res3x3.yx = res.M[3];  res3x3.yy = res.M[4];  res3x3.yz = res.M[5];
  res3x3.zx = res.M[6];  res3x3.zy = res.M[7];  res3x3.zz = res.M[8];

  return res3x3;
}


MATRIX exp_(MATRIX& x, double dt){
/**
  This function computes exp(x*dt) for a given matrix x via:  C * x * C.T() = D

  \param[in] x input matrix
  \param[in] dt scaling factor

*/

  if(x.n_cols != x.n_rows){
    cout<<"Error in libspecialfunctions::exp_ : the input matrix is not square\n"; exit(0); 
  }

  int i,j;
 
  // Let us first diagonalize the input matrix x
  int sz = x.n_cols;  
  MATRIX* evec; evec = new MATRIX(sz, sz);  *evec = 0.0;
  MATRIX* eval; eval = new MATRIX(sz, sz);  *eval = 0.0;
  MATRIX res(sz,sz);

  // Find the eigenvalues of the the S matrix
  libmeigen::solve_eigen(x, *eval, *evec, 0);  // x * evec = evec * eval  ==>  x = evec * eval * evec.T()

  for(i=0;i<sz;i++){ res.M[i*sz+i]= std::exp(dt * eval->get(i,i)); }

  // Convert to the original basis
  res = (*evec) * res * ((*evec).T());

  delete eval;
  delete evec;

  return res;
}


CMATRIX exp_(CMATRIX& x, complex<double> dt){
/**
  This function computes exp(x*dt) for a given matrix x via:  C.T() * x * C = D

  \param[in] x input matrix
  \param[in] dt scaling factor

*/

  if(x.n_cols != x.n_rows){
    cout<<"Error in libspecialfunctions::exp_ : the input matrix is not square\n"; exit(0); 
  }

  int i,j;
 
  // Let us first diagonalize the input matrix x
  int sz = x.n_cols;  
  CMATRIX* evec; evec = new CMATRIX(sz, sz);  *evec = complex<double>(0.0, 0.0);
  CMATRIX* eval; eval = new CMATRIX(sz, sz);  *eval = complex<double>(0.0,0.0);
  CMATRIX res(sz,sz);

  // Find the eigenvalues of the the S matrix
  libmeigen::solve_eigen(x, *eval, *evec, 0);  // x * evec = evec * eval  ==>  x = evec * eval * evec.H()

  for(i=0;i<sz;i++){ res.M[i*sz+i]= std::exp(dt * eval->get(i,i)); }

  // Convert to the original basis
  res = (*evec) * res * ((*evec).H());

  delete eval;
  delete evec;

  return res;
}



MATRIX exp_2(MATRIX& x, double dt, int nterms, double max_tol){
/**
  This function computes exp(x*dt) for a given matrix x via the Taylor series sum

  \param[in] x input matrix
  \param[in] dt scaling factor

*/

  if(x.n_cols != x.n_rows){
    cout<<"Error in libspecialfunctions::exp_ : the input matrix is not square\n"; exit(0); 
  }

  int sz = x.n_cols;  
  MATRIX tmp(sz,sz);
  MATRIX res(sz,sz);

  res.identity();
  tmp.identity();
  double err = 1e+10;
  
  int n = 1;
  while(n < (nterms+1) && err > max_tol){

    tmp = x * tmp;
    tmp *= (dt/(double)n);
    res += tmp;

    err = tmp.max_elt();
    n++;
  }


  if(err > max_tol){
      cout<<"Error in exp_2(MATRIX...): The the sum has not converged to the defined tolerance of "<<max_tol
          <<" in "<<nterms<<" iterations.\nExiting now...\n";
      exit(0);
  }

  return res;
}


MATRIX exp_2(MATRIX& x, double dt, int nterms){

  return exp_2(x, dt, nterms, 0.0);

}

MATRIX exp_2(MATRIX& x, double dt){

  return exp_2(x, dt, 1000, 0.0);

}




CMATRIX exp_2(CMATRIX& x, complex<double> dt, int nterms, double max_tol){
/**
  This function computes exp(x*dt) for a given matrix x via the Taylor series sum

  \param[in] x input matrix
  \param[in] dt scaling factor

*/

  if(x.n_cols != x.n_rows){
    cout<<"Error in libspecialfunctions::exp_ : the input matrix is not square\n"; exit(0); 
  }

  int sz = x.n_cols;  
  CMATRIX tmp(sz,sz);
  CMATRIX res(sz,sz);

  res.identity();
  tmp.identity();
  double err = 1e+10;
  
  int n = 1;
  while(n < (nterms+1) && err > max_tol){

    tmp = x * tmp;
    tmp *= (dt/(double)n);
    res += tmp;

    err = abs(tmp.max_elt());
    n++;
  }

  if(err > max_tol){
      cout<<"Error in exp_2(CMATRIX...): The the sum has not converged to the defined tolerance of "<<max_tol
          <<" in "<<nterms<<" iterations.\nExiting now...\n";
      exit(0);
  }


  return res;
}


CMATRIX exp_2(CMATRIX& x, complex<double> dt, int nterms){

  return exp_2(x, dt, nterms, 0.0);

}

CMATRIX exp_2(CMATRIX& x, complex<double> dt){

  return exp_2(x, dt, 1000, 0.0);

}





MATRIX3x3 exp1_(MATRIX3x3& x,double dt){
  // This function calculates exp(x*dt)*sinh(x*t)/(x*t)
  // where 1/x is x^-1
  // C.T() * x * C = D
  MATRIX x3x3(3,3);
  x3x3.set(0,0, x.xx); x3x3.set(0,1, x.xy); x3x3.set(0,2, x.xz);
  x3x3.set(1,0, x.yx); x3x3.set(1,1, x.yy); x3x3.set(1,2, x.yz);
  x3x3.set(2,0, x.zx); x3x3.set(2,1, x.zy); x3x3.set(2,2, x.zz);

  MATRIX res(3,3);
  MATRIX C(3,3);
  MATRIX D(3,3);

  libmeigen::solve_eigen(x3x3, D, C, 0);


  D.M[0] = exp(dt*D.M[0])*sinh_(dt*D.M[0]);
  D.M[4] = exp(dt*D.M[4])*sinh_(dt*D.M[4]);
  D.M[8] = exp(dt*D.M[8])*sinh_(dt*D.M[8]);


  res = C * D * C.T();
  MATRIX3x3 res3x3;
  res3x3.xx = res.M[0];  res3x3.xy = res.M[1];  res3x3.xz = res.M[2];
  res3x3.yx = res.M[3];  res3x3.yy = res.M[4];  res3x3.yz = res.M[5];
  res3x3.zx = res.M[6];  res3x3.zy = res.M[7];  res3x3.zz = res.M[8];

  return res3x3;

}


MATRIX exp1_(MATRIX& x,double dt){
  // This function calculates exp(x*dt)*sinh(x*t)/(x*t)
  // where 1/x is x^-1
  // C.T() * x * C = D
  int sz = x.n_cols;

  MATRIX res(sz,sz);
  MATRIX C(sz,sz);
  MATRIX D(sz,sz);

  libmeigen::solve_eigen(x, D, C, 0);


  D.M[0] = exp(dt*D.M[0])*sinh_(dt*D.M[0]);
  D.M[4] = exp(dt*D.M[4])*sinh_(dt*D.M[4]);
  D.M[8] = exp(dt*D.M[8])*sinh_(dt*D.M[8]);

  res = C * D * C.T();

  return res;
}


//============================ Random number generators =================
double RANDOM(double a,double b){

        double u;
               u = rand();
                   u =u/RAND_MAX;

        return (a+(b-a)*u);
}


int merge_sort(vector< pair<int,double> >& in, vector< pair<int,double> >& out){

  if(out.size()>0){ out.clear(); }
  int sz = in.size();
  if(sz==0){ }
  else if(sz==1){ out = in; }
  else{
    // Divide in into 2 blocks of approximately same size each
    int half = sz/2;
    vector< pair<int,double> > in1,in2,out1,out2;
    for(int i=0;i<half;i++){ in1.push_back(in[i]); }
    merge_sort(in1,out1);
       
    for(int i=half;i<sz;i++){ in2.push_back(in[i]); }
    merge_sort(in2,out2); 
      
    // Now merge two parts
    int cl,cr; cl = 0; cr = half;
    while((cl<half) && (cr<sz)){
 
    /// EXTREMELY IMPORTANT !!!  This simple, slight difference - the use of < or <= makes HUGE differnece
    /// The "good" version maximally preserves the ordering of orbitals, so one does not run into trouble of alternating
    /// charges - this also leads to symmetric charge distribution in unrestricted formulations even without population smearing
    /// I think it is even more than that - this can lead to convergence (or faster convergence), while the wrong
    /// method may lead to either non-convergent scheme or to sifnificantly slower convergenc.

    if(out1[cl].second<=out2[cr-half].second){ out.push_back(out1[cl]); cl++; }  ///< <-- This is good
//    if(out1[cl].second<out2[cr-half].second){ out.push_back(out1[cl]); cl++; }  ///< <-- Try old

    /// The "bad" version will alternate order of nearby orbitals
    /// It is here only for the purpose of "demonstration of pathological implementation"
//    if(out1[cl].second<out2[cr-half].second){ out.push_back(out1[cl]); cl++; } <-- This is BAD
      else{ out.push_back(out2[cr-half]); cr++; }
    } 
    while(cl<half){ out.push_back(out1[cl]); cl++;}
    while(cr<sz)  { out.push_back(out2[cr-half]); cr++;}
  }
  return 0;

}// int merge_sort(vector< pair<int,double> >& in, vector< pair<int,double> >& out)


boost::python::list merge_sort(boost::python::list inp){
  int i;
  int sz = len(inp);
  vector< pair<int,double> > inp1, out;

  for(i=0;i<sz;i++){
     int indx = boost::python::extract<int>(inp[i][0]);
     double val = boost::python::extract<double>(inp[i][1]);

     inp1.push_back(pair<int, double>(indx, val));
  }


  merge_sort(inp1, out);


  boost::python::list res;

  for(i=0;i<sz;i++){
    boost::python::list ri;
    ri.append(out[i].first);
    ri.append(out[i].second);
  
    res.append(ri);
  }
 
  return res;

}



int merge_sort(vector< double >& in, vector< double >& out){

  if(out.size()>0){ out.clear(); }
  int sz = in.size();
  if(sz==0){ }
  else if(sz==1){ out = in; }
  else{

    out = vector<double>(sz, 0.0);

    // Divide in into 2 blocks of approximately same size each
    int half = sz/2;
    vector< double > in1(half, 0.0);
    vector< double > out1(half, 0.0);
    vector< double > in2(sz-half, 0.0);
    vector< double > out2(sz-half, 0.0);

    for(int i=0;i<half;i++){ in1[i] = in[i]; }
    merge_sort(in1,out1);
       
    for(int i=half;i<sz;i++){ in2[i-half] = in[i]; }
    merge_sort(in2,out2); 
      
    // Now merge two parts
    int cl,cr, indx; cl = 0; cr = half; indx = 0;
    while((cl<half) && (cr<sz)){
 
    /// EXTREMELY IMPORTANT !!!  This simple, slight difference - the use of < or <= makes HUGE differnece
    /// The "good" version maximally preserves the ordering of orbitals, so one does not run into trouble of alternating
    /// charges - this also leads to symmetric charge distribution in unrestricted formulations even without population smearing
    /// I think it is even more than that - this can lead to convergence (or faster convergence), while the wrong
    /// method may lead to either non-convergent scheme or to sifnificantly slower convergenc.

      if(out1[cl] <= out2[cr-half] ){ out[indx] = out1[cl]; cl++; indx++; }  ///< <-- This is good

    /// The "bad" version will alternate order of nearby orbitals
    /// It is here only for the purpose of "demonstration of pathological implementation"
      else{ out[indx] = out2[cr-half]; cr++; indx++; }
    } 
    while(cl<half){ out[indx] = out1[cl]; cl++; indx++; }
    while(cr<sz)  { out[indx] = out2[cr-half]; cr++; indx++; }
  }
  return 0;

}// int merge_sort(vector< double >& in, vector< double >& out)


vector< vector<int> > permutations_reiteration(vector<int>& given_list, int size, int num_elements, vector< vector<int> >& list_of_permutations){
    /*
    This function generates all permutations of a given list and outputs as a list of lists.
    The heart of this algorithm is inspired by:
    https://www.geeksforgeeks.org/heaps-algorithm-for-generating-permutations/
    */

    vector<int> one_permutation;
    int i;

    if(size == 1){
        one_permutation.clear();
        for(i = 0; i<num_elements; i++)
            one_permutation.push_back(given_list[i]);
        list_of_permutations.push_back(one_permutation);
        return list_of_permutations;
    } // if

    for(i=0;i<size;i++){
        list_of_permutations = permutations_reiteration(given_list, size-1, num_elements, list_of_permutations);

        if(size%2==0) // if even
            swap(given_list[0], given_list[size-1]);

        else // if odd
            swap(given_list[i], given_list[size-1]);
} // for loop
    return list_of_permutations;

}

vector< vector<int> > compute_all_permutations(vector<int>& given_list){
    /*
    This function acts as a wrapper for the permutations_reiteration function,
    simplifying the required user input to only a list.
    */
    int num_elements;
    vector< vector<int> > list_of_permutations;
    list_of_permutations.empty();
    num_elements= given_list.size();

    return permutations_reiteration(given_list, num_elements, num_elements, list_of_permutations);
}

}// namespace libspecialfunctions
}// namespace liblibra



