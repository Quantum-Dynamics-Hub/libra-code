#include "PrimitiveS.h"

//----------------------- Members of PrimitiveS class -----------------------
// Copy constructor
PrimitiveS::PrimitiveS(const PrimitiveS& g){

  // Default parameters
  // Simple S-type function
  n = 0;         is_n = 1;
  l = 0;         is_l = 1;
  m = 0;         is_m = 1;
  alpha = 1.0;   is_alpha = 1;
                 is_center = 0;

  if(g.is_n)   {  n = g.n; is_n = 1; }
  if(g.is_l)   {  l = g.l; is_l = 1; }
  if(g.is_m)   {  m = g.m; is_m = 1; }
  if(g.is_alpha) {  alpha = g.alpha; is_alpha = 1; }
  if(g.is_center){
    center = new VECTOR;  
    center = g.center; is_center = 1;
  }

}


PrimitiveS& PrimitiveS::operator=(const PrimitiveS& g){

  // Default parameters
  // Simple S-type function
  n = 0;         is_n = 1;
  l = 0;         is_l = 1;
  m = 0;         is_m = 1;
  alpha = 1.0;   is_alpha = 1;
                 is_center = 0;

  if(g.is_n)   {  n = g.n; is_n = 1; }
  if(g.is_l)   {  l = g.l; is_l = 1; }
  if(g.is_m)   {  m = g.m; is_m = 1; }
  if(g.is_alpha) {  alpha = g.alpha; is_alpha = 1; }
  if(g.is_center){
//    center = new VECTOR;  
    center = g.center; is_center = 1;
  }

  return *this;

}

void PrimitiveS::show_info(){

  std::cout<<"PrimitiveS properties:"<<std::endl;

  if(is_n)   {std::cout<<"n = "<<n<<" unitless"<<std::endl;   }
  if(is_l)   {std::cout<<"l = "<<l<<" unitless"<<std::endl;   }
  if(is_m)   {std::cout<<"m = "<<m<<" unitless"<<std::endl;   }
  if(is_alpha) {std::cout<<"alpha = "<<alpha<<" Bohr^-1"<<std::endl;   }
  if(is_center){std::cout<<"center = "<<*center<<" Bohr"<<std::endl;   }


  std::cout<<std::endl;

}

double PrimitiveS::norm(){
  double res;
/*
  Meaning: if N - is a result of this function and
  S(n,l,m,alpha) = (r-center)^(n-1) * exp(-alpha * (r - center)) * Y_lm(theta,phi) is a STO
  then
  s = N * S(l,m,n,alp) - is normaized: integral(s,s) = 1.0
*/

  res = sqrt(2.0*alpha/FACTORIAL(2*n))*FAST_POW(2.0*alpha,n);  // this is normalization factor


  return res;
}


