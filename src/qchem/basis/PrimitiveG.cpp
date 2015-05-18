#include "PrimitiveG.h"

//----------------------- Members of PrimitiveG class -----------------------
// Copy constructor
PrimitiveG::PrimitiveG(const PrimitiveG& g){

  // Default parameters
  x_exp = 0.0;         is_x_exp = 1;
  y_exp = 0.0;         is_y_exp = 1;
  z_exp = 0.0;         is_z_exp = 1;
  G_alpha = 1.0;       is_G_alpha = 1;
                       is_G_center = 0;
                       is_G_value = 0;

  if(g.is_x_exp)   {  x_exp = g.x_exp; is_x_exp = 1; }
  if(g.is_y_exp)   {  y_exp = g.y_exp; is_y_exp = 1; }
  if(g.is_z_exp)   {  z_exp = g.z_exp; is_z_exp = 1; }
  if(g.is_G_alpha) {  G_alpha = g.G_alpha; is_G_alpha = 1; }
  if(g.is_G_center){
    G_center = new VECTOR;  
    *G_center = *g.G_center; is_G_center = 1;
  }
  if(g.is_G_value) {  G_value = g.G_value; is_G_value = 1; }

}

PrimitiveG& PrimitiveG::operator=(const PrimitiveG& g){

  // Default parameters
  x_exp = 0.0;         is_x_exp = 1;
  y_exp = 0.0;         is_y_exp = 1;
  z_exp = 0.0;         is_z_exp = 1;
  G_alpha = 1.0;       is_G_alpha = 1;
                       is_G_center = 0;
                       is_G_value = 0;
 
  if(g.is_x_exp)   {  x_exp = g.x_exp; is_x_exp = 1; }
  if(g.is_y_exp)   {  y_exp = g.y_exp; is_y_exp = 1; }
  if(g.is_z_exp)   {  z_exp = g.z_exp; is_z_exp = 1; }
  if(g.is_G_alpha) {  G_alpha = g.G_alpha; is_G_alpha = 1; }
  if(g.is_G_center){  *G_center = *g.G_center; is_G_center = 1; }
  if(g.is_G_value) {  G_value = g.G_value; is_G_value = 1;}

  return *this;

}

void PrimitiveG::show_info(){

  std::cout<<"PrimitiveG properties:"<<std::endl;

  if(is_x_exp)   {std::cout<<"x_exp = "<<x_exp<<" unitless"<<std::endl;   }
  if(is_y_exp)   {std::cout<<"y_exp = "<<y_exp<<" unitless"<<std::endl;   }
  if(is_z_exp)   {std::cout<<"z_exp = "<<z_exp<<" unitless"<<std::endl;   }
  if(is_G_alpha) {std::cout<<"G_alpha = "<<G_alpha<<" Bohr^-1"<<std::endl;   }
  if(is_G_center){std::cout<<"G_center = "<<*G_center<<" Bohr"<<std::endl;   }
  if(is_G_value) {std::cout<<"G_value = "<<G_value<<" unitless"<<std::endl;   }

  std::cout<<std::endl;

}


double PrimitiveG::Evaluate(VECTOR& pos){

  VECTOR r; r = pos - *G_center;
  double r2 = r.length2();
  G_value = FAST_POW(r.x,x_exp)*FAST_POW(r.y,y_exp)*FAST_POW(r.z,z_exp)*exp(-G_alpha * r2);
  is_G_value = 1;

  return G_value;
}

double PrimitiveG::norm(){
  double res;
/*
  double nom = pow((8.0*G_alpha),(x_exp+y_exp+z_exp))*FACTORIAL(x_exp)*FACTORIAL(y_exp)*FACTORIAL(z_exp);
  double denom = FACTORIAL(2*x_exp)*FACTORIAL(2*y_exp)*FACTORIAL(2*z_exp);
  res = pow((2.0*G_alpha/M_PI),0.75)*sqrt(nom/denom);

  Meaning: if N - is a result of this function and
  G(l,m,n,alp) = x^l * y^m * z^n * exp(-alp*r^2)  a primitive Gaussian
  then
  g = N * G(l,m,n,alp) - is normaized: integral(g,g) = 1.0
*/

//  OVERLAP_INTEGRAL(primitives[i],primitives[j],1)


  double x = FAST_POW(G_alpha,(2*(x_exp+y_exp+z_exp)+3));
         x = pow(x,0.25);
         x = x * FAST_POW(2,(x_exp+y_exp+z_exp));
  double y = sqrt(DFACTORIAL(2*x_exp-1)*DFACTORIAL(2*y_exp-1)*DFACTORIAL(2*z_exp-1));


  res = pow((2.0/M_PI),0.75) *(x/y);

  return res;
}


