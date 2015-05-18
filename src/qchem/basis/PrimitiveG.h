#ifndef PrimitiveG_H
#define PrimitiveG_H


#include "Mathematics.h"

//============ Primitive Gaussian functions class ====================

class PrimitiveG{
  // f = (x-center.x)^{x_exp} * (y-center.y)^{y_exp} * (z-center.z)^(z_exp) * exp(-alpha * (r - center)^2)
public:
  // Members
  int x_exp;          int is_x_exp;
  int y_exp;          int is_y_exp;
  int z_exp;          int is_z_exp;
  double G_alpha;     int is_G_alpha;
  VECTOR* G_center;   int is_G_center;
 
  // Function value
  double G_value;     int is_G_value;

    
  //----------------------------------------------
  // Constructor
  PrimitiveG(){
    // Simple S-type function
    x_exp = 0.0;         is_x_exp = 1;
    y_exp = 0.0;         is_y_exp = 1;
    z_exp = 0.0;         is_z_exp = 1;
    G_alpha = 1.0;       is_G_alpha = 1;
                         is_G_center = 0;
                         is_G_value = 0;
  }
  // Arbitrary L function constructor
  PrimitiveG(int l,int m,int n,double alp,VECTOR& center){ 
    x_exp = l;           is_x_exp = 1;
    y_exp = m;           is_y_exp = 1;
    z_exp = n;           is_z_exp = 1;
    G_alpha = alp;       is_G_alpha = 1;
    G_center = &center;  is_G_center = 1;
                         is_G_value = 0;
  }

  // Copy constructor
  PrimitiveG(const PrimitiveG&);  

  void init(int l,int m,int n,double alp,VECTOR& center){
    x_exp = l;           is_x_exp = 1;
    y_exp = m;           is_y_exp = 1;
    z_exp = n;           is_z_exp = 1;
    G_alpha = alp;       is_G_alpha = 1;
    G_center =&center;   is_G_center = 1;
                         is_G_value = 0;
  }
  void init(int l,int m,int n,double alp,VECTOR* center){
    x_exp = l;           is_x_exp = 1;
    y_exp = m;           is_y_exp = 1;
    z_exp = n;           is_z_exp = 1;
    G_alpha = alp;       is_G_alpha = 1;
    G_center = center;   is_G_center = 1;
                         is_G_value = 0;
  }


  double Evaluate(VECTOR&);
  double norm();

  PrimitiveG& operator=(const PrimitiveG&);
  void show_info();
      

};


#endif // PrimitiveG_H
