#ifndef PrimitiveS_H
#define PrimitiveS_H


#include "Mathematics.h"

//============ Primitive Slater functions class ====================

class PrimitiveS{
  // f = (r-center)^(n-1) * exp(-alpha * (r - center)) * Y_lm(theta,phi)
  // Y_lm(theta,phi) = P_l^m(cos(theta)) * FHI_m(phi)  - sperical harmonics
  // P_l^m(cos(theta)) - normalized associated Legendre functions
public:
  // Members
  int n;            int is_n;
  int l;            int is_l;
  int m;            int is_m;
  double alpha;     int is_alpha;
  VECTOR* center;   int is_center;
 
    
  //----------------------------------------------
  // Constructor
  PrimitiveS(){
    // Simple S-type function
    n = 0;         is_n = 1;
    l = 0;         is_l = 1;
    m = 0;         is_m = 1;
    alpha = 1.0;   is_alpha = 1;
                   is_center = 0;

  }
  // Arbitrary L function constructor
  PrimitiveS(int _n,int _l,int _m,double _alp,VECTOR& _center){ 
    n = _n;           is_n = 1;
    l = _l;           is_l = 1;
    m = _m;           is_m = 1;
    alpha = _alp;     is_alpha = 1;
    center = &_center; is_center = 1;

  }

  // Copy constructor
  PrimitiveS(const PrimitiveS&);  

  void init(int _n,int _l,int _m,double _alp,VECTOR& _center){
    n = _n;           is_n = 1;
    l = _l;           is_l = 1;
    m = _m;           is_m = 1;
    alpha = _alp;     is_alpha = 1;
    center = &_center; is_center = 1;

  }
  void init(int _n,int _l,int _m,double _alp,VECTOR* _center){
    n = _n;           is_n = 1;
    l = _l;           is_l = 1;
    m = _m;           is_m = 1;
    alpha = _alp;     is_alpha = 1;
    center = _center; is_center = 1;

  }


//  double Evaluate(VECTOR&);
  double norm();

  PrimitiveS& operator=(const PrimitiveS&);
  void show_info();
      

};


#endif // PrimitiveS_H
