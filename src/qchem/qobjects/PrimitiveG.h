#ifndef PrimitiveG_H
#define PrimitiveG_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libqobjects{


//============ Primitive Gaussian functions class ====================

class PrimitiveG{
  // f = (x-X)^{x_exp} * (y-Y)^{y_exp} * (z-Z)^(z_exp) * exp(-alpha * (r - R)^2)
  // R = (X, Y, Z)^T

public:
  // Members data
  int x_exp;          int is_x_exp;
  int y_exp;          int is_y_exp;
  int z_exp;          int is_z_exp;
  double alpha;       int is_alpha;
  VECTOR R;           int is_R; 
  // Function value
  double value;       int is_value;


  // Getters-setters
  int get_x_exp();
  int get_y_exp();
  int get_z_exp();
  double get_alpha();
  VECTOR get_R();
  double get_value();
                                   
  void set_x_exp(int _x);
  void set_y_exp(int _y);
  void set_z_exp(int _z);
  void set_alpha(double _alp);
  void set_R(VECTOR& _R);



  //----------------------------------------------
  // Member functions
  void init(int l,int m,int n,double alp,VECTOR& center);  /// general initialization/setup

  PrimitiveG::PrimitiveG();   /// Default c-tor
  PrimitiveG::PrimitiveG(int l,int m,int n,double alp,VECTOR& center); /// general c-tor
  PrimitiveG(const PrimitiveG&);  /// Copy c-tor
  PrimitiveG& operator=(const PrimitiveG&);  /// assignment operator


  // Computations and print
  double compute(VECTOR&);
  double norm2(); 
  double norm1(); 
  double normalization_factor(); 
  void show_info();

  // Transformations
  void shift_position(const VECTOR&);
     
};

// Molints functions exported for use with Gaussians objects:
double gaussian_overlap
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
);

double gaussian_overlap
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB
);

boost::python::list gaussian_overlap
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs
);

double gaussian_overlap( PrimitiveG& GA, PrimitiveG& GB,int is_normalize);

double gaussian_overlap( PrimitiveG& GA, PrimitiveG& GB);





}// namespace libqobjects
}// namespace libqchem



#endif // PrimitiveG_H
