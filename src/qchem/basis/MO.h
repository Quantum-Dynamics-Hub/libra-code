#ifndef MO_H
#define MO_H

#include "AO.h"

class MO{

// Objects of this class are 1-electron molecular orbitals
// MO = linear combination of AOs

public:

  // Members
  int spin;                       int is_spin;       // may be 1 (alpha), -1(beta)

  int expansion_size;             int is_expansion_size;
  vector<AO> ao;                  int is_ao;
  vector<double> coefficients;    int is_coefficients;


  // Constructor
  MO(){
    expansion_size = 0;           is_expansion_size = 1;
                                  is_ao = 0;
                                  is_coefficients = 0;
  }

  void add_ao(double c,AO g){
    coefficients.push_back(c);  is_coefficients = 1;
    ao.push_back(g);    is_ao  = 1;
    expansion_size++;
  }

  void show_info();
  double Evaluate(VECTOR&);
  void normalize();
  double norm2();



  MO operator+(MO ob);
  void operator+=(MO ob);
//  MO operator=(MO ob);
  friend MO operator*(const double& f,  const MO& m1);  // Multiplication of MO and double;

/*
  friend int operator == (const MATRIX& m1, const MATRIX& m2);  // Are matrices equal;
  friend int operator != (const MATRIX& m1, const MATRIX& m2);  // Are matrices not equal;

  friend double operator%(MATRIX& m1, MATRIX& m2);  // Scalar product of matrices;
  friend MATRIX operator^(VECTOR& v1,VECTOR& v2);   // Tensor product of two vectors
  friend VECTOR operator*(const MATRIX& m,  const VECTOR& v);   // Multiplication of vector and matrix;
//  friend VECTOR operator*(const VECTOR& v,  const MATRIX& m);   // Multiplication of vector and matrix; 
  friend MATRIX operator*(const double& f,  const MATRIX& m1);  // Multiplication of matrix and double;
  friend MATRIX operator*(const MATRIX &m1, const double  &f);  // Multiplication of matrix and double;
*/

};


#endif // MO_H
