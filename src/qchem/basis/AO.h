#ifndef AO_H
#define AO_H

#include "PrimitiveG.h"
#include "PrimitiveS.h"

class AO{

// Objects of this class are atomic orbitals
// AO = linear combination of primitive gaussians (contraction)

public:

  // Members
  std::string element;              int is_element;
  std::string ao_shell;             int is_ao_shell;      // e.g. 1s, 2p, 3p, 3d, etc.
  std::string ao_shell_type;        int is_ao_shell_type; // e.g. s, p, d,  etc.
  std::string ao_name;              int is_ao_name;       // e.g. 1s, 2px, 3dxy,  etc.

  int x_exp;                        int is_x_exp;
  int y_exp;                        int is_y_exp;
  int z_exp;                        int is_z_exp;


  int expansion_size;               int is_expansion_size;
  vector<PrimitiveG> primitives;    int is_primitives;
  vector<PrimitiveS> s_primitives;  int is_s_primitives; //
  vector<double> coefficients;      int is_coefficients;

  // System-specific properties
  int at_indx;                      int is_at_indx;   // index of the atom on which this AO is centered

  // Constructor
  AO(){
                                    is_x_exp = 0;
                                    is_y_exp = 0;
                                    is_z_exp = 0;
                                    is_element = 0;
                                    is_ao_shell = 0;
                                    is_ao_shell_type = 0;
                                    is_ao_name = 0;
      expansion_size = 0;           is_expansion_size = 1;
                                    is_primitives = 0;
                                    is_s_primitives = 0;
                                    is_coefficients = 0;
                                    is_at_indx = 0;
  }
  AO(const AO&); // Copy constructor
  AO& operator=(const AO&);
  
  void clear(){
    if(primitives.size()>0) { primitives.clear(); }
    if(s_primitives.size()>0) { s_primitives.clear(); }
    if(coefficients.size()>0){ coefficients.clear(); }
                                    is_x_exp = 0;
                                    is_y_exp = 0;
                                    is_z_exp = 0;
                                    is_element = 0;
                                    is_ao_shell = 0;
                                    is_ao_name = 0;
      expansion_size = 0;           is_expansion_size = 1;
                                    is_primitives = 0;
                                    is_s_primitives = 0;
                                    is_coefficients = 0;
    
  }

  void add_primitive(double c,PrimitiveG g){   
    coefficients.push_back(c);  is_coefficients = 1;
    primitives.push_back(g);    is_primitives   = 1;
    expansion_size++;
  }

  void add_primitive(double c,PrimitiveS g){   
    coefficients.push_back(c);  is_coefficients = 1;
    s_primitives.push_back(g);    is_s_primitives   = 1;
    expansion_size++;
  }

  void show_info();
  void normalize();
  double Evaluate(VECTOR&);
  double norm();
  double norm2();

  void move(const VECTOR&);
  void set_position(const VECTOR&);

  AO operator+(AO ob);
  void operator+=(AO ob);
//  AO operator=(AO ob);
  friend AO operator*(const double& f,  const AO& m1);  // Multiplication of MO and double;



};

#endif // AO_H
