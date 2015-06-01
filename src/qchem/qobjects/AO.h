#ifndef AO_H
#define AO_H

#include "../../mmath/libmmath.h"
using namespace libmmath;

#include "PrimitiveG.h"

namespace libqchem{
namespace libqobjects{


class AO{

// Objects of this class are atomic orbitals
// AO = linear combination of primitive gaussians (contraction)

public:

  // Members
  std::string element;              int is_element;
  std::string ao_shell;             int is_ao_shell;      // e.g. 1s, 2p, 3p, 3d, etc.
  std::string ao_shell_type;        int is_ao_shell_type; // e.g. s, p, d,  etc.
  std::string ao_name;              int is_ao_name;       // e.g. 1s, 2px, 3dxy,  etc.  int at_indx;

  // Angular momentum is defined for AOs, but each AO may have many Gaussians with different
  // alpha constants
  int x_exp;                        int is_x_exp;
  int y_exp;                        int is_y_exp;
  int z_exp;                        int is_z_exp;
  int expansion_size;               int is_expansion_size;
  vector<PrimitiveG> primitives;    int is_primitives;
  vector<double> coefficients;      int is_coefficients;  // these coefficients correspond to normalized primitive Gaussians



  //----------------- Function members --------------------
  AO();  ///   c-tor
  AO(const AO&); /// cc-tor
  AO& operator=(const AO&);

  // Construction and print  
  void init();
  void clear();
  void add_primitive(double c,PrimitiveG g);
  void show_info();

  // Computations
  double compute(VECTOR&);
  double norm2();
  double norm1();
  double normalization_factor();
  void normalize();


  // Transformations
  void shift_position(const VECTOR&);

};

}// namespace libqobjects
}// namespace libqchem



#endif // AO_H
