#ifndef NUCLEAR_H
#define NUCLEAR_H

#include "../../mmath/libmmath.h"
using namespace libmmath;

namespace libdyn{
namespace libnuclear{


class Nuclear{
  // This clas represents a point in a N-dimensional phase space, N = 6*Nnucl
  // plus some properties are associated with each projection (coordinate)

  public:

  //------------- Data members ---------------------
  int n_nucl;                     // number of nuclear DOFs

  // Atomic properties
  vector<double> mass;            // (generalized) masses of particles

  // Dynamic variables
  vector<double> q;               // nuclear DOFs
  vector<double> p;               // conjugate momenta
  vector<double> f;               // forces


  //------------- Constructors ----------------------
  Nuclear();
  Nuclear(int _n_nucl);
  
};


} // libnuclear
} // libdyn


#endif // NUCLEAR_H
