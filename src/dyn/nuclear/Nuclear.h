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
  int nnucl;                      // number of nuclear DOFs

  // Atomic properties
  vector<double> mass;            // (generalized) masses of particles

  // Dynamic variables
  vector<double> q;               // nuclear DOFs
  vector<double> p;               // conjugate momenta
  vector<double> f;               // forces
  vector<int> ctyp;               // coordinate type, for instance one can set 0,1,2 to be x,y,z projections of all atoms, so
                                  // in this case ctyp will be: [0, 1, 2, 0, 1, 2, 0, 1, 2. ..]
                                  // or one can set 0 to be the x coordinates of atoms in given fragment, so the ordering may look like:
                                  // [0, 0, 0, 0, 0, ... 1, 1, 1, 1, 1, ... ]


  //------------- Constructors ----------------------
  Nuclear();
  Nuclear(int _nnucl);
  Nuclear(const Nuclear& );

  void propagate_p(int i,double dt);
  void propagate_p(double dt);
  void propagate_p(double dt,vector<int>& active);

  void scale_p(int i,double scl);
  void scale_p(double scl);
  void scale_p(double scl,vector<int>& active);

  void propagate_q(int i,double dt);
  void propagate_q(double dt);
  void propagate_q(double dt,vector<int>& active);

  void scale_q(int i,double scl);
  void scale_q(double scl);
  void scale_q(double scl,vector<int>& active);


  
};



} // libnuclear
} // libdyn


#endif // NUCLEAR_H
