#ifndef MULTIPOLES_H
#define MULTIPOLES_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{

///=========== 1D =====================
double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs, double& dI_dXa, double& dI_dXb,
  vector<double*>& aux,int n_aux 
);

double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs, double& dI_dXa, double& dI_dXb
);

boost::python::list transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize, int is_derivs
);

double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb,
  int is_normalize
);

double transition_dipole_moment
( int nxa,double alp_a, double Xa, 
  int nxb,double alp_b, double Xb
);





///=========== 3D =====================
VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
);

VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
);

boost::python::list transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize,int is_derivs
);

VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
  int is_normalize
);

VECTOR transition_dipole_moment
( int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
  int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb
);





}// namespace libmolint
}// namespace libqchem


#endif // MULTIPOLES_H
