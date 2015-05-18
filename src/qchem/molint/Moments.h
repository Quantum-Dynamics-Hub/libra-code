#ifndef MOMENTS_H
#define MOMENTS_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{


// Essentially the generalized overlaps
// 1D version
double gaussian_moment(int nx, double alp, double X,
                       int nxa,double alp_a, double Xa,
                       int nxb,double alp_b, double Xb
                      );

// 3D version
double gaussian_moment(int nx, int ny,  int nz,  double alp, VECTOR& R,
                       int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                       int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                      );

double gaussian_moment(int nx, int ny,  int nz,  double alp, VECTOR& R,
                       int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                       int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                       int is_normalize, 
                       VECTOR& dIdA, VECTOR& dIdB
                      );
double gaussian_moment(int nx, int ny,  int nz,  double alp, VECTOR& R,
                       int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                       int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                       int is_normalize, 
                       VECTOR& dIdA, VECTOR& dIdB,
                       vector<double*>& auxd,int n_aux
                      );


}// namespace libmolint
}// namespace libqchem


#endif // MOMENTS_H
