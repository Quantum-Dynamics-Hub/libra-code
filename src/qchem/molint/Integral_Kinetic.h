#ifndef INTEGRAL_KINETIC_H
#define INTEGRAL_KINETIC_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{


// 1D Gaussian kinetic energy
double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb );

// 3D Gaussians kinetic energy
double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                       );
double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
                        VECTOR& dIdA, VECTOR& dIdB
                       );
double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
                        VECTOR& dIdA, VECTOR& dIdB,
                        vector<double*>& auxd,int n_aux
                       );



}// namespace libmolint
}// namespace libqchem


#endif // INTEGRAL_KINETIC_H
