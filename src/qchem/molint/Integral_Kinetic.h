#ifndef INTEGRAL_KINETIC_H
#define INTEGRAL_KINETIC_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{


// 1D Gaussian kinetic energy
double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize, int is_derivs, double& dI_dxa,double& dI_dxb,
                        vector<double*>& aux,int n_aux );
double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize, int is_derivs, double& dI_dxa,double& dI_dxb
                       );
boost::python::list kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                                     int is_normalize, int is_derivs );
double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize );
double kinetic_integral(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb );




// 3D Gaussians kinetic energy
double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
                        int is_derivs,
                        VECTOR& dIdA, VECTOR& dIdB,
                        vector<double*>& auxd,int n_aux
                       );
double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, int is_derivs,
                        VECTOR& dIdA, VECTOR& dIdB
                       );
boost::python::list kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                     int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                     int is_normalize, int is_derivs
                                    );
double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize
                       );
double kinetic_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                       );









}// namespace libmolint
}// namespace libqchem


#endif // INTEGRAL_KINETIC_H
