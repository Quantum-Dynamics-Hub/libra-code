#ifndef OVERLAPS_H
#define OVERLAPS_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{


// Basic functions
double gaussian_int(int n, double alp);
double gaussian_norm(int n,double alp);


// Basic overlaps
// 1D Gaussians
double gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb );

// 3D Gaussians overlaps
double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                       );
double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
                        VECTOR& dIdA, VECTOR& dIdB
                       );
double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
                        VECTOR& dIdA, VECTOR& dIdB,
                        vector<double*>& auxd,int n_aux
                       );



}// namespace libmolint
}// namespace libqchem


#endif // OVERLAPS_H
