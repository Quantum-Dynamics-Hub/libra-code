#ifndef OVERLAPS_H
#define OVERLAPS_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{



// Basic overlaps
// 1D Gaussians
double gaussian_overlap_ref(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb );



double gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize, int is_derivs, double& dI_dxa,double& dI_dxb,
                        vector<double*>& aux,int n_aux );


double gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize, int is_derivs, double& dI_dxa,double& dI_dxb
                       );
boost::python::list gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                                     int is_normalize, int is_derivs );



double gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize );

double gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb );







// 3D Gaussians overlaps
double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
                        int is_derivs,
                        VECTOR& dIdA, VECTOR& dIdB,
                        vector<double*>& auxd,int n_aux
                       );

double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, int is_derivs,
                        VECTOR& dIdA, VECTOR& dIdB
                       );

boost::python::list gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                     int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                     int is_normalize, int is_derivs
                                    );

double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize
                       );

double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                       );




// 1D STOs
double sto_norm(int n, double alp);
double sto_overlap(int na, int la, int ma, double alp_a, int nb, int lb, int mb, double alp_b, 
                   double R, double R_cutoff);
double sto_overlap_fast(int na, int la, int ma, double alp_a, int nb, int lb, int mb, double alp_b, 
                   double R, double R_cutoff);





}// namespace libmolint
}// namespace libqchem


#endif // OVERLAPS_H
