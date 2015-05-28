#ifndef PSEUDOPOTENTIAL_H
#define PSEUDOPOTENTIAL_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{


double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB,
                   vector<double*>& auxd,int n_aux
                  );
double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB
                  );

boost::python::list pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                                int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                                int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb,
                                int is_normalize, int is_derivs
                               );

double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, const VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, const VECTOR& Rb
                  );



}// namespace libmolint
}// namespace libqchem


#endif // PSEUDOPOTENTIAL_H
