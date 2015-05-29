#ifndef ERI_H
#define ERI_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{

double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize, 
   int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD,
    vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
);

double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize, 
   int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD
);

boost::python::list electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize, int is_derivs
);

double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd,
   int is_normalize
);

double electron_repulsion_integral
(
   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
   int nxc,int nyc, int nzc, double alp_c, VECTOR& Rc,
   int nxd,int nyd, int nzd, double alp_d, VECTOR& Rd
);






}// namespace libmolint
}// namespace libqchem

#endif //ERI_H
