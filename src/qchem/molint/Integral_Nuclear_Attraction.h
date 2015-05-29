#ifndef INTEGRAL_NUCLEAR_ATTRACTION_H
#define INTEGRAL_NUCLEAR_ATTRACTION_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{


double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                   VECTOR& Rc,int is_normalize, 
                                   int is_derivs, VECTOR& DA,VECTOR& DB, VECTOR& DC,
                                   vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv                                   
                                  );

double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                   VECTOR& Rc,int is_normalize, 
                                   int is_derivs, VECTOR& DA, VECTOR& DB, VECTOR& DC
                                  );

boost::python::list nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                                     int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                                     VECTOR& Rc, int is_normalize, int is_derivs
                                    );

double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        VECTOR& Rc, int is_normalize
                       );

double nuclear_attraction_integral(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb, VECTOR& Rc
                       );


}// namespace libmolint
}// namespace libqchem


#endif // INTEGRAL_NUCLEAR_ATTRACTION_H
