#ifndef PSEUDOPOTENTIAL_H
#define PSEUDOPOTENTIAL_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{


double pseudopot02(double C0, double C2, double alp, VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                  );

double pseudopot02(double C0, double C2, double alp, VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                   int is_normalize, 
                   VECTOR& dIdA, VECTOR& dIdB
                  );

double pseudopot02(double C0, double C2, double alp, VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                   int is_normalize, 
                   VECTOR& dIdA, VECTOR& dIdB,
                   vector<double*>& auxd,int n_aux
                  );


}// namespace libmolint
}// namespace libqchem


#endif // PSEUDOPOTENTIAL_H
