#ifndef POTENTIALS_FRAG_H
#define POTENTIALS_FRAG_H

#include "../mmath/libmmath.h"
using namespace libmmath;
using namespace libmmath::liblinalg;

namespace libpot{


//------------------------- Fragment-Fragment potentials ------------------------------------
double Gay_Berne(VECTOR& ri,VECTOR& rj,VECTOR& ui,VECTOR& uj,          /*Inputs*/
                 VECTOR& fi,VECTOR& fj,VECTOR& ti,VECTOR& tj,          /*Outputs*/
                 double di, double dj,double li,double lj,
                 double e0,double rat,double dw,double mu,double nu);  /*Parameters*/

double Girifalco12_6(VECTOR& ri,VECTOR& rj,
                     VECTOR& fi,VECTOR& fj,
                     double a,double alp,double bet
                    );

}//namespace libpot

#endif // POTENTIALS_FRAG_H
