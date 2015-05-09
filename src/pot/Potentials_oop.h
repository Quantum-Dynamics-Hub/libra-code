#ifndef POTENTIALS_OOP_H
#define POTENTIALS_OOP_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libpot{

//------------------ Out-of-plane (OOP) potentials -------------------------------

double OOP_Fourier(VECTOR& r1,VECTOR& r2,VECTOR& r3,VECTOR& r4,        /*Inputs*/
                   VECTOR& f1,VECTOR& f2,VECTOR& f3,VECTOR& f4,        /*Outputs*/
                   double Kijkl,double C0,double C1,double C2,int opt  /*Parameters*/
                   );

double OOP_Wilson(VECTOR& r1,VECTOR& r2,VECTOR& r3,VECTOR& r4,   /*Inputs*/
                  VECTOR& f1,VECTOR& f2,VECTOR& f3,VECTOR& f4,   /*Outputs*/
                  double Kijkl,double xi_0                       /*Parameters*/
                   );

double OOP_Harmonic(VECTOR& r0,VECTOR& r1,VECTOR& r2,VECTOR& r3,  /*Inputs*/
                    VECTOR& f0,VECTOR& f1,VECTOR& f2,VECTOR& f3,  /*Outputs*/
                    double Kijkl                                  /*Parameters*/
                   );


}// namespace libpot

#endif //POTENTIALS_OOP_H
