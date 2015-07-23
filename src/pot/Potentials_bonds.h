#ifndef POTENTIALS_BONDS_H
#define POTENTIALS_BONDS_H

#include "../mmath/libmmath.h"
using namespace libmmath;
using namespace libmmath::liblinalg;


namespace libpot{

//------------------ Bond potentials ------------------------------

double Bond_Harmonic(VECTOR& ri,VECTOR& rj,  /*Inputs*/
                     VECTOR& fi,VECTOR& fj,  /*Outputs*/
                      double K, double r0);  /*Parameters*/
boost::python::list Bond_Harmonic(VECTOR ri,VECTOR rj,    /*Inputs*/
                                  double K, double r0);   /*Parameters*/


double Bond_Quartic(VECTOR& ri,VECTOR& rj,  /*Inputs*/
                    VECTOR& fi,VECTOR& fj,  /*Outputs*/
                    double K, double r0);   /*Parameters*/
boost::python::list Bond_Quartic(VECTOR ri,VECTOR rj,    /*Inputs*/
                                 double K, double r0);   /*Parameters*/


double Bond_Morse(VECTOR& ri,VECTOR& rj,            /*Inputs*/
                  VECTOR& fi,VECTOR& fj,            /*Outputs*/
                  double D, double r0,double alp);  /*Parameters*/
boost::python::list Bond_Morse(VECTOR ri,VECTOR rj,              /*Inputs*/
                               double D, double r0,double alp);  /*Parameters*/



}// namespace libpot

#endif // POTENTIALS_BONDS_H
