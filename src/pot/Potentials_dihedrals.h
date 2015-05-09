#ifndef POTENTIALS_DIHEDRALS_H
#define POTENTIALS_DIHEDRALS_H

#include "../mmath/libmmath.h"
using namespace libmmath;


namespace libpot{




//------------------ Dihedral/Torsion potentials ------------------------------

double Dihedral_General(VECTOR& ri,VECTOR& rj,VECTOR& rk,VECTOR& rl, /*Inputs*/
                        VECTOR& fi,VECTOR& fj,VECTOR& fk,VECTOR& fl, /*Outputs*/
                        double Vphi,double phi0,int n,int opt        /*Parameters*/
                        );
double Dihedral_Fourier(VECTOR& ri,VECTOR& rj,VECTOR& rk,VECTOR& rl,    /*Inputs*/
                        VECTOR& fi,VECTOR& fj,VECTOR& fk,VECTOR& fl,    /*Outputs*/
                        double Vphi1,double Vphi2,double Vphi3,int opt  /*Parameters*/
                        );

}// namespace libpot

#endif //POTENTIALS_DIHEDRALS_H
