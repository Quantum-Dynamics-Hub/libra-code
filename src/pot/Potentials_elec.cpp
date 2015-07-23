#include "Potentials_elec.h"

using namespace libmmath;
using namespace libmmath::liblinalg;


namespace libpot{

double Elec_Coulomb(VECTOR& ri,VECTOR& rj,     /*Inputs*/             
                    VECTOR& fi,VECTOR& fj,     /*Outputs*/
                    double qi,double qj,
                    double eps,double delta){  /*Parameters*/
//****************** double Lennard-Jones potential **************************
//*                                                                          *
//*       E = electric*qi*qj/(eps*|rij+delta|)                               *
//*                                                                          *
//****************************************************************************
  double energy,r2,r6,r12,d1,d2;
  VECTOR rij = ri - rj;
  d1 = rij.length();
  d2 = d1 + delta;
  energy = electric*(qi*qj/(eps*d2));
  fi = (energy/d2)*(rij/d1);
  fj = -fi;
  return energy;
}


}// namespace libpot


