#include "Potentials.h"

namespace libpot{

double Girifalco12_6(VECTOR& ri,VECTOR& rj,
                     VECTOR& fi,VECTOR& fj,
                     double a,double alp,double bet
                    ){
/**********************************************************
  This is the potential of interaction of 2 C60 molecules,
  with atoms interaction via LJ 12-6 potential:
  E(C-C) = B/r^12 - A/r^6
  The parameters alp and bet are related to A and B as:
  alp = A*N^2/12*(2*a)^6, bet = B*N^2/90*(2*a)^12
  N = is the number of atom on the sphere of radius a

  Reference:
  Girifalco, L. A. "Molecular Properties of C60 in Gas and 
  Solid Phases" J. Chem. Phys. 1992, 96, 858-861
**********************************************************/
  VECTOR rij = ri - rj;
  double d = rij.length();
  double s = 0.5*(d/a);

  double s2,s4,s8,s10;
  double sp,sp3,sp9;
  double sm,sm3,sm9;
  double energy;
  double f_mod;

  s2 = (1.0/s); s2 = s2 * s2;
  s4 = s2*s2;
  s8 = s4*s4;
  s10 = s2*s8;

  sp = (s+1.0); sp = (1.0/sp);
  sp3 = sp*sp*sp;
  sp9 = sp3*sp3*sp3;
  
  sm = (s-1.0); sm = (1.0/sm);
  sm3 = sm*sm*sm;
  sm9 = sm3*sm3*sm3;

  // Energy
  energy = -alp*((sm3/s) + (sp3/s) - 2.0*s4) 
           +bet*((sm9/s) + (sp9/s) - 2.0*s10);

  // Force
  f_mod = -alp*((sm3*s2) + (3.0*sm3*sm/s) + (sp3*s2) + (3.0*sp3*sp/s) - 8.0*(s4/s)) 
          +bet*((sm9*s2) + (9.0*sm9*sm/s) + (sp9*s2) + (9.0*sp9*sp/s) -20.0*(s10/s));

  f_mod = 0.5*(f_mod/a);

  fi = f_mod*(rij/d);
  fj = -fi;

  return energy;
}


}//namespace libpot

