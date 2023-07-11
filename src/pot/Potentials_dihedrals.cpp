/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Potentials_dihedrals.h"

/// liblibra namespace
namespace liblibra{


namespace libpot{

double Dihedral_General(VECTOR& ri,VECTOR& rj,VECTOR& rk,VECTOR& rl, /*Inputs*/
                        VECTOR& fi,VECTOR& fj,VECTOR& fk,VECTOR& fl, /*Outputs*/
                        double Vphi,double phi0,int n,int opt        /*Parameters*/
                        ){
//******************** double GENERAL_TORSIONS *******************************************
//*                                                                                      *
//*  opt 0 and 1:  u_i     = Vphi*(1-cos(n*phi0)*cos(n*phi))     General1                *
//*                                                                                      *
//*  opt 2 and 3:  u_i     = Vphi*(1-cos(n*(phi-phi0)))          General2                *
//*                                                                                      *
//*  Force fields: UFF, DREIDING                                                         *
//*                                                                                      *
//****************************************************************************************
  // Derivatives are calculated using a method from:
  // Journal of Computational Chemistry, 1992, V. 13, P. 585-594
  double direction,modt,modu,cos_phi,phi,K1,K2,energy;
  energy = 0.0;
  VECTOR t,u,rp;
  if((opt==0)||(opt==2)){  direction = 1.0; } // This defines TORSION angle
  if((opt==1)||(opt==3)){  direction = -1.0;} // This defined DIHEDRAL angle
     
  VECTOR rij = ri - rj;
  VECTOR rkj = rk - rj;
  VECTOR rlk = rl - rk;
  t.cross(rij, rkj);  t = direction*t;
  u.cross(rlk, rkj);
  modt = t.length();
  modu = u.length();
  fi = fj = fk = fl = 0.0;

  if((modt!=0.0) && (modu!=0.0)){
    // Calculate angle
    cos_phi = (t*u)/(modt*modu);
    if(cos_phi>= 1.0){cos_phi=1.0;}
    else if(cos_phi<=-1.0){cos_phi=-1.0;}
    rp.cross(t,u);
    double sgn = SIGN(rkj*rp);
    sgn = (sgn==0.0)?1.0:sgn;
    phi = sgn*acos(cos_phi);

    // Energy and forces
    if((opt==0)||(opt==1)){
      K1 = Vphi*(1.0-cos(n*phi0)*cos(n*phi));
      K2 = Vphi*n*cos(n*phi0)*sin(n*phi);      // d(K1)/d(phi))
    }
    else if((opt==2)||(opt==3)){
      K1 = Vphi*(1.0-cos(n*(phi-phi0)));
      K2 = Vphi*n*sin(n*(phi-phi0));      // d(K1)/d(phi))
    }
    energy = K1;

    // Forces
    MATRIX3x3  Ti, Tj, Tk, Tl;
    MATRIX3x3  Ui, Uj, Uk, Ul;

    Ti.skew(direction*rkj);  // dt/dri
    Uj.skew(rlk);            // du/drj
    Tk.skew(-direction*rij); // dt/drk
    Ul.skew(rkj);            // du/drl
    Tj = -(Ti + Tk);
    Uk = -(Uj + Ul);

    // d(cos)/dxi = (d(cos)/dt,dt/dxi) + (d(cos)/du,du/dxi)
    // --- formula (38) divided by (-1/sin) (see formula 30) ---
    VECTOR  dfdt,dfdu;
    VECTOR tunit = t.unit();
    VECTOR uunit = u.unit();
    VECTOR rkjunit = rkj.unit();
    dfdt.cross(tunit,rkjunit);
    dfdt = (dfdt/modt);
    dfdu.cross(uunit,rkjunit);
    dfdu = -(dfdu/modu);


    //-----------------------------------------------------------
    // d(cos)/dxi = (d(cos)/dt,dt/dxi) + (d(cos)/du,du/dxi)  -see for example formula 32
    // minus added to make gradients to be forces
    // all combined - using formula 30
    fi = -K2*Ti*dfdt;
    fj = -K2*(Tj*dfdt + Uj*dfdu);
    fk = -K2*(Tk*dfdt + Uk*dfdu);
    fl = -K2*Ul*dfdu;

  }//if modt!=0.0 && modu!=0.0

  return energy;
}

double Dihedral_Fourier(VECTOR& ri,VECTOR& rj,VECTOR& rk,VECTOR& rl,    /*Inputs*/
                        VECTOR& fi,VECTOR& fj,VECTOR& fk,VECTOR& fl,    /*Outputs*/
                        double Vphi1,double Vphi2,double Vphi3,int opt  /*Parameters*/
                        ){
//******************** double FOURIER_TORSIONS *******************************************
//*                                                                                      *
//*  opt 0 and 1: u_i = V_1 *(1+cos(phi)) + V_2*(1-cos(2*phi))+ V_3*(1+cos(3*phi))       *
//*                                                                                      *
//****************************************************************************************
  // Derivatives are calculated using a method from:
  // Journal of Computational Chemistry, 1992, V. 13, P. 585-594
  double direction,modt,modu,cos_phi,phi,K1,K2,energy;
  energy = 0.0;
  VECTOR t,u,rp;
  if(opt==0){  direction = 1.0; } // This defines TORSION angle
  if(opt==1){  direction =-1.0; } // This defines DIHEDRAL angle

  VECTOR rij = ri - rj;
  VECTOR rkj = rk - rj;
  VECTOR rlk = rl - rk;
  t.cross(rij, rkj);  t = direction*t;
  u.cross(rlk, rkj);
  modt = t.length();
  modu = u.length();
  fi = fj = fk = fl = 0.0;

  if((modt!=0.0) && (modu!=0.0)){
    // Calculate angle
    cos_phi = (t*u)/(modt*modu);
    if(cos_phi>= 1.0){cos_phi=1.0;}
    else if(cos_phi<=-1.0){cos_phi=-1.0;}
    rp.cross(t,u);
    double sgn = SIGN(rkj*rp);
    sgn = (sgn==0.0)?1.0:sgn;
    phi = sgn*acos(cos_phi);

    // Energy and forces
    K1 = Vphi1*(1.0+cos(phi)) + Vphi2*(1.0-cos(2.0*phi)) + Vphi3*(1.0+cos(3.0*phi));
    K2 = -Vphi1*sin(phi) + 2.0*Vphi2*sin(2.0*phi) - 3.0*Vphi3*sin(3.0*phi);      // d(K1)/d(phi))
    
    energy = K1;

    // Forces
    MATRIX3x3  Ti, Tj, Tk, Tl;
    MATRIX3x3  Ui, Uj, Uk, Ul;

    Ti.skew(direction*rkj);  // dt/dri
    Uj.skew(rlk);            // du/drj
    Tk.skew(-direction*rij); // dt/drk
    Ul.skew(rkj);            // du/drl
    Tj = -(Ti + Tk);
    Uk = -(Uj + Ul);

    // d(cos)/dxi = (d(cos)/dt,dt/dxi) + (d(cos)/du,du/dxi)
    // --- formula (38) divided by (-1/sin) (see formula 30) ---
    VECTOR  dfdt,dfdu;
    VECTOR tunit = t.unit();
    VECTOR uunit = u.unit();
    VECTOR rkjunit = rkj.unit();
    dfdt.cross(tunit,rkjunit);
    dfdt = (dfdt/modt);
    dfdu.cross(uunit,rkjunit);
    dfdu = -(dfdu/modu);


    //-----------------------------------------------------------
    // d(cos)/dxi = (d(cos)/dt,dt/dxi) + (d(cos)/du,du/dxi)  -see for example formula 32
    // minus added to make gradients to be forces
    // all combined - using formula 30
    fi = -K2*Ti*dfdt;
    fj = -K2*(Tj*dfdt + Uj*dfdu);
    fk = -K2*(Tk*dfdt + Uk*dfdu);
    fl = -K2*Ul*dfdu;

  }//if modt!=0.0 && modu!=0.0

  return energy;
}


}// namespace libpot
}// liblibra

