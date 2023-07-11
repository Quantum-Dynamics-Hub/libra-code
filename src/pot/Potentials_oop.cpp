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

#include "Potentials_oop.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libpot{


double OOP_Fourier(VECTOR& r1,VECTOR& r2,VECTOR& r3,VECTOR& r4,        /*Inputs*/
                   VECTOR& f1,VECTOR& f2,VECTOR& f3,VECTOR& f4,        /*Outputs*/
                   double Kijkl,double C0,double C1,double C2,int opt  /*Parameters*/
                  ){
//******************** double FOURIER_OUT_OF_PLANE *******************************
//*                                                                              *
//*          u_i     = Kijkl*(C0 + C1*cos(gamma_ijkl) + C2*cos(2*gamma_ijkl))    *
//*                                                                              *
//*  Force fields: UFF                                                           *
//*  In input the atom 2 (r2) is central atom, connected to all other 3 atoms !  *
//*                                                                              *
//********************************************************************************
  // Derivatives are calculated using a method from: 
  // adapted to the case of out-of-plane angle
  // J. Comput. Chem. 1996, V. 17, P. 1132-1141

  VECTOR ri,rj,rk,rl,fi,fj,fk,fl;
  double energy = 0.0;
  f1 = f2 = f3 = f4 = 0.0;

  Kijkl /= 3.0; // because 3 inversions


  for(int perm=0;perm<3;perm++){
    if(perm==0)     { ri = r4; rj = r1; rk = r2; rl = r3; } // Angle LIJK
    else if(perm==1){ ri = r2; rj = r1; rk = r3; rl = r4; } // Angle JIKL
    else if(perm==2){ ri = r3; rj = r1; rk = r4; rl = r2; } // Angle KILJ

    VECTOR F = ri - rj;
    VECTOR G = rj - rk;
    VECTOR H = rj - rl;

    VECTOR A; A.cross(F,G); 
    VECTOR B; B.cross(H,G);
    //VECTOR B = rl - rj;

    double mA = A.length();
    double mB = B.length();
    double mG = G.length();

    fi = fj = fk = fl = 0.0;

    if((mA>0.0) && (mB>0.0)){

      // Calculate angle
      double cos_phi = (A*B)/(mA*mB);    
      VECTOR tmp; tmp.cross(B,A);
      double sin_phi = (tmp*G)/(mA*mB*mG);
      double phi = atan2(sin_phi,cos_phi);


      // Energy and forces
      double K1 = Kijkl*(C0 +C1*cos(phi) + C2*cos(2.0*phi));
      double K2 = -Kijkl*(C1*sin(phi) + 2.0*C2*sin(2.0*phi));      // d(K1)/d(phi)) 
      energy += K1;


      VECTOR dphi_dF = -(mG/(mA*mA))*A;
      VECTOR dphi_dG = ((G*F)/(mG*mA*mA))*A - ((H*G)/(mB*mB*mG))*B;
      VECTOR dphi_dH = (mG/(mB*mB))*B;
      //VECTOR dphi_dB; dphi_dB.cross(B,G); dphi_dB /=  (mB*mB*mG);
      

      //-----------------------------------------------------------
      fi = -K2*dphi_dF;
      fj = -K2*(-dphi_dF + dphi_dG + dphi_dH );
      fk = -K2*(-dphi_dG);
      fl = -K2*(-dphi_dH);

    }//if modt!=0.0 && modu!=0.0

  if(perm==0)     { f4 += fi; f1 += fj; f2 += fk; f3 += fl;}  // Angle LIJK
  else if(perm==1){ f2 += fi; f1 += fj; f3 += fk; f4 += fl; } // Angle JIKL
  else if(perm==2){ f3 += fi; f1 += fj; f4 += fk; f2 += fl; } // Angle KILJ
 
  }// for all permutations

  return energy;
}




double OOP_Fourier_old(VECTOR& r1,VECTOR& r2,VECTOR& r3,VECTOR& r4,        /*Inputs*/
                       VECTOR& f1,VECTOR& f2,VECTOR& f3,VECTOR& f4,        /*Outputs*/
                       double Kijkl,double C0,double C1,double C2,int opt  /*Parameters*/
                      ){
//******************** double FOURIER_OUT_OF_PLANE *******************************
//*                                                                              *
//*          u_i     = Kijkl*(C0 + C1*cos(gamma_ijkl) + C2*cos(2*gamma_ijkl))    *
//*                                                                              *
//*  Force fields: UFF                                                           *
//*  In input the atom 2 (r2) is central atom, connected to all other 3 atoms !  *
//*                                                                              *
//********************************************************************************
  // Derivatives are calculated using a method from:
  // Journal of Computational Chemistry, 1992, V. 13, P. 585-594
  double direction,modt,modu,cos_phi,phi,K1,K2,energy;
  VECTOR ri,rj,rk,rl,fi,fj,fk,fl;
  energy = 0.0;
  f1 = f2 = f3 = f4 = 0.0;
  VECTOR t,u,rp;
  if(opt==0){  direction = 1.0; } // This defines TORSION angle
  if(opt==1){  direction =-1.0; } // This defines DIHEDRAL angle

  for(int perm=0;perm<3;perm++){
    if(perm==0)     { ri = r4; rj = r1; rk = r2; rl = r3; } // Angle LIJK
    else if(perm==1){ ri = r2; rj = r1; rk = r3; rl = r4; } // Angle JIKL
    else if(perm==2){ ri = r3; rj = r1; rk = r4; rl = r2; } // Angle KILJ

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
      phi = SIGN(rkj*rp)*acos(cos_phi);

      // Energy and forces
      K1 = Kijkl*(C0 +C1*cos(phi) + C2*cos(2.0*phi));
      K2 = -Kijkl*(C1*sin(phi) + 2.0*C2*sin(2.0*phi));      // d(K1)/d(phi))
      energy += K1;

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

  if(perm==0)     { f4 += fi; f1 += fj; f2 += fk; f3 += fl;}  // Angle LIJK
  else if(perm==1){ f2 += fi; f1 += fj; f3 += fk; f4 += fl; } // Angle JIKL
  else if(perm==2){ f3 += fi; f1 += fj; f4 += fk; f2 += fl; } // Angle KILJ
 
  }// for all permutations

  return energy;
}

double OOP_Wilson(VECTOR& r1,VECTOR& r2,VECTOR& r3,VECTOR& r4,   /*Inputs*/
                  VECTOR& f1,VECTOR& f2,VECTOR& f3,VECTOR& f4,   /*Outputs*/
                  double Kijkl,double xi_0                       /*Parameters*/
                   ){
//******************** double Wilson angle (wag) *********************************
//*                                                                              *
//*          u_i     = Kijkl*(xi_ijkl-xi_0)^2                                    *
//*                                                                              *
//*  Force fields: MMFF94                                                        *
//*                                                                              *
//* In input the atom 1 (r1) is a central atom, connected to other 3 atoms       *
//********************************************************************************
// Forces are derived in
//  http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
 
  double modt,modu,cos_phi,sin_phi,cos_teta_ijk,sin_teta_ijk,tan_teta_ijk,phi,K1,K2,energy,dij,dkj,dlj;
  VECTOR ri,rj,rk,rl,fi,fj,fk,fl;
  energy = 0.0;
  f1 = f2 = f3 = f4 = 0.0;
  VECTOR t,u,rp;

  for(int perm=0;perm<3;perm++){
    if(perm==0)     { ri = r4; rj = r1; rk = r2; rl = r3;}  // Angle LIJK
    else if(perm==1){ ri = r2; rj = r1; rk = r3; rl = r4; } // Angle JIKL
    else if(perm==2){ ri = r3; rj = r1; rk = r4; rl = r2; } // Angle KILJ

    VECTOR rij = ri - rj;  dij = rij.length(); rij.normalize();
    VECTOR rkj = rk - rj;  dkj = rkj.length(); rkj.normalize();
    VECTOR rlj = rl - rj;  dlj = rlj.length(); rlj.normalize();

    // Calculate angle cosines and sines
    cos_teta_ijk = rij * rkj;
    t.cross(rij, rkj);
    sin_teta_ijk = t.length();
    tan_teta_ijk = sin_teta_ijk / cos_teta_ijk;

    u = t.unit();
    sin_phi = u * rlj;
    u.cross(u,rlj);
    cos_phi = u.length();
    phi = atan2(sin_phi,cos_phi);

    // Energy and forces
    phi = phi - xi_0;
    K1 = Kijkl*phi*phi;
    K2 = 2.0*Kijkl*phi/cos_phi;       // d(K1)/d(phi)) * d(phi)/d(sin(phi))

    energy += K1;
    // - d(sin(phi))/dri , drj, drk, drl
    t.cross(rkj,rlj);
    fi = ( t + (-rij + rkj * tan_teta_ijk * sin_phi) )/( dij * sin_teta_ijk );
    t.cross(rlj,rij);
    fk = ( t + (-rkj + rij * tan_teta_ijk * sin_phi) )/( dkj * sin_teta_ijk );
    t.cross(rij,rlj);
    fl = ( t/sin_teta_ijk - rlj * sin_phi )/dlj;
    fj = -(fi + fk + fl);
    // final forces are:
    // fi = - dU/dri = dU/dsin(phi) * (-dsin(phi)/dri)
    fi *= K2;
    fj *= K2;
    fk *= K2;
    fl *= K2;

  if(perm==0)     { f4 += fi; f1 += fj; f2 += fk; f3 += fl;}  // Angle LIJK
  else if(perm==1){ f2 += fi; f1 += fj; f3 += fk; f4 += fl; } // Angle JIKL
  else if(perm==2){ f3 += fi; f1 += fj; f4 += fk; f2 += fl; } // Angle KILJ

  }// for all permutations

  return energy;
}


double OOP_Harmonic(VECTOR& r0,VECTOR& r1,VECTOR& r2,VECTOR& r3,  /*Inputs*/
                    VECTOR& f0,VECTOR& f1,VECTOR& f2,VECTOR& f3,  /*Outputs*/
                    double Kijkl                                  /*Parameters*/
                   ){
//******************** double FOURIER_OUT_OF_PLANE *******************************
//*                                                                              *
//*          u_i     = Kijkl*h^2                                                 *
//*                                                                              *
//*  Native Force fields: Tripos5.2                                              *
//*  In input the atom 0 (r0) is central atom, connected to all other 3 atoms !  *
//*                                                                              *
//********************************************************************************
  // Derivatives done by me
  double modN,energy;
  VECTOR N,dN2dr1,dN2dr3;
  f0 = f1 = f2 = f3 = 0.0;

  VECTOR r12 = r1 - r2;
  VECTOR r32 = r3 - r2;
  VECTOR r02 = r0 - r2;

  N.cross(r12, r32); 
  modN = N.length();

  double h = (r02*N)/modN;
  double Kh = Kijkl * h;

  energy = Kh*h;


  dN2dr1.x = (2.0*r12.x - r12.y)*r32.z*r32.z + (r32.y - 2.0*r32.x)*r12.z*r32.z +
             (2.0*r12.x - r12.z)*r32.y*r32.y + (r32.z - 2.0*r32.x)*r12.y*r32.y +
             N.y*r32.y - N.z*r32.z;

  dN2dr1.y = (2.0*r12.y - r12.x)*r32.z*r32.z + (r32.x - 2.0*r32.y)*r12.z*r32.z +
             (2.0*r12.y - r12.z)*r32.x*r32.x + (r32.z - 2.0*r32.y)*r12.x*r32.x +
             N.z*r32.z - N.x*r32.x;

  dN2dr1.z = (2.0*r12.z - r12.x)*r32.y*r32.y + (r32.x - 2.0*r32.z)*r12.y*r32.y +
             (2.0*r12.z - r12.y)*r32.x*r32.x + (r32.y - 2.0*r32.z)*r12.x*r32.x +
             N.x*r32.x - N.y*r32.y;



  dN2dr3.x = (2.0*r32.x - r32.y)*r12.z*r12.z + (r12.y - 2.0*r12.x)*r32.z*r12.z +
             (2.0*r32.x - r32.z)*r12.y*r12.y + (r12.z - 2.0*r12.x)*r32.y*r12.y +
             -N.y*r12.y + N.z*r12.z;

  dN2dr3.y = (2.0*r32.y - r32.x)*r12.z*r12.z + (r12.x - 2.0*r12.y)*r32.z*r12.z +
             (2.0*r32.y - r32.z)*r12.x*r12.x + (r12.z - 2.0*r12.y)*r32.x*r12.x +
             -N.z*r12.z + N.x*r12.x;

  dN2dr3.z = (2.0*r32.z - r32.x)*r12.y*r12.y + (r12.x - 2.0*r12.z)*r32.y*r12.y +
             (2.0*r32.z - r32.y)*r12.x*r12.x + (r12.y - 2.0*r12.z)*r32.x*r12.x +
             -N.x*r12.x + N.y*r12.y;

  VECTOR p1 = (-0.5/(modN*modN*modN))*dN2dr1;
  VECTOR p3 = (-0.5/(modN*modN*modN))*dN2dr3;

  MATRIX3x3 dndr1,dndr3,tmp1,tmp3;
  tmp1.skew(r32); tmp3.skew(r12);
  dndr1.tensor_product(p1,N);  dndr1 =( dndr1 + (tmp1/modN) );
  dndr3.tensor_product(p3,N);  dndr3 =( dndr3 - (tmp3/modN) );

  f0 = 2.0*Kh*(N/modN);
  f1 = 2.0*Kh*dndr1*r02;
  f3 = 2.0*Kh*dndr3*r02;
  f2 = -(f0+f1+f3);

  return energy;
}


}//namespace libpot
}// liblibra

