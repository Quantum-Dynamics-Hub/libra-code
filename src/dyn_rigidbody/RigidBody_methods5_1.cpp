/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file RigidBody_methods5_1.cpp
  \brief The file implements the exact propagation method for free RB, based on Euler elliptic integrals

*/

#include "RigidBody.h"

/// liblibra namespace
namespace liblibra{


/// librigidbody namespace
namespace librigidbody{


void RigidBody::initialize_exact_rb(){
/**
  \brief Initialize auxiliary variables for exact free RB integration

  This function pre-computes some auxiliary quantities which are used
  later in exact solution to the free rigid-body problem in terms
  of Jacobi elliptic functions as described by:
  1) Hernandez de la Pena, L.; van Zon, R.; Schofield, J.; Opps, S. B.
  "Discontinuous molecular dynamics for semiflexible and rigid bodies"
  J. Chem. Phys. 2007, 126, 074105
  2) van Zon, R.; Schofield, J. "Numerical implementation of the exact
  dynamics of free rigid bodies" J. Comp. Phys. 2007, 225, 145-164
*/

  // Here we need to check if the necessary quantities has been calculated
  // and calculate them if necessary
//  I    = rb_I_e;
  invI = rb_invI_e;

  top_l = rb_l_e;
  top_w = rb_invI_e * top_l;

  // Calculate the rotational energy, total angular momentum and its square
  double Er = 0.5*top_l * top_w;
  double L2 = top_l.length2();
  double L  = sqrt(L2);
  MATRIX3x3 I_tmp;
  double Ix,Iy,Iz,sgnx,sgny,sgnz,eta,ksi;
  int status;
  int n =1;
  double dA2;
  double q2n, ksin,base;
  double r0,i0;

  if(L2>0.0){

  //================== Special case ================================
  if((top_l.x!=0.0 && top_l.y==0.0 && top_l.z==0.0 ) ||
     (top_l.x==0.0 && top_l.y!=0.0 && top_l.z==0.0 ) ||
     (top_l.x==0.0 && top_l.y==0.0 && top_l.z!=0.0 )
    ){ ;; }
  //================= General case =================================
  else{

  // Choose the ordering
  status = 0;

  //----------------- Spherical top -----------------------------------
  if((rb_A==rb_B) && (rb_B==rb_C)){   status = 1;  permutindx = 0;   }

  //------------------ Symmetric top ----------------------------------
  else if(rb_A==rb_B || rb_A==rb_C || rb_B==rb_C){
    for(int permut=0;permut<6;permut++){
      MATRIX3x3 invU; invU = U[permut].inverse();
      I_tmp = U[permut]*rb_I_e*U[permut];
      Ix = I_tmp.xx;
      Iy = I_tmp.yy;
      Iz = I_tmp.zz;
      // Find such permutation that Ix = Iy
      if( Ix == Iy ){ invpermutindx = permutindx = permut;
                      if(permutindx==3){ invpermutindx = 4; }
                      else if(permutindx==4){ invpermutindx = 3; }
                      status = 1;  }
    }// for all permutations

    if(!status){
      cout<<"In symmetric top case. The necessary permutation can not be found"<<endl;
      exit(1);
    }else{
    // Update all variables to suit for the ordering just found
//      MATRIX invU; U[permutindx].Inverse(invU);
      top_w = U[permutindx] * top_w;
      top_l = U[permutindx] * top_l;
      Iint  = U[permutindx] * rb_I_e * U[permutindx];
      invI  = U[permutindx] * rb_invI_e * U[permutindx];
      I1    = Iint.xx;
      I2    = Iint.yy;
      I3    = Iint.zz;
    }
    // Calculate precession frequency
      wp = (1.0-(I3/I1))*top_w.z;
  } // if symmetric top

  //----------------- Asymmetric top ---------------------------------
  else{

    for(int permut=0;permut<6;permut++){
      int invpermut = permut;
      if(permut==3){ invpermut = 4; }
      else if(permut==4){ invpermut = 3; }


      I_tmp = U[permut]*rb_I_e*U[invpermut]; //cout<<"permut = "<<permut<<"I_tmp = "<<I_tmp<<endl;
      Ix = I_tmp.xx;
      Iy = I_tmp.yy;
      Iz = I_tmp.zz;

      // Find Jacobi ordering:
      if( ((Ix>=Iy)&&(Iy>=Iz) && (Er>=0.5*L2/Iy))  ||
          ((Ix<=Iy)&&(Iy<=Iz) && (Er<0.5*L2/Iy))
          ){ invpermutindx = permutindx = permut;
             if(permutindx==3){ invpermutindx = 4; }
             else if(permutindx==4){ invpermutindx = 3; }
             status = 1; break; }
    }// for all permutations

    if(!status){
      cout<<"In asymmetric top case. The necessary permutation can not be found"<<endl;
      exit(1);
    }else{
    top_w = U[permutindx] * top_w;
    top_l = U[permutindx] * top_l;
    Iint  = U[permutindx] * rb_I_e* U[invpermutindx];
    I1    = Iint.xx;
    I2    = Iint.yy;
    I3    = Iint.zz;
    }

//------------------------------------------------------------------

    // Calculate the amplitudes of the angular velocity
    sgnx = (top_w.x==0.0)?1.0:SIGN(top_w.x);
    sgnz = (top_w.z==0.0)?1.0:SIGN(top_w.z);
    sgny = (top_w.y==0.0)?1.0:SIGN(top_w.y);
    sgny = -1.0*sgnx;


    double L12 = L2 - 2.0*Er*I3;
    double L23 = L2 - 2.0*Er*I1;
    double I13 = I1 - I3;
    double I23 = I2 - I3;
    double omega1,omega2,omega3;
    omega1 = top_w.x; omega2 = top_w.y; omega3 = top_w.z;
    double omega1m = copysign(sqrt(L12/I1/I13),omega1);
    double omega2m = -copysign(sqrt(omega2*omega2+I1*I13*omega1*omega1/I2/I23),omega1);
    double omega3m = copysign(sqrt(-L23/I3/I13),omega3);
    double omegap  = I23*copysign(sqrt(L23/(-I23)/I1/I2/I3), omega3);


    top_wm.x =  sgnx * sqrt((L2 - 2.0*I3*Er)/(I1*(I1-I3)));
    top_wm.y =  sgny * sqrt((L2 - 2.0*I3*Er)/(I2*(I2-I3)));
    top_wm.z =  sgnz * sqrt((L2 - 2.0*I1*Er)/(I3*(I3-I1)));

    // Precession frequency
    wp   = SIGN(I2-I3)*sgnz*sqrt((L2-2.0*I1*Er)*(I3-I2)/(I1*I2*I3));
    wp = omegap;

    // Elliptic parameter
    m    = (L2-2.0*I3*Er)*(I1-I2)/((L2-2.0*I1*Er)*(I3-I2));

    // Elliptic integral
    Ellint(m,(top_w.y/top_wm.y),1e-12,K,eps);  




    double f;
    Ellint((1.0-m),(I3*top_wm.z/L),1e-12,Kcompl,f); 
    //cout<<"Combined version:\n";
    //cout<<"f = "<<f<<" Kcompl = "<<Kcompl<<endl;

    // ...instead of:
    //Kcompl = Km((1.0-m),1e-12);  //boost::math::ellint_1<double>(sqrt(1.0-m));
    q      = exp(-M_PI*Kcompl/K);              // nome
    //eta    = SIGN(top_w.z)*Kcompl - boost::math::ellint_1<double,double>(sqrt(1.0-m),asin(I3*top_wm.z/L)); //Ellipe(asin(I3*top_wm.z/L),sqrt(1-m),IntN);
    //f = boost::math::ellint_1<double,double>(sqrt(1.0-m),asin(I3*top_wm.z/L));
    //cout<<"Original version version:\n";
    //cout<<"f = "<<f<<" Kcompl = "<<Kcompl<<endl;

    eta = SIGN(top_w.z)*Kcompl - f;
    ksi    = exp(M_PI*eta/K);

//------------------------------------------------------------------------

    A2     = L/I1 + M_PI*wp*(1.0+ksi)/(2.0*K*(ksi-1.0));

    base = -(M_PI*wp/K);
    n = 1;
    q2n  = q*q;
    ksin = ksi;
    do{
        dA2 = (q2n/(1.0-q2n))*(ksin-(1.0/ksin));
        A2 += base*dA2;
        q2n *= (q*q);
        ksin*=  ksi;
        n++;
    } while(fabs(dA2)>MACHPREC);

//------------------------------------------------------

    NT = (int)(log(MACHPREC)/log(q)+0.5);

    r0 = i0 = 0.0;
    if(cr!=NULL){ delete [] cr; }
    if(ci!=NULL){ delete [] ci; }
    cr = new double[NT];
    ci = new double[NT];

    base = 2.0*sqrt(sqrt(1.0));// q - according to formula, 1 - according to code
    double q_n = 1.0;
    double q_2n = 1.0;
    double q_n_2 = 1.0;

    for(n=0;n<NT;n++){
        cr[n] =  base*q_n*q_n_2*cosh(0.5*(2.0*n+1.0)*M_PI*eta/K);
        ci[n] = -base*q_n*q_n_2*sinh(0.5*(2.0*n+1.0)*M_PI*eta/K);

        r0   += cr[n]*sin(0.5*(2.0*n+1.0)*M_PI*eps/K);
        i0   += ci[n]*cos(0.5*(2.0*n+1.0)*M_PI*eps/K);
        q_n  *= (-q);
        q_2n *= (q*q);
        q_n_2 *= (q_2n / q);

    }

    // Alternative implementation
    A1 = atan2(i0,r0);

//--------------------------------------------------------------

  }// another general case - non-spherical top
  //------------------------------------------------------------------
  }// general case
  //==================================================================
  }// if L2>0.0

}

void RigidBody::propagate_exact_rb(double dt){
/**
  \brief Exact propagation of the free RB dynamics

  \param[in] dt The propagation duration

  This function provides the exact solution to the free rigid-body
  problem in terms of Jacobi elliptic functions as described by:
  1) Hernandez de la Pena, L.; van Zon, R.; Schofield, J.; Opps, S. B.
  "Discontinuous molecular dynamics for semiflexible and rigid bodies"
  J. Chem. Phys. 2007, 126, 074105
  2) van Zon, R.; Schofield, J. "Numerical implementation of the exact
  dynamics of free rigid bodies" J. Comp. Phys. 2007, 225, 145-164

*/
  // Pre-compute auxiliary parameters
  initialize_exact_rb();

  double t;
  double t0 = -(eps/wp);

  double u, am, sn, cn, dn;
  VECTOR nz(0.0,0.0,1.0);
  VECTOR omegat,dir,L;
  MATRIX3x3 T1,T2,P,Rz,W;
  double psi,cos_psi,sin_psi;
  double Rev,Imv,Mod,Km_,arg,C,S;
  double Lpn, Lmodn,L2;
  double s1,s2,s3,c1,c2,c3;

  double Lp; Lp = sqrt(top_l.x*top_l.x + top_l.y*top_l.y);
  double Lmod = top_l.length();

  if(Lmod>0.0){

  //=================== Special case =========================
  if((top_l.x!=0.0 && top_l.y==0.0 && top_l.z==0.0 ) ||
     (top_l.x==0.0 && top_l.y!=0.0 && top_l.z==0.0 ) ||
     (top_l.x==0.0 && top_l.y==0.0 && top_l.z!=0.0 )
    ){
      exit(1);
      T2.Rotation(-t*top_w);
      rb_A_I_to_e = T2 * rb_A_I_to_e;
  }
  //================== General case ==========================
  else{

  //----------------- Spherical top case ---------------------
  if((I1==I2) && (I2==I3)){
    exit(1);
    T2.Rotation(-t*top_w);
    rb_A_I_to_e = T2 * rb_A_I_to_e;
  }
  //----------------- Symmetric top case ---------------------
  else if(I1==I2) {
    exit(1);
    dir = top_l/Lmod;
    psi = t*(top_l.length()/I1);
    cos_psi = cos(psi);
    sin_psi = sin(psi);
    T2.Rotation(-t*(top_l/I1));

    dir = nz;
    psi = t * wp;
    cos_psi = cos(psi);
    sin_psi = sin(psi);
    T1.Rotation(-t*wp*nz);

    Rz.xx = cos_psi;     Rz.xy = sin_psi;    Rz.xz = 0.0;
    Rz.yx =-sin_psi;     Rz.yy = cos_psi;    Rz.yz = 0.0;
    Rz.zx = 0.0;         Rz.zy = 0.0;        Rz.zz = 1.0;

    Rz = U[permutindx] * Rz * U[permutindx]; // This accounts for intermediate permutations
                                                 // of vector l_e
    P = T1 * T2;
    P = U[permutindx] * P * U[permutindx]; // This accounts for intermediate permutations
                                               // of matrix A_I_to_e
    rb_A_I_to_e = P * rb_A_I_to_e;
    rb_l_e      = Rz * rb_l_e;


  }
  //------------------- Asymmetric top case ------------------------
  else{

    t = dt;
    u = wp*t + eps;
    Jacobi_Elliptic(u,m,tol,am,sn,cn,dn);
    // Calculate angular velocity
    omegat.x = top_wm.x * cn;
    omegat.y = top_wm.y * sn;
    omegat.z = top_wm.z * dn;

    // Calculate the rotation angle psi (via sin and cos)
    Rev = Imv = 0.0;
    t = U[permutindx].Determinant()*dt;
    u = wp*t + eps;
    Km_ = Km(m,1e-12); //boost::math::ellint_1<double>(sqrt(m)); //Km(m,tol);
    for(int n=0;n<NT;n++){
      arg = 0.5*(2.0*n+1.0)*M_PI*(u)/Km_;
      Rev += cr[n]*sin(arg);
      Imv += ci[n]*cos(arg);
    }
    Mod = sqrt(Rev*Rev + Imv*Imv);

    C = cos(A1 + A2*t);
    S = sin(A1 + A2*t);
    cos_psi = (C*Rev + S*Imv)/Mod;
    sin_psi = (S*Rev - C*Imv)/Mod;

//-----------------------------------------------------------------------

    // Calculate orientation matrix propagation and final angular momenta
    L = Iint * omegat;

    // Calculate transformation matrix T1 (with ordered variables)
    Lpn = sqrt(L.x*L.x + L.y*L.y);
    Lmodn = L.length();
    L2 = Lmodn * Lmodn;

    MATRIX3x3 T1t0,T1t;
    s1 = (Lp/Lmod);    c1 = (top_l.z/Lmod);
    s2 = (top_l.y/Lp); c2 = (top_l.x/Lp);

    /*****************************************
            c1  0  -s1           c2  s2  0
    T1t =   0   1   0      *     -s2 c2  0
            s1  0   c1           0   0   1
    *****************************************/
    T1t0.xx = c1*c2;  T1t0.xy = c1*s2;  T1t0.xz = -s1;
    T1t0.yx = -s2;    T1t0.yy = c2;     T1t0.yz = 0.0;
    T1t0.zx = s1*c2;  T1t0.zy = s1*s2;  T1t0.zz = c1;

    s1 = (Lpn/Lmodn); c1 = (L.z/Lmodn);
    s2 = (L.y/Lpn);   c2 = (L.x/Lpn);
    /*****************************************
            c1  0  -s1           c2  s2  0
    T1t =   0   1   0      *     -s2 c2  0
            s1  0   c1           0   0   1
    *****************************************/
    T1t.xx = c1*c2;  T1t.xy = c1*s2;  T1t.xz = -s1;
    T1t.yx = -s2;    T1t.yy = c2;     T1t.yz = 0.0;
    T1t.zx = s1*c2;  T1t.zy = s1*s2;  T1t.zz = c1;

    MATRIX3x3 T2p;
    T2p.xx = cos_psi;   T2p.xy = sin_psi;  T2p.xz = 0.0;
    T2p.yx = -sin_psi;  T2p.yy = cos_psi;  T2p.yz = 0.0;
    T2p.zx = 0.0;       T2p.zy = 0.0;      T2p.zz = 1.0;

    P = T1t.T() * T2p * T1t0;
    P = U[invpermutindx] * P * U[permutindx]; // This accounts for intermediate permutations
                                              // of matrix A_I_to_e
    rb_A_I_to_e = P * rb_A_I_to_e;

    rb_l_e = U[invpermutindx] * L;
    rb_w_e = U[invpermutindx] * omegat;

   }// non-spherical top case
   //-------------------------------------------------------------------
   }// general case
   //====================================================================
   }// if Lmod>0.0

   // Update dependent variables
   set_angular_momentum(rb_l_e);
   set_orientation(rb_A_I_to_e);

}


}// namespace librigidbody
}// liblibra
