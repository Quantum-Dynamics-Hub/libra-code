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

#include "Potentials_mb_elec.h"

/// liblibra namespace
namespace liblibra{


namespace libpot{




double Elec_Ewald3D(vector<VECTOR*>& r, vector<double>& q, MATRIX3x3& box,    /* Inputs */ 
                    vector<VECTOR*>& f, MATRIX3x3& at_stress,  /* Outputs*/
                    vector< std::pair<int, int> >& exclusions, 
                    vector< VECTOR* >& excl_vector1,
                    vector< VECTOR* >& excl_vector2,
                    vector<double>& excl_scales,
                    int rec_deg,int pbc_deg, double etha, double R_on, double R_off    /* Parameters */
                   ){

// This is gonna be a cleaner and simpler version of the previous one
// this time, we are going for simplicity and modularity


//********************* double ELECTROSTATIC *****************************
//* 3D Ewald summation given by:                                         *
//*  Karasawa, N.; Goddard III,W. A. "Acceleration of Convergence for    *
//*  Lattice Sums" J.Phys.Chem. 1989, 93,7320-7327  (/Theory/Ewald1.pdf) *
//* Exclusions and formula clarification are given by:                   *
//* 1) Procacci, P.; Marchi, M. "Taming the Ewald sum in molecular       *
//*  dynamics simulations of solvated proteins via a multiple time step  *
//*  algorithm" J.Chem.Phys. 1996, 104, 3003-3012 (/Theory/Ewald2.pdf)   *
//* 2) Komeiji, Y. "Ewald summation and multiple time step methods for   *
//*  molecular dynamics simulation of biological molecules"              *
//*  J. Mol. Struct. (Theochem) 2000, 530, 237-243 (/Theory/Ewald3.pdf)  *
//*                                                                      *
//*  u = S1 + S2 + S3                                                    *
//*  S1 = (1/2etha)*SUMM { Qij *erfc(a)/a }                              *
//*                L,i,j                                                 *
//*                                                                      *
//*  S2 = (2pi/Omega)*SUMM'{S(h)S(-h)h^-2 exp(-b*b)}                     *
//*                     h                                                *
//*                                                                      *
//*  S3 = -(1/sqrt(pi)*etha)*SUMM {Qii}                                  *
//*                           i                                          *
//*  a = |ri - rj -RL|/etha                                              *
//*  b = 1/2*h*etha                                                      *
//************************************************************************

  int Lx,Ly,Lz;
  int gi,gj;
  double a,b,hmod,fact;
  double Qii,Qij,r2;
  double etha_13;
  double SW;
  VECTOR dSW,f_mod;
  double E1,E2,E3,energy;
  VECTOR f_mod1,f_mod2,dij,rij;
  VECTOR h,h1,h2,h3,t,tv1,tv2,tv3,av,tv,g1,g2,g3;
  MATRIX3x3 tp;
  MATRIX3x3 I;  I.identity();
  double omega,erfc1,s1mod;
  double const1,const2;
  int xshift,yshift,zshift;
  int i,j;
  int sz = r.size();

  etha_13 = (1.0/(etha*etha*etha));
  
  // Reciprocal vectors
  box.get_vectors(tv1,tv2,tv3);
  box.inverse().T().get_vectors(g1,g2,g3);
  t.cross(tv2,tv3);    h1 = 2.0*M_PI*t/(tv1*t);
  t.cross(tv3,tv1);    h2 = 2.0*M_PI*t/(tv2*t);
  t.cross(tv1,tv2);    h3 = 2.0*M_PI*t/(tv3*t);

  omega = box.Determinant(); 
  const1 = (2.0*M_PI*electric/omega);
  const2 = 2.0/sqrt(M_PI);

  //------------------ Initialize forces and stress -----------------
  energy = 0.0;
//  for(i=0;i<sz;i++){ *f[i] = 0.0; }
//  at_stress = 0.0;


  // =============== S1 (exclude i==j pairs for central cell)=========
  // This is original version
  E1 = 0.0;
  for(int na=-pbc_deg;na<=pbc_deg;na++){
    for(int nb=-pbc_deg;nb<=pbc_deg;nb++){
      for(int nc=-pbc_deg;nc<=pbc_deg;nc++){
  
        //  Note: It is very important that the translation vectors are symmetric, leading to:
        //  summ (over tv)  is equivalent to summ (over -tv)
 
        tv = (na*tv1 + nb*tv2 + nc*tv3);

        for(int i=0;i<sz;i++){
          for(int j=0;j<sz;j++){
             
            rij = *r[i] - *r[j] - tv;
            r2 = rij.length2();
 
            if(na==0 && nb==0 && nc==0 && i==j){ } // skip this
            else{

            SW = 1.0; dSW = 0.0;
            VECTOR rj = *r[j] + tv;
            SWITCH(*r[i],rj,R_on,R_off,SW,dSW); 

            if(SW>0.0){
              a = sqrt(r2)/etha;
              Qij = q[i]*q[j]*electric;
              erfc1 = ERFC(a)/a;
              E1 += SW*erfc1*Qij;

              s1mod = etha_13*SW*(Qij/(a*a))*(erfc1 + const2*exp(-a*a) ); 
              f_mod = 0.5*(s1mod*rij - (erfc1/etha)*Qij*dSW);
              *f[i] += f_mod;
              *f[j] -= f_mod;

              tp.tensor_product(rij,f_mod);   at_stress += tp;

            }// if SW>0.0
            }//else
              
          }// for j
        }// for i

      }// for nc
    }// for nb
  }// for na


  //======== intra (excluded pairs are included in reciprocal part - exclude it) ======
  double E1_excl = 0.0;

  int nexcl = exclusions.size();

  for(int e=0;e<nexcl;e++){
    i = exclusions[e].first; 
    j = exclusions[e].second;
    VECTOR rij; rij = (*r[i] + *excl_vector1[e]) - (*r[j] + *excl_vector2[e]);
    double scl = excl_scales[e];
    

    if(i!=j){

      r2 = rij.length2();
      double d1 = sqrt(r2);
      double en = electric*q[i]*q[j]/d1;

      VECTOR fi = (en/r2)*(rij/d1);

      // Note that here we are subtracting both energy and forces
      // Scaling works differently: we already have included the full interaction in
      // the reciprocal space, so we are going to exclude the entire one..
      // The scaling returns part of it. E.g. if scaling is 0.2 we are going to add overall -0.8
      E1_excl += (scl - 1.0)*en;
      *f[i] += (scl - 1.0)*fi;
      *f[j] -= (scl - 1.0)*fi;

      tp.tensor_product(rij,fi);   at_stress -= (scl - 1.0)*tp;

    }
  }// for e - all exclusions



  // ================== S2 =======================
  // Reciprocal summ contribution
  E2 = 0.0;
  for(Lx=-rec_deg;Lx<=rec_deg;Lx++){
    for(Ly=-rec_deg;Ly<=rec_deg;Ly++){
      for(Lz=-rec_deg;Lz<=rec_deg;Lz++){
  
        //  Note: It is very important that the translation vectors are symmetric, leading to:
        //  summ (over h)  is equivalent to summ (over -h)
  
        h = Lx*h1 + Ly*h2 + Lz*h3;
        hmod = h.length();
        b = 0.5*hmod*etha;

        if((Lx==0)&&(Ly==0)&&(Lz==0)){ }//Exclude this case
        else{
          fact = (exp(-b*b)/(hmod*hmod));
          double sum1 = 0.0;
          double sum2 = 0.0;
          VECTOR vsum1, vsum2; 
          vsum1 = 0.0;
          vsum2 = 0.0;
          for(i=0;i<sz;i++){
            sum1 += q[i]*cos(h* *r[i]);
            sum2 += q[i]*sin(h* *r[i]);
          }// for i
        double sum3 = fact*(sum1*sum1 + sum2*sum2);
        E2 += sum3;

        //----------- Forces ---------------------
        for(i=0;i<sz;i++){
          f_mod2 = h*fact*(sin(h* *r[i])*sum1 - cos(h* *r[i])*sum2);
          f_mod  = 2.0*const1*q[i]*f_mod2;
          *f[i] += f_mod;
        }

        // Tensors
        // This is an all-atomic pressure tensor due to Ewald sum
          tp.tensor_product(h,h); tp = sum3*const1*(I - 2.0*((1.0 + b*b)/(h*h))*tp );  
          at_stress += tp;

        }// else
      }// for Lz
    }// for Ly
  }// for Lx

  energy += const1*E2;
  //==============================================

  //=============== Additive constant to energy (self-interactions) ===========
  // it does not contribute to forces
  E3 = 0.0;
  for(i=0;i<sz;i++){  E3 += electric*(q[i]*q[i]);  }
  energy -= (0.5*const2/etha)*E3;


  return energy;
}


double Elec_Ewald3D(vector<VECTOR>& r, vector<double>& q, MATRIX3x3& box, double epsilon,  /* Inputs */ 
                    vector<VECTOR>& f, MATRIX3x3& at_stress,  /* Outputs*/
                    int rec_deg,int pbc_deg, double etha, double R_on, double R_off    /* Parameters */
                   ){
/**
  Simplest, Python-friendly Ewald sum function - no exclusions

  This function takes coordinates in a.u. (Bohrs) and returns the energy in a.u. (Hatree)
*/


//********************* double ELECTROSTATIC *****************************
//* 3D Ewald summation given by:                                         *
//*  Karasawa, N.; Goddard III,W. A. "Acceleration of Convergence for    *
//*  Lattice Sums" J.Phys.Chem. 1989, 93,7320-7327  (/Theory/Ewald1.pdf) *
//* Exclusions and formula clarification are given by:                   *
//* 1) Procacci, P.; Marchi, M. "Taming the Ewald sum in molecular       *
//*  dynamics simulations of solvated proteins via a multiple time step  *
//*  algorithm" J.Chem.Phys. 1996, 104, 3003-3012 (/Theory/Ewald2.pdf)   *
//* 2) Komeiji, Y. "Ewald summation and multiple time step methods for   *
//*  molecular dynamics simulation of biological molecules"              *
//*  J. Mol. Struct. (Theochem) 2000, 530, 237-243 (/Theory/Ewald3.pdf)  *
//*                                                                      *
//*  u = S1 + S2 + S3                                                    *
//*  S1 = (1/2etha)*SUMM { Qij *erfc(a)/a }                              *
//*                L,i,j                                                 *
//*                                                                      *
//*  S2 = (2pi/Omega)*SUMM'{S(h)S(-h)h^-2 exp(-b*b)}                     *
//*                     h                                                *
//*                                                                      *
//*  S3 = -(1/sqrt(pi)*etha)*SUMM {Qii}                                  *
//*                           i                                          *
//*  a = |ri - rj -RL|/etha                                              *
//*  b = 1/2*h*etha                                                      *
//************************************************************************

  MATRIX3x3 tp;
  MATRIX3x3 I;  I.identity();
  int i,j;
  int sz = r.size();


  // Reciprocal vectors
  VECTOR tv1,tv2,tv3; // unit cell vectors
  VECTOR g1,g2,g3, t; // auxiliary vectors
  VECTOR h1,h2,h3;    // reciprocal space vectors

  box.get_vectors(tv1,tv2,tv3);
  box.inverse().T().get_vectors(g1,g2,g3);
  t.cross(tv2,tv3);    h1 = 2.0*M_PI*t/(tv1*t);
  t.cross(tv3,tv1);    h2 = 2.0*M_PI*t/(tv2*t);
  t.cross(tv1,tv2);    h3 = 2.0*M_PI*t/(tv3*t);

  double omega = box.Determinant(); 
  double const1 = (2.0*M_PI/omega)/epsilon;
  double const2 = (2.0/sqrt(M_PI))/epsilon;
  double const3 = (0.5/etha);

  //------------------ Initialize forces and stress -----------------
  double energy = 0.0;
  for(i=0;i<sz;i++){ f[i] = 0.0; }
  at_stress = 0.0;


  // =============== S1 (exclude i==j pairs for central cell)=========
  // This is original version
  double E1 = 0.0;
  double SW;  // switching function
  VECTOR dSW; // and its derivative
  VECTOR tv;  // lattice translation vectors
  VECTOR derfc1_dri, dEdri, rj;
  VECTOR rij; // distance between particles, including lattice translation


  for(int na=-pbc_deg;na<=pbc_deg;na++){
    for(int nb=-pbc_deg;nb<=pbc_deg;nb++){
      for(int nc=-pbc_deg;nc<=pbc_deg;nc++){
  
        //  Note: It is very important that the translation vectors are symmetric, leading to:
        //  summ (over tv)  is equivalent to summ (over -tv)
 
        tv = (na*tv1 + nb*tv2 + nc*tv3);

        for(int i=0;i<sz;i++){
          for(int j=0;j<sz;j++){
             
            rij = r[i] - r[j] - tv;
            double rijt = rij.length2();
 
            if(na==0 && nb==0 && nc==0 && i==j){ } // skip this
            else{

              SW = 1.0; dSW = 0.0;
              rj = r[j]; rj += tv;
              SWITCH(r[i],rj,R_on,R_off,SW,dSW); 

              if(SW>0.0){
                /**
                       E = (0.5/etha) * q[i]*q[j]/epsilon * ERFC(a)/a * SW
                */

                double a = sqrt(rijt)/etha;
                double Qij = q[i]*q[j]/epsilon;
                double erfc1 = ERFC(a)/a;
                E1 += Qij*erfc1*SW;


                /**
                       dE/dr_i = (0.5/etha) * (q[i]*q[j]/epsilon) * (ERFC(a)/a) * dSW/dr_i +
                                     SW *  q[i]*q[j]/epsilon * (0.5/etha) * d (ERFC(a)/a)/da * da/dr_i 

                       d(ERFC(a)/a)/da =   ((-2*exp(-a**2) / sqrt(pi) ) * a - ERFC(a)  )/a**2 = 

                       = -(  (2*exp(-a**2)/ sqrt(pi) )  + (ERFC(a)/a)  ) /a 
 
                       da/dr_i  = (rij / sqrt(rijt) ) / etha = rij / (sqrt(rijt) * etha) = rij / (a * etha**2)

                      So, the last term together:
 
                 - SW *  q[i]*q[j]/epsilon * (0.5/etha**3) * rij * (  (2*exp(-a**2)/ sqrt(pi) )  + (ERFC(a)/a)  ) /a**2 

                Note that: a*a*etha*etha = rijt 

                */

               VECTOR derfc1_dri = -((epsilon*const2*exp(-a*a) + erfc1)/rijt ) * rij;
               VECTOR dEdri  = const3 * Qij * ( derfc1_dri * SW + erfc1 * dSW);

                f[i] -= dEdri;
                f[j] += dEdri;

                tp.tensor_product(rij,dEdri);   at_stress -= tp;

              }// if SW>0.0
            }//else
              
          }// for j
        }// for i

      }// for nc
    }// for nb
  }// for na

  energy += const3*E1;


  // ================== S2 =======================
  // Reciprocal sum contribution
  double E2 = 0.0;
  double sum1 = 0.0;
  double sum2 = 0.0;
  double fact;
  VECTOR h;  // reciprocal space lattice translation vectors
  VECTOR f_mod; // force

  for(int Lx=-rec_deg;Lx<=rec_deg;Lx++){
    for(int Ly=-rec_deg;Ly<=rec_deg;Ly++){
      for(int Lz=-rec_deg;Lz<=rec_deg;Lz++){
  
        //  Note: It is very important that the translation vectors are symmetric, leading to:
        //  summ (over h)  is equivalent to summ (over -h)
  
        h = Lx*h1 + Ly*h2 + Lz*h3;
        double hmod = h.length();
        double b = 0.5*hmod*etha;

        if((Lx==0)&&(Ly==0)&&(Lz==0)){ }//Exclude this case
        else{
          fact = (exp(-b*b)/(hmod*hmod));

          sum1 = 0.0;   sum2 = 0.0; 
          for(i=0;i<sz;i++){

              sum1 += q[i]*cos(h*r[i]);
              sum2 += q[i]*sin(h*r[i]);

          }// for i

          double sum3 = fact*(sum1*sum1 + sum2*sum2);
          E2 += sum3;

          //----------- Forces ---------------------
          for(i=0;i<sz;i++){

              f_mod = h;
              f_mod  *= const1 * 2.0*fact*q[i]*(sin(h* r[i])*sum1 - cos(h* r[i])*sum2); // -dsum3/dr_i 
              f[i] += f_mod;  
          }

          // Tensors
          // This is an all-atomic pressure tensor due to Ewald sum
          tp.tensor_product(h,h); tp = sum3*const1*(I - 2.0*((1.0 + b*b)/(h*h))*tp );  
          at_stress += tp;

        }// else
      }// for Lz
    }// for Ly
  }// for Lx

  energy += const1*E2;
  //==============================================

  //=============== Additive constant to energy (self-interactions) ===========
  // it does not contribute to forces
  double E3 = 0.0;
  for(i=0;i<sz;i++){  E3 += (q[i]*q[i]);  }
  energy -= (0.5*const2/etha)*E3;


  return energy;

}




double Elec_Ewald3D(VECTOR* r,                                               /* Inputs */ 
                    VECTOR* g,
                    VECTOR* m,
                    VECTOR* f, 
                    MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
                    int sz,double* q,
                    int nexcl, int* excl1, int* excl2, double* scale,
                    MATRIX3x3* box,int rec_deg,int pbc_deg,
                    double etha,int is_cutoff, double R_on, double R_off,int& time,
                    vector< vector<triple> >& images, vector<triple>& central_translation,
                    double* dr2,double dT_2, int& update_displ2     /* Parameters */
                   ){
//********************* double ELECTROSTATIC *****************************
//* 3D Ewald summation given by:                                         *
//*  Karasawa, N.; Goddard III,W. A. "Acceleration of Convergence for    *
//*  Lattice Sums" J.Phys.Chem. 1989, 93,7320-7327  (/Theory/Ewald1.pdf) *
//* Exclusions and formula clarification are given by:                   *
//* 1) Procacci, P.; Marchi, M. "Taming the Ewald sum in molecular       *
//*  dynamics simulations of solvated proteins via a multiple time step  *
//*  algorithm" J.Chem.Phys. 1996, 104, 3003-3012 (/Theory/Ewald2.pdf)   *
//* 2) Komeiji, Y. "Ewald summation and multiple time step methods for   *
//*  molecular dynamics simulation of biological molecules"              *
//*  J. Mol. Struct. (Theochem) 2000, 530, 237-243 (/Theory/Ewald3.pdf)  *
//*                                                                      *
//*  u = S1 + S2 + S3                                                    *
//*  S1 = (1/2etha)*SUMM { Qij *erfc(a)/a }                              *
//*                L,i,j                                                 *
//*                                                                      *
//*  S2 = (2pi/Omega)*SUMM'{S(h)S(-h)h^-2 exp(-b*b)}                     *
//*                     h                                                *
//*                                                                      *
//*  S3 = -(1/sqrt(pi)*etha)*SUMM {Qii}                                  *
//*                           i                                          *
//*  a = |ri - rj -RL|/etha                                              *
//*  b = 1/2*h*etha                                                      *
//************************************************************************

  int Lx,Ly,Lz;
  int gi,gj;
  double a,b,hmod,fact;
  double Qii,Qij,r2;
  double etha_13;
  double SW;
  VECTOR dSW,f_mod;
  double E1,E2,E3,energy;
  VECTOR f_mod1,f_mod2,dij,rij,gij,mij;
  VECTOR h,h1,h2,h3,t,tv1,tv2,tv3,av,tv,g1,g2,g3;
  MATRIX3x3 tp;
  MATRIX3x3 I;  I.identity();
  double omega,erfc1,s1mod;
  double const1,const2;
  int xshift,yshift,zshift;
  int i,j;

  etha_13 = (1.0/(etha*etha*etha));
  
  // Reciprocal vectors
  box->get_vectors(tv1,tv2,tv3);
  box->inverse().T().get_vectors(g1,g2,g3);
  t.cross(tv2,tv3);    h1 = 2.0*M_PI*t/(tv1*t);
  t.cross(tv3,tv1);    h2 = 2.0*M_PI*t/(tv2*t);
  t.cross(tv1,tv2);    h3 = 2.0*M_PI*t/(tv3*t);

  omega = box->Determinant(); 
  const1 = (2.0*M_PI*electric/omega);
  const2 = 2.0/sqrt(M_PI);

  //------------------ Initialize forces and stress -----------------
  energy = 0.0;
  for(i=0;i<sz;i++){ f[i] = 0.0; }
  at_stress = 0.0;
  fr_stress = 0.0;
  ml_stress = 0.0;

/*
  // =============== S1 (exclude i==j pairs for central cell)=========
  // This is original version
  E1 = 0.0;
  for(int na=-pbc_deg;na<=pbc_deg;na++){
    for(int nb=-pbc_deg;nb<=pbc_deg;nb++){
      for(int nc=-pbc_deg;nc<=pbc_deg;nc++){
  
  //  Note: It is very important that the translation vectors are symmetric, leading to:
  //  summ (over tv)  is equivalent to summ (over -tv)
 
        tv = (na*tv1 + nb*tv2 + nc*tv3);

        for(int i=0;i<sz;i++){
          for(int j=0;j<sz;j++){
             
            rij = r[i] - r[j] - tv;
            gij = g[i] - g[j] - tv;
            mij = m[i] - m[j] - tv;
            r2 = rij.length2();
 
            if(na==0 && nb==0 && nc==0 && i==j){ } // skip this
            else{

            SW = 1.0; dSW = 0.0;
            VECTOR rj = r[j]+tv;
            if(is_cutoff){ SWITCH(r[i],rj,R_on,R_off,SW,dSW); }
            if(SW>0.0){
              a = sqrt(r2)/etha;
              //erfc1 = (boost::math::erfc<double>(a))/a;
              Qij = q[i]*q[j]*electric;
              erfc1 = ERFC(a)/a;
              E1 += SW*erfc1*Qij;

              s1mod = etha_13*SW*(Qij/(a*a))*(erfc1 + const2*exp(-a*a) ); 
              f_mod = 0.5*(s1mod*rij - (erfc1/etha)*Qij*dSW);
              f[i] += f_mod;
              f[j] -= f_mod;

              tp.tensor_product(rij,f_mod);   at_stress += tp;
              tp.tensor_product(gij,f_mod);   fr_stress += tp;
              tp.tensor_product(mij,f_mod);   ml_stress += tp;

            }// if SW>0.0
            }//else
              
          }// for j
        }// for i

      }// for nc
    }// for nb
  }// for na
*/

  double R_skin2 = 4.0;
  Cell cl(tv1,tv2,tv3,R_off+2.0);
  triple central_translation_ij;
  vector<triple> images_ij;

  update_displ2 = 0;

  // =============== S1 (exclude i==j pairs for central cell)=========
  // This is new version - to accomodate for new treatment of PBC
  int count = 0;
  int time_ij = time%5;
  E1 = 0.0;
/*
  if(time_ij==0){
    if(images.size()>0){ images.clear(); }
    if(central_translation.size()>0) { central_translation.clear(); }
  }
*/
 int is_new = 0;
 if(images.size()==0){ is_new = 1; }
 double scl = 2.0; // use 1 if full i,j loops are used: (i in [1,sz]) and (j in [1,sz])
                   // use 2 if half i,j loops are used: (i in [1,sz]) and (j in [i,sz])

  for(i=0;i<sz;i++){
    for(j=i;j<sz;j++){

      rij = r[i] - r[j];


// No Verlet list
//    cl.calculate(rij,images_ij,central_translation_ij);

// With Verlet list
/*
      if(time_ij==0){
        cl.calculate(rij,images_ij,central_translation_ij);  
        images.push_back(images_ij);
        central_translation.push_back(central_translation_ij);
      }else{
        images_ij = images[count];
        central_translation_ij = central_translation[count];
      }      
*/

// Even more efficient version  
     int is_update = 0;
     if(is_new) { is_update = 1; }
     else{
       if(images[count].size()==0){ is_update = 1; }
       else{       
         // Calculate new central translation:
         int nxshift = floor(rij*g1+0.5);
         int nyshift = floor(rij*g2+0.5);
         int nzshift = floor(rij*g3+0.5);
         if((nxshift!=central_translation[count].n1)||
            (nyshift!=central_translation[count].n2)||
            (nzshift!=central_translation[count].n3)
           ){ is_update = 1; }
         else{          
           double displ2  = dr2[i] + dr2[j] + dT_2 ;
           if(displ2>R_skin2){ is_update = 1;  }
         }// else
       }// else
     }// else

     is_update = 1; // this is no Verlet list option - for checking

      if(is_update && is_new){
        // Update list
        cl.calculate(rij,images_ij,central_translation_ij);
        images.push_back(images_ij);
        central_translation.push_back(central_translation_ij);
        update_displ2 = 1;
      }
      else if(is_update && !is_new){
        cl.calculate(rij,images_ij,central_translation_ij);
        images[count] = images_ij;
        central_translation[count]  = central_translation_ij;
        update_displ2 = 1;
      }
      else{
        // Use existing list
        images_ij = images[count];
        central_translation_ij = central_translation[count];        
      }     


      count++;
      int n_images = images_ij.size();
      xshift = central_translation_ij.n1;
      yshift = central_translation_ij.n2;
      zshift = central_translation_ij.n3;

//      xshift = floor(rij*g1+0.5);
//      yshift = floor(rij*g2+0.5);
//      zshift = floor(rij*g3+0.5);


      for(int im=0;im<n_images;im++){
//      for(int na=(-pbc_deg+xshift);na<=(pbc_deg+xshift);na++){
//        for(int nb=(-pbc_deg+yshift);nb<=(pbc_deg+yshift);nb++){
//          for(int nc=(-pbc_deg+zshift);nc<=(pbc_deg+zshift);nc++){
         int na = images_ij[im].n1;
         int nb = images_ij[im].n2;
         int nc = images_ij[im].n3;

  //  Note: It is very important that the translation vectors are symmetric, leading to:
  //  summ (over tv)  is equivalent to summ (over -tv)

            tv = (na*tv1 + nb*tv2 + nc*tv3);

            if((na==xshift) && (nb==yshift) && (nc==zshift) && (i==j)){ 
            // skip this - gamma point
            }
            else{

              rij = r[i] - r[j] - tv;
              gij = g[i] - g[j] - tv;
              mij = m[i] - m[j] - tv;
              r2 = rij.length2();

              SW = 1.0; dSW = 0.0;
              VECTOR rj = r[j]+tv;
              if(is_cutoff){ SWITCH(r[i],rj,R_on,R_off,SW,dSW); }
              if(SW>0.0){
                a = sqrt(r2)/etha;
                //erfc1 = (boost::math::erfc<double>(a))/a;
                Qij = q[i]*q[j]*electric;
                erfc1 = ERFC(a)/a;
                E1 += SW*erfc1*Qij;

                s1mod = etha_13*SW*(Qij/(a*a))*(erfc1 + const2*exp(-a*a) );
                f_mod = 0.5*(s1mod*rij - (erfc1/etha)*Qij*dSW);
                f_mod *= scl;
                f[i] += f_mod;
                f[j] -= f_mod;

                tp.tensor_product(rij,f_mod);   at_stress += tp;
                tp.tensor_product(gij,f_mod);   fr_stress += tp;
                tp.tensor_product(mij,f_mod);   ml_stress += tp;

              }// if SW>0.0
            }//else  - not gamma point

//          }// for nc
//        }// for nb
//      }// for na
      }// for im
    }// for j
  }// for i

  time_ij++;
  time = time_ij;

  energy += (0.5/etha)*E1*scl;

  //cout<<"(0.5/etha)*E1*scl = "<<(0.5/etha)*E1*scl<<endl;

  //======== intra (excluded pairs are included in reciprocal part - exclude it) ======
  E1 = 0.0;
  E2 = 0.0;
  for(int e=0;e<nexcl;e++){
    i = excl1[e]; 
    j = excl2[e];

//    cout<<"e = "<<e<<" i = "<<i<<" j = "<<j<<" a = "<<a<<endl;
    if(i!=j){

      rij = r[i] - r[j];
      xshift = floor(rij*g1+0.5);
      yshift = floor(rij*g2+0.5);
      zshift = floor(rij*g3+0.5);
      tv = (xshift*tv1 + yshift*tv2 + zshift*tv3);

      Qij = electric*q[i]*q[j];
      rij = r[i] - r[j] - tv;
      gij = g[i] - g[j] - tv;
      mij = m[i] - m[j] - tv;

      r2 = rij.length2();
      a = sqrt(r2)/etha;

      //cout<<"e = "<<e<<" i = "<<i<<" j = "<<j<<" a = "<<a<<endl;


      //================== direct space contributions (which should be excluded) =====================
      SW = 1.0; dSW = 0.0;
      VECTOR rj; rj = r[j]+tv;
      if(is_cutoff){ SWITCH(r[i],rj,R_on,R_off,SW,dSW); }

      erfc1 = ERFC(a)/a;
      E2 += SW*erfc1*Qij;
      s1mod = etha_13*SW*(Qij/(a*a))*(erfc1 + const2*exp(-a*a) );
      f_mod = (s1mod*rij - (erfc1/etha)*Qij*dSW);
      f[i] -= f_mod;
      f[j] += f_mod;

      tp.tensor_product(rij,f_mod);   at_stress -= tp;
      tp.tensor_product(gij,f_mod);   fr_stress -= tp;
      tp.tensor_product(mij,f_mod);   ml_stress -= tp;


//    s1mod = (Qij/(a*a))*(erfc1 + const2*exp(-a*a) );
//    f_mod = (etha_13*s1mod*rij);
//    f[i] -= f_mod;
//    f[j] += f_mod;


      //================== reciprocal space contributions =========================
      erfc1 = ERF(a)/a;
      E1 += erfc1*Qij;

      s1mod = (Qij/(a*a))*(erfc1 - const2*exp(-a*a) );
      f_mod = (etha_13*s1mod*rij);

      f[i] -= f_mod;
      f[j] += f_mod;

      tp.tensor_product(rij,f_mod);   at_stress -= tp;
      tp.tensor_product(gij,f_mod);   fr_stress -= tp; 
      tp.tensor_product(mij,f_mod);   ml_stress -= tp;

      //=========================================================================

    }// i!=j

  }
  energy -= ((E1+E2)/etha);
//  cout<<"E1 = "<<E1<<endl;
//  cout<<"E2 = "<<E2<<endl;
//  cout<<"etha = "<<etha<<endl;
//  cout<<"-((E1+E2)/etha) = "<<-((E1+E2)/etha)<<endl;

  // ================== S2 =======================
  // Reciprocal summ contribution
  E2 = 0.0;
  for(Lx=-rec_deg;Lx<=rec_deg;Lx++){
    for(Ly=-rec_deg;Ly<=rec_deg;Ly++){
      for(Lz=-rec_deg;Lz<=rec_deg;Lz++){
  
  //  Note: It is very important that the translation vectors are symmetric, leading to:
  //  summ (over h)  is equivalent to summ (over -h)
  

        h = Lx*h1 + Ly*h2 + Lz*h3;
        hmod = h.length();
        b = 0.5*hmod*etha;

        if((Lx==0)&&(Ly==0)&&(Lz==0)){ }//Exclude this case
        else{
          fact = (exp(-b*b)/(hmod*hmod));
          double sum1 = 0.0;
          double sum2 = 0.0;
          VECTOR vsum1, vsum2; 
          vsum1 = 0.0;
          vsum2 = 0.0;
          for(i=0;i<sz;i++){
            sum1 += q[i]*cos(h*r[i]);
            sum2 += q[i]*sin(h*r[i]);
          }// for i
        double sum3 = fact*(sum1*sum1 + sum2*sum2);
        E2 += sum3;

        //----------- Forces ---------------------
        for(i=0;i<sz;i++){
          f_mod2 = h*fact*(sin(h*r[i])*sum1 - cos(h*r[i])*sum2);
          f_mod  = 2.0*const1*q[i]*f_mod2;
          f[i] += f_mod;
          // Here we need correction due to rigid-body constraints (reaction)
          // Which is not included in Ewald sum
          tp.tensor_product((r[i]-g[i]),-f_mod); fr_stress += tp;
          tp.tensor_product((r[i]-m[i]),-f_mod); ml_stress += tp;
        }

        // Tensors
        // This is a all-atomic pressure tensor due to Ewald sum
          tp.tensor_product(h,h); tp = sum3*const1*(I - 2.0*((1.0 + b*b)/(h*h))*tp );  
          at_stress += tp;
          fr_stress += tp;
          ml_stress += tp;

        }// else
      }// for Lz
    }// for Ly
  }// for Lx

  //cout<<"const1*E2 = "<<const1*E2<<endl;

  energy += const1*E2;
  //==============================================

  //=============== Additive constant to energy (self-interactions) ===========
  E3 = 0.0;
  for(i=0;i<sz;i++){  E3 += electric*(q[i]*q[i]);  }
  energy -= (0.5*const2/etha)*E3;
  //cout<<"-(0.5*const2/etha)*E3 = "<<-(0.5*const2/etha)*E3<<endl;

  //cout<<"energy = "<<energy<<endl;
  return energy;
}


}// namespace libpot
}// liblibra

