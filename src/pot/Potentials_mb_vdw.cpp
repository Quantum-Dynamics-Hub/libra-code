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
/**
  \file Potentials_mb_vdw.cpp
  This file implements the many-body potentials (with few-body potentials as sepcial case) involving
  vdW interactions. Something like the lattice sums, or just summing all the 2-body pairs without creating
  large number of auxiliary data.
*/

#include "Potentials_mb_vdw.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libpot{



double VdW_Ewald3D(vector<VECTOR>& r, vector<int>& types, int max_type, vector<double>& Bij, MATRIX3x3& box, /* Inputs */ 
                   vector<VECTOR>& f, MATRIX3x3& at_stress,  /* Outputs*/
                   int rec_deg,int pbc_deg, double etha, double R_on, double R_off    /* Parameters */                   
                   ){
/**
  Simplest, Python-friendly Ewald sum function - no exclusions

  This function takes coordinates in a.u. (Bohrs) and returns the energy in a.u. (Hatree)

  Here we use the following convention:
  Bij[a] contains the parameter for the pair indexed a = map(type_i,type_j), with the
  maping:
  a = type_i * max_type + type_j, that is the size of the vector Bij must be max_type**2
   
*/

//********************* double dispersion ********************************
//* 3D Ewald summation given by:                                         *
//*  Karasawa, N.; Goddard III,W. A. "Acceleration of Convergence for    *
//*  Lattice Sums" J.Phys.Chem. 1989, 93,7320-7327  (/Theory/Ewald1.pdf) *
//* Exclusions and formula clarification are given by:                   *
//* 1) Procacci, P.; Marchi, M. "Taming the Ewald sum in molecular       *
//*  dynamics simulations of solvated proteins via a multiple time step  *
//*  algorithm" J.Chem.Phys. 1996, 104, 3003-3012 (/Theory/Ewald2.pdf)   *
//*                                                                      *
//*  E = -1/2   sum    {  B_ij / |R_i - R_j -L|^6 }                      *
//*            i,j, L                                                    *
//*                                                                      *
//*     This is a general function, with any Bij (not necessarily with ) *
//*     the geometric mean rule                                          *
//*                                                                      *
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
  double sqrt_M_PI = sqrt(M_PI);
  double etha3 = etha * etha * etha;
  double etha6 = etha3 * etha3;
  double const1 = 0.5/etha6;
  double const2 = M_PI*sqrt_M_PI/(24.0*omega);

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

                double a = sqrt(rijt)/etha;
                double a_minus2 = 1.0/(a*a);
                double a_minus4 = a_minus2*a_minus2;
                double a_minus6 = a_minus4*a_minus2;
                double Qij = Bij[ types[i] * max_type + types[j] ]; // q[i]*q[j];

                double pref = Qij*(a_minus6 + a_minus4 + 0.5*a_minus2);
                double expa2 = exp(-a*a);
                VECTOR dpref = -(Qij*(6.0*a_minus6 + 4.0*a_minus4 + a_minus2) /(a* a * etha*etha))  *  rij;
                VECTOR dexpa2 = -(2.0*expa2/(etha*etha)) *  rij;
                
                double en = pref*expa2;
                VECTOR den_dri = (pref * dexpa2 + dpref * expa2);

                E1 -= en*SW;
                VECTOR dEdri = const1 * ( den_dri * SW + en * dSW);  // -dE1/dr[i]

                f[i] += dEdri;
                f[j] -= dEdri;

                tp.tensor_product(rij,dEdri);   at_stress += tp;

              }// if SW>0.0
            }//else
              
          }// for j
        }// for i

      }// for nc
    }// for nb
  }// for na

  energy += const1*E1;


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

        if((Lx==0)&&(Ly==0)&&(Lz==0)){    }
        else{
          fact =  ( (exp(-b*b))*(0.5/(b*b) - 1.0)/b  + sqrt_M_PI*ERFC(b) ) /(hmod*hmod*hmod);

          sum1 = 0.0;  
         
          for(i=0;i<sz;i++){
            for(j=0;j<sz;j++){

              sum1 += Bij[ types[i] * max_type + types[j] ]*cos(h*(r[i]-r[j]));

            }// for j
          }// for i

          E2 -= fact*sum1;


          //----------- Forces ---------------------
          for(i=0;i<sz;i++){
            for(j=0;j<sz;j++){

              f_mod = h;
              f_mod *= -const2 * Bij[ types[i] * max_type + types[j] ] * sin(h*(r[i]-r[j])); // -dE2/dr_i 
              f[i] += f_mod;  
              f[j] -= f_mod;  

          // Tensors
          // This is an all-atomic pressure tensor due to Ewald sum
//          tp.tensor_product(h,h); tp = sum3*const1*(I - 2.0*((1.0 + b*b)/(h*h))*tp );  
//          at_stress += tp;

            }
          }


        }// else
      }// for Lz
    }// for Ly
  }// for Lx

  energy += const2*E2;
  //==============================================

  //=============== Additive constant to energy (self-interactions) ===========
  // it does not contribute to forces
  double E3 = 0.0;
  for(i=0;i<sz;i++){  E3 += Bij[ types[i] * max_type + types[i] ];  }
  energy += (const1/6.0)*E3;

  E3 = 0.0;
  for(i=0;i<sz;i++){
    for(j=0;j<sz;j++){
      E3 += Bij[ types[i] * max_type + types[j] ];;
    }
  }
  E3 *= (4.0*const2/etha3);
  energy -= E3; 


  return energy;

}




double VdW_Ewald3D(vector<VECTOR>& r, vector<double>& q, MATRIX3x3& box, /* Inputs */ 
                   vector<VECTOR>& f, MATRIX3x3& at_stress,  /* Outputs*/
                   int rec_deg,int pbc_deg, double etha, double R_on, double R_off    /* Parameters */
                   ){
/**
  Simplest, Python-friendly Ewald sum function - no exclusions

  This function takes coordinates in a.u. (Bohrs) and returns the energy in a.u. (Hatree)

   
*/


//********************* double dispersion ********************************
//* 3D Ewald summation given by:                                         *
//*  Karasawa, N.; Goddard III,W. A. "Acceleration of Convergence for    *
//*  Lattice Sums" J.Phys.Chem. 1989, 93,7320-7327  (/Theory/Ewald1.pdf) *
//* Exclusions and formula clarification are given by:                   *
//* 1) Procacci, P.; Marchi, M. "Taming the Ewald sum in molecular       *
//*  dynamics simulations of solvated proteins via a multiple time step  *
//*  algorithm" J.Chem.Phys. 1996, 104, 3003-3012 (/Theory/Ewald2.pdf)   *
//*                                                                      *
//*  E = -1/2   sum    {  B_ij / |R_i - R_j -L|^6 }                      *
//*            i,j, L                                                    *
//*                                                                      *
//*     We assume the geometric rule: B_ij = sqrt(B_ii * B_jj)           *
//*     The meaning of the "charges" is : q[i] = sqrt(B_ii)              *
//*                                                                      *
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
  double sqrt_M_PI = sqrt(M_PI);
  double etha3 = etha * etha * etha;
  double etha6 = etha3 * etha3;
  double const1 = 0.5/etha6;
  double const2 = M_PI*sqrt_M_PI/(24.0*omega);

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

                double a = sqrt(rijt)/etha;
                double a_minus2 = 1.0/(a*a);
                double a_minus4 = a_minus2*a_minus2;
                double a_minus6 = a_minus4*a_minus2;
                double Qij = q[i]*q[j];

                double pref = Qij*(a_minus6 + a_minus4 + 0.5*a_minus2);
                double expa2 = exp(-a*a);
                VECTOR dpref = -(Qij*(6.0*a_minus6 + 4.0*a_minus4 + a_minus2) /(a* a * etha*etha))  *  rij;
                VECTOR dexpa2 = -(2.0*expa2/(etha*etha)) *  rij;
                
                double en = pref*expa2;
                VECTOR den_dri = (pref * dexpa2 + dpref * expa2);

                E1 -= en*SW;
                VECTOR dEdri  = const1 * ( den_dri * SW + en * dSW);  // -dE1/dr[i]

                f[i] += dEdri;
                f[j] -= dEdri;

                tp.tensor_product(rij,dEdri);   at_stress += tp;

              }// if SW>0.0
            }//else
              
          }// for j
        }// for i

      }// for nc
    }// for nb
  }// for na

  energy += const1*E1;


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

        if((Lx==0)&&(Ly==0)&&(Lz==0)){    }
        else{
          fact =  ( (exp(-b*b))*(0.5/(b*b) - 1.0)/b  + sqrt_M_PI*ERFC(b) ) /(hmod*hmod*hmod);

          sum1 = 0.0;   sum2 = 0.0; 
          for(i=0;i<sz;i++){

              sum1 += q[i]*cos(h*r[i]);
              sum2 += q[i]*sin(h*r[i]);

          }// for i

          E2 -= fact*(sum1*sum1 + sum2*sum2);


          //----------- Forces ---------------------
          for(i=0;i<sz;i++){

              f_mod = h;
              f_mod *= const2 * (2.0)*fact*q[i]*(sin(h* r[i])*sum1 - cos(h* r[i])*sum2); // -dsum3/dr_i 
              f[i] += f_mod;  
          }

          // Tensors
          // This is an all-atomic pressure tensor due to Ewald sum
//          tp.tensor_product(h,h); tp = sum3*const1*(I - 2.0*((1.0 + b*b)/(h*h))*tp );  
//          at_stress += tp;

        }// else
      }// for Lz
    }// for Ly
  }// for Lx

  energy += const2*E2;
  //==============================================

  //=============== Additive constant to energy (self-interactions) ===========
  // it does not contribute to forces
  double E3 = 0.0;
  for(i=0;i<sz;i++){  E3 += (q[i]*q[i]);  }
  energy += (const1/6.0)*E3;

  E3 = 0.0;
  for(i=0;i<sz;i++){
    for(j=0;j<sz;j++){
      E3 += q[i]*q[j];
    }
  }
  E3 *= (4.0*const2/etha3);
  energy -= E3; 


  return energy;

}




double Vdw_LJ(VECTOR* r,                                               /* Inputs */
              VECTOR* g,
              VECTOR* m,
              VECTOR* f,
              MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
              int sz,double* epsilon, double* sigma,
              int nexcl, int* excl1, int* excl2, double* scale,
              MATRIX3x3* box,int rec_deg,int pbc_deg,
              double etha,int is_cutoff, double R_on, double R_off,
              int& time, vector< vector<triple> >& images, vector<triple>& central_translation,
              double* dr2, double dT2, int& is_update     /* Parameters */
             ){

//  cout<<"In MB-LJ potential!!!\n";

  int count = 0; // for triples and central translations
  int excl = 0;  // for exclusions
  int i,j;
  double SW,sig,eps,en;
  VECTOR dSW,fmod;
  int na,nb,nc;
  int xshift,yshift,zshift;
  double r2,energy;
  VECTOR dij,rij,gij,mij,f1,f2,f12;
  VECTOR tv1,tv2,tv3,g1,g2,g3,tv;
  MATRIX3x3 tp;
  MATRIX3x3 I;  I.identity();
  double fscl = 1.0;
  double tscale = 1.0;

  // Reciprocal vectors
  box->get_vectors(tv1,tv2,tv3);
  box->inverse().T().get_vectors(g1,g2,g3);

  //------------------ Initialize forces and stress -----------------
  energy = 0.0;
  for(i=0;i<sz;i++){ f[i] = 0.0; }
  at_stress = 0.0;
  fr_stress = 0.0;
  ml_stress = 0.0;

  double R_skin = 2.0;
  double R_skin2 = 4.0;
  Cell cl(tv1,tv2,tv3,R_off+R_skin);
  triple central_translation_ij;
  vector<triple> images_ij;
  int update_displ2 = 0;

  int is_new = 0;
  if(images.size()==0){ is_new = 1; }
  double scl_const = 1.0; // use 1 if full i,j loops are used: (i in [1,sz]) and (j in [1,sz])
                          // use 2 if half i,j loops are used: (i in [1,sz]) and (j in [i,sz])
  double scl;
  double scl1;

//  cout<<"In MB-VDW\n";
//  cout<<"sz = "<<sz<<endl;

  for(i=0;i<sz;i++){
    for(j=i;j<sz;j++){

//     cout<<"i = "<<i<<" j = "<<j;

     // Here we use the ordering of the excl1 and excl2 arrays which corresponds to 
     // given way of double loop organization
     if(((i==excl1[excl])&&(j==excl2[excl]))||
        ((j==excl1[excl])&&(i==excl2[excl]))  
       ){ scl = scl_const * scale[excl]; excl++; }
     else{ scl = scl_const; }

//    cout<<" excl = "<<excl<<" scl = "<<scl;

// Even more efficient version

     int is_update = 0;
/*
     if(is_new) { is_update = 1; }
     else{
       if(images[count].size()==0){ is_update = 1; }
       else{
         // Calculate new central translation:
//         int nxshift = floor(rij*g1+0.5);
//         int nyshift = floor(rij*g2+0.5);
//         int nzshift = floor(rij*g3+0.5);
//         if((nxshift!=central_translation[count].n1)||
//            (nyshift!=central_translation[count].n2)||
//            (nzshift!=central_translation[count].n3)
//           ){ is_update = 1; }
//         else{
           double displ2  = dr2[i] + dr2[j] + dT2 ;
           if(displ2>R_skin2){ is_update = 1;  }
//         }// else
       }// else
     }// else
*/
 
     is_update = 1; // this is no Verlet list option - for checking

      if(is_update && is_new){
        // Update list
        rij = r[i] - r[j];
        cl.calculate(rij,images_ij,central_translation_ij);
        images.push_back(images_ij);
        central_translation.push_back(central_translation_ij);
        update_displ2 = 1;
      }
      else if(is_update && !is_new){
        rij = r[i] - r[j]; 
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

//      cout<<" n_images = "<<n_images;


      for(int im=0;im<n_images;im++){
//        cout<<" ( im="<<im;
        na = images_ij[im].n1;
        nb = images_ij[im].n2;
        nc = images_ij[im].n3;
//        cout<<" : "<<na<<","<<nb<<","<<nc;

        if((na==xshift) && (nb==yshift) && (nc==zshift) && (i==j)){ /* skip this - atom interacts with itself*/   }
        else{
          if((na==xshift) && (nb==yshift) && (nc==zshift)){ scl1 = scl; }
          else{ scl1 = scl_const; }

//          cout<<", scl1="<<scl1;

          if(scl1>0.0){

          tv = (na*tv1 + nb*tv2 + nc*tv3);

          rij = r[i] - r[j] - tv;
          gij = g[i] - g[j] - tv;
          mij = m[i] - m[j] - tv;
          r2 = rij.length2();

          SW = 1.0; dSW = 0.0;
          VECTOR rj = r[j]+tv;
          if(is_cutoff){ SWITCH(r[i],rj,R_on,R_off,SW,dSW); }
          if(SW>0.0){
            f1 = f2 = 0.0;
            sig = (sigma[i]*sigma[j]);
            eps = (epsilon[i]*epsilon[j]);
            en = Vdw_LJ(r[i],rj,f1,f2,sig,scl1*eps);
            energy += SW*en;
            f12 = (SW*f1 - en*dSW);
            f[i] += f12;
            f[j] -= f12; 

             tp.tensor_product(rij , f12);   at_stress += tscale*tp;
             tp.tensor_product(gij , f12);   fr_stress += tscale*tp;
           }
           }// scl1>0.0
         }//else  - not gamma point
//        cout<<")";
      }// for im
//      cout<<endl;
    }// for j
  }// for i

//  time_ij++;
//  time = time_ij;

  return energy;
}


double Vdw_LJ1(VECTOR* r,                                               /* Inputs */
               VECTOR* g,
               VECTOR* m,
               VECTOR* f,
               MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
               int sz,double* epsilon, double* sigma,
               int nexcl, int* excl1, int* excl2, double* scale,
               MATRIX3x3* box,int rec_deg,int pbc_deg,
               double etha,int is_cutoff, double R_on, double R_off,
               int& time, vector< vector<quartet> >& at_neib, vector<triple>& central_translation,
               double* dr2, double dT2, int& is_update     /* Parameters */
              ){

//  This version is based on lists for each atom

  int count = 0; // for triples and central translations
  int excl = 0;  // for exclusions
  int i,j,k;
  double SW,sig,eps,en;
  VECTOR dSW,fmod;
  int na,nb,nc;
  int xshift,yshift,zshift;
  double r2,energy;
  VECTOR dij,rij,gij,mij,f1,f2,f12;
  VECTOR tv1,tv2,tv3,g1,g2,g3,tv;
  MATRIX3x3 tp;
  MATRIX3x3 I;  I.identity();
  double fscl = 1.0;
  double tscale = 1.0;

  // Reciprocal vectors
  box->get_vectors(tv1,tv2,tv3);
  box->inverse().T().get_vectors(g1,g2,g3);

  //------------------ Initialize forces and stress -----------------
  energy = 0.0;
  for(i=0;i<sz;i++){ f[i] = 0.0; }
  at_stress = 0.0;
  fr_stress = 0.0;
  ml_stress = 0.0;

  double R_skin = 2.0;
  double R_skin2 = 4.0;
// Cell version
//  Cell cl(tv1,tv2,tv3,R_off+R_skin);
//  cl.update_vlist(sz,r,at_neib,central_translation);

// Nlist version

  double cellx,celly,cellz;
  cellx = 1.0*(R_off+R_skin);
  celly = 1.0*(R_off+R_skin);
  cellz = 1.0*(R_off+R_skin);
//  vector< vector<quartet> > nlist;

//  cout<<"time = "<<time<<endl;
//  cout<<"time%4 ="<<time%4<<endl;
//  if(time==0||time==20){
//    cout<<"update nlist\n";
    if(at_neib.size()>0) { at_neib.clear(); }
    make_nlist_auto(sz,r,*box,cellx,celly,cellz,R_off+R_skin,at_neib);
//    time=0;
//  }
//  time++;


//  cout<<"update_vlist is done\n";
//  cout<<"at_neib.size() = "<<at_neib.size()<<endl;
//  cout<<"central_translation.size() = "<<central_translation.size()<<endl;
//  cout<<"sz = "<<sz<<endl;

//  exit(0);

  int update_displ2 = 0;

  double scl_const = 1.0; // use 1 if full i,j loops are used: (i in [1,sz]) and (j in [1,sz])
                          // use 2 if half i,j loops are used: (i in [1,sz]) and (j in [i,sz])
  double scl;
  double scl1;
  int is_next_j;
  int excl_ij;

//  j = -1;
  for(i=0;i<sz;i++){
    int sz1 = at_neib[i].size();
//    cout<<" i = "<<i<<" sz1 = "<<sz1<<endl;
    j = -1;
    int cnt = 0;
    for(k=0;k<sz1;k++){

//       xshift = central_translation[count].n1;
//       yshift = central_translation[count].n2;
//       zshift = central_translation[count].n3;
       // Here the ordering is important 
//       if(at_neib[i][k].j!=j){
//         count++; // this variable numerates the pairs of atom to consider
//       }
       is_next_j = (at_neib[i][k].j!=j);
       j = at_neib[i][k].j;
       cnt++;

//     cout<<"i = "<<i<<" k = "<<k<<" j = "<<j<<endl;      
       if(is_next_j){
         int is_scaled = 0;
         for(excl=0;excl<nexcl;excl++){ 
           if( ((i==excl1[excl]) && (j==excl2[excl]))||
               ((j==excl1[excl]) && (i==excl2[excl]))
             ){ is_scaled = 1; excl_ij = excl; }
         }// for excl

         if(is_scaled){ scl = scl_const*scale[excl_ij]; }
         else{ scl = scl_const; }

//       cout<<"i = "<<i<<" k = "<<k<<" j = "<<j<<" scl = "<<scl<<" excl_ij = "<<excl_ij<<" cnt="<<cnt<<endl;
        
       cnt  = 0;
       }// is_next_j
       else{  // The j is the same - only different translation => keep scaling the same
       }

/*
     // Here we use the ordering of the excl1 and excl2 arrays which corresponds to
     // given way of double loop organization
     if(((i==excl1[excl])&&(j==excl2[excl]))
//        ((j==excl1[excl])&&(i==excl2[excl]))
       ){ 
       scl = scl_const * scale[excl];
       if(is_next_j){  excl++;}
     }
     else{ scl = scl_const; }

//     cout<<"i = "<<i<<" k = "<<k<<" j = "<<j<<" scl = "<<scl<<" excl = "<<excl;

*/
// Even more efficient version
        int is_update = 1;

        na = at_neib[i][k].n1;
        nb = at_neib[i][k].n2;
        nc = at_neib[i][k].n3;

//        if((na==xshift) && (nb==yshift) && (nc==zshift) && (i==j)){ /* skip this - atom interacts with itself*/   }
        if((at_neib[i][k].is_central==1) && (i==j)){ /* skip this - singular case*/ }
        else{
//          if((na==xshift) && (nb==yshift) && (nc==zshift)){ scl1 = scl; }
          if(at_neib[i][k].is_central){ scl1 = scl; }
          else{ scl1 = scl_const; }

//          cout<<" scl1 = "<<scl1<<" na = "<<na<<" nb = "<<nb<<" nc = "<<nc;

          if(scl1>0.0){

          tv = (na*tv1 + nb*tv2 + nc*tv3);

          rij = r[i] - r[j] - tv;
          gij = g[i] - g[j] - tv;
          mij = m[i] - m[j] - tv;
          r2 = rij.length2();

          SW = 1.0; dSW = 0.0;
          VECTOR rj = r[j]+tv;
          if(is_cutoff){ SWITCH(r[i],rj,R_on,R_off,SW,dSW); }
          if(SW>0.0){
            f1 = f2 = 0.0;
            sig = (sigma[i]*sigma[j]);
            eps = (epsilon[i]*epsilon[j]);
            en = Vdw_LJ(r[i],rj,f1,f2,sig,scl1*eps);
            energy += SW*en;
            f12 = (SW*f1 - en*dSW);
            f[i] += f12;
            f[j] -= f12;

             tp.tensor_product(rij , f12);   at_stress += tscale*tp;
             tp.tensor_product(gij , f12);   fr_stress += tscale*tp;
           }
           }// scl1>0.0
         }//else  - not self-interaction
//         cout<<endl;
      }// for k ~ (j,n1,n2,n3)
  }// for i

  return energy;
}


double Vdw_LJ2_no_excl(VECTOR* r,                                               /* Inputs */
                       VECTOR* g,
                       VECTOR* m,
                       VECTOR* f,
                       MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
                       int sz,double* epsilon, double* sigma,
                       int nexcl, int* excl1, int* excl2, double* scale,
                       MATRIX3x3* box,int rec_deg,int pbc_deg,
                       double etha,int is_cutoff, double R_on, double R_off,        
                       int& tim, vector< vector<excl_scale> >& excl_scales
                      ){

//  This version is based on lists for each atom
// Moreover we do not store the lists of each atom - do not story any
// big lists - to get speed up
// In other words we open the make_nlist_auto function here - inside 
// the execution loops
//
// excl_scales - the scaling constants for exclusion interactions for all atoms
// excl_scales[0] - contains exclusion scales for pairs 0-1_0, 0-2_0, ... 0-n_0, where n_0 - is a number of exclusions with first atom at_indx1
// excl_scales[1] - ...
// ...                  
// excl_scales[N_atoms]

/* Debug with this:
  try{   }
  catch(char *e){ printf("Exception Caught: %s\n",e); exit(0);   }
*/

//  int count = 0; // for triples and central translations
  int excl = 0;  // for exclusions
  int i,j,a,b,c,k;
  double SW,sig,eps,en;
  VECTOR dSW,fmod;
  double r2,energy;
  VECTOR dij,rij,gij,mij,f1,f2,f12;
  VECTOR t1,t2,t3,g1,g2,g3,tv;
  MATRIX3x3 tp;
  MATRIX3x3 I;  I.identity();
  double fscl = 1.0;
  double tscale = 1.0;

  // Reciprocal vectors
  box->get_vectors(t1,t2,t3);
  box->inverse().T().get_vectors(g1,g2,g3);

  //>..................... From make_nlist_auto part 1 ...................<
  // Calculate constants
  double Roff2 = R_off * R_off;

  double cellx,celly,cellz;
  cellx = 1.1*(R_off);
  celly = 1.1*(R_off);
  cellz = 1.1*(R_off);

  int Nat = sz;
  vector<VECTOR> R;
  vector<quartet> initT; // initial translations
  VECTOR* tmp; tmp = new VECTOR[Nat];


  apply_pbc(*box,Nat,r,tmp,initT); // now all tmp in simulation cell, no outside
                                   // r = tmp + (*box)*initT, that is initT are the integer vectors (in units of box)
                                   // to restore original r from tmp


  triple minb,maxb;

  find_min_shell(t1,t2,t3,g1,g2,g3,R_off,minb,maxb);




  // Bounding box is based on maximal atomic displacements
  int indx = 0;
  vector<quartet> globqt;
  VECTOR minr,maxr; minr = maxr = r[0];
  for(a=minb.n1;a<=maxb.n1;a++){
    for(b=minb.n2;b<=maxb.n2;b++){
      for(c=minb.n3;c<=maxb.n3;c++){
        for(i=0;i<Nat;i++){

          VECTOR t = tmp[i] + a*t1 + b*t2 + c*t3;  R.push_back(t);
          if(t.x<=minr.x){ minr.x = t.x; } if(t.x>=maxr.x){ maxr.x = t.x; }
          if(t.y<=minr.y){ minr.y = t.y; } if(t.y>=maxr.y){ maxr.y = t.y; }
          if(t.z<=minr.z){ minr.z = t.z; } if(t.z>=maxr.z){ maxr.z = t.z; }

          quartet locqt;
          locqt.j = i; locqt.n1 = a; locqt.n2 = b; locqt.n3 = c;
          globqt.push_back(locqt);  // R[indx] = tmp[i] + (*box)*(globqt[indx]), that is globqt[indx] is an integer vector
                                    // in units of box to produce images R[indx] from the folded coordinates
          indx++;
        }// for i
      }// for c
    }// for b
  }// for a

  // Number of partitions in corresponding direction
  VECTOR maxdr; maxdr = maxr - minr;
  int Nx = (floor(maxdr.x/cellx)+1);if(Nx<1){ Nx = 1; cout<<"Error: Nx<1\n"; exit(0);}
  int Ny = (floor(maxdr.y/celly)+1);if(Ny<1){ Ny = 1; cout<<"Error: Ny<1\n"; exit(0);}
  int Nz = (floor(maxdr.z/cellz)+1);if(Nz<1){ Nz = 1; cout<<"Error: Nz<1\n"; exit(0);}
  int Ncells = Nx*Ny*Nz;

  vector<int> dummy;
  vector<quartet> qdummy;
  vector<int> at2cell_indx(indx,-1);            // is a complex index of the given atom
  vector< vector<int> > cell2at(Ncells,dummy);// cell2at[i] - contains complex indexes of all atoms in cell i

  // Calculate neighbors of each cell (sub-cell)
  // we use serial indexes of both central cell and its neighbors
  vector< vector<int> > neibc(Ncells,dummy);     // indexes of neighboring cells for given cell index

  for(c=0;c<Ncells;c++){
    form_neibc(c,neibc[c],Nx,Ny,Nz,cellx,celly,cellz,R_off);
  }

  // Calculate mappings between atom indexes and cell (sub-cell) indexes
  for(i=0;i<indx;i++){
    VECTOR diff = R[i] - minr;  // position of
    triple trp;
    trp.n1 = floor(diff.x/cellx); if(trp.n1<0.0){ trp.n1 = 0; } if(trp.n1>Nx){ trp.n1 = Nx;}
    trp.n2 = floor(diff.y/celly); if(trp.n2<0.0){ trp.n2 = 0; } if(trp.n2>Ny){ trp.n2 = Ny;}
    trp.n3 = floor(diff.z/cellz); if(trp.n3<0.0){ trp.n3 = 0; } if(trp.n3>Nz){ trp.n3 = Nz;}
    c = Nz*Ny*trp.n1 + Nz*trp.n2 + trp.n3;
    if(c>=Ncells){ cout<<"This can not be!!! Dynamical error: c>=Ncells - array not allocated.\n"; exit(0); }

    at2cell_indx[i] = c;
    cell2at[c].push_back(i);

  }

  int cc = (maxb.n3-minb.n3+1)*(maxb.n2-minb.n2+1)*(-minb.n1) + (maxb.n3-minb.n3+1)*(-minb.n2) + (-minb.n3); // index of central image/cell (not sub-cell)

  //>..................... From make_nlist_auto part 1 ...................<

  //------------------ Initialize forces and stress -----------------
  energy = 0.0;
  for(i=0;i<sz;i++){ f[i] = 0.0; }
  at_stress = 0.0;
  fr_stress = 0.0;
  ml_stress = 0.0;

  int update_displ2 = 0;
  double scl_const = 1.0; // use 1 if full i,j loops are used: (i in [1,sz]) and (j in [1,sz])
                          // use 2 if half i,j loops are used: (i in [1,sz]) and (j in [i,sz])
  double scl;
  double scl1;
  int is_next_j;
  int excl_ij;


  //>..................... From make_nlist_auto part 2 ...................<

  int shift = Nat*cc;

  for (int at_indx1=0;at_indx1<Nat;at_indx1++){

    int i1  = shift+at_indx1; //complex index of real atom at_indx1
    int ci1 = at2cell_indx[i1]; // complex index of cell to which atom i belongs
    int sz1 = neibc[ci1].size();// number of neighbor cells of the cell with index ci1

    int excl_indx = -1;
    // find the index of the array of exclusions involving atom at_indx1 as a first atom
    for(int excl=0;excl<excl_scales.size();excl++){
      if(excl_scales[excl][0].at_indx1==at_indx1){ excl_indx = excl; break; }
    }



    for(c=0;c<sz1;c++){
      int ci2 = neibc[ci1][c]; // one of the neighboring cells of the cell l
      int sz2 = cell2at[ci2].size();// number of atoms in the cell with index ci2

      for(int a=0;a<sz2;a++){ // iterations over atoms in cell ci2
        int i2 = cell2at[ci2][a]; // complex index of atom a of the cell ci2
        int at_indx2 = i2 % Nat;  // real index of atom, corresponding to the atom with complex index i2

        if(at_indx2>=at_indx1){
        VECTOR dR = R[i1] - R[i2];

        double modR = dR.x*dR.x;
        if(modR<=Roff2){
          modR += dR.y*dR.y;
          if(modR<=Roff2){
            modR += dR.z*dR.z;
            if(modR<=Roff2){

              int n1,n2,n3,is_central; 
// This is correct - no worries!!!
              // globqt[i2] - initT[at_indx2] = is the integer vector (in box units) to produce given image R[i2]
              // from the original coordinate of the atom i :  r[i]
              // similarly we have for globqt[i1] - initT[at_indx2], however because we currently considering atom i1 as a 
              // center of our coordinate system we have globqt[i1] = (0,0,0)
              // The difference between two translation vectors is given below and has a meaning of the relative displacement
              // of atom r[at_indx2] with respect to r[at_indx1] to give the same vector as between R[i2] and R[i1] found in close contact
              n1 = globqt[i2].n1 + (initT[at_indx1].n1 - initT[at_indx2].n1);
              n2 = globqt[i2].n2 + (initT[at_indx1].n2 - initT[at_indx2].n2);
              n3 = globqt[i2].n3 + (initT[at_indx1].n3 - initT[at_indx2].n3);
// And this is too!!!
              is_central = 0;
              if(n1==0 && n2==0 && n3==0) { is_central = 1; }

              //============ Calculate scaling ========================  
              int is_scaled = 0; 
              scl = scl_const;
              if(excl_indx>-1){ 
                for(excl=0;excl<excl_scales[excl_indx].size();excl++){
                  if(excl_scales[excl_indx][excl].at_indx2==at_indx2){
                    is_scaled = 1; scl = scl_const * excl_scales[excl_indx][excl].scale; break;
                  }
                }
              }


              //============= Calculation part =========================
              if((is_central==1) && (at_indx1==at_indx2)){ /* skip this - singular case*/ }
              else{
              // go here in cases: a) is_central == 0  (both at_indx1==at_indx2 or at_indx1!=at_indx2)
              //                   b) is_central == 1  and at_indx1!=at_indx2

//                if(is_central==0){ // Comment this line to get full treatment (include central exclusions)
                // The above line (with corresponding closing one) should remain commented to include major contribution
                // non-excluded pair vdw interactions within central cell, exclusions are simply scaled.

                // The reason why the above line can be uncommented is for the RESPA - separating time scales
                // so slowly-varying non-central part (this whole function) can be done not so often, while
                // the central part (see below function) needs to be updated more often

                if(is_central){ scl1 = scl; }
                else{ scl1 = scl_const; }

                if(scl1*scl1>0.0){

                  tv = (n1*t1 + n2*t2 + n3*t3);
                  gij = g[at_indx1] - g[at_indx2] - tv;

                  SW = 1.0; dSW = 0.0;
                  VECTOR rj = r[at_indx2]+tv;
                  if(is_cutoff){ SWITCH(r[at_indx1],rj,R_on,R_off,SW,dSW); }
                  if(SW>0.0){
                    f1 = f2 = 0.0;
                    sig = (sigma[at_indx1]*sigma[at_indx2]);
                    eps = (epsilon[at_indx1]*epsilon[at_indx2]);
                    en = Vdw_LJ(r[at_indx1],rj,f1,f2,sig,scl1*eps);
                    energy += SW*en;
                    //cout<<" at_indx1 = "<<at_indx1<<" at_indx2 = "<<at_indx2<<" sig = "<<sig<<" eps = "<<eps<<" en = "<<en<<" SW = "<<SW<<endl;
                    f12 = (SW*f1 - en*dSW);
                    f[at_indx1] += f12;
                    f[at_indx2] -= f12;

                    tp.tensor_product(rij , f12);   at_stress += tscale*tp;
                    tp.tensor_product(gij , f12);   fr_stress += tscale*tp;
                  }


                }// scl1>0.0
//                }// is_central==0

              }//else  - not self-interaction
         
            }//zik
          }//yik
        }//xik
        }// at_indx2>=at_indx1
      }// for a

    }// for c
  }// for at_indx1

  delete [] tmp;


  return energy;
}

double Vdw_LJ2_excl(VECTOR* r,                                               /* Inputs */
                    VECTOR* g,
                    VECTOR* m,
                    VECTOR* f,
                    MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
                    int sz,double* epsilon, double* sigma,
                    int nexcl, int* excl1, int* excl2, double* scale,
                    MATRIX3x3* box,int rec_deg,int pbc_deg,
                    double etha,int is_cutoff, double R_on, double R_off,
                   int& tim, vector< vector<excl_scale> >& excl_scales
                   ){
/**********************************************************************************
 This function calculates only exclusions interacting via vdw LJ potential with
 appropriately scaled parameters

**********************************************************************************/

  double SW,sig,eps,en;
  VECTOR dSW,fmod;
  double r2,energy;
  VECTOR dij,rij,gij,mij,f1,f2,f12;
  MATRIX3x3 tp;
  double tscale = 1.0;


  //>..................... From make_nlist_auto part 1 ...................<
  //------------------ Initialize forces and stress -----------------
  energy = 0.0;
  for(int i=0;i<sz;i++){ f[i] = 0.0; }
  at_stress = 0.0;
  fr_stress = 0.0;
  ml_stress = 0.0;



  for(int excli=0;excli<excl_scales.size();excli++){
     for(int exclij=0;exclij<excl_scales[excli].size();exclij++){

      int at_indx1 = excl_scales[excli][exclij].at_indx1;
      int at_indx2 = excl_scales[excli][exclij].at_indx2;
      double scl1 = excl_scales[excli][exclij].scale;

      //============= Calculation part =========================
      if(at_indx1==at_indx2){} // skip this - singular case
      else{
        if(scl1>0.0){          
          VECTOR gij = g[at_indx1] - g[at_indx2];
          SW = 1.0; dSW = 0.0;
          if(is_cutoff){ SWITCH(r[at_indx1],r[at_indx2],R_on,R_off,SW,dSW); }
          if(SW>0.0){
            f1 = f2 = 0.0;
            sig = (sigma[at_indx1]*sigma[at_indx2]);
            eps = (epsilon[at_indx1]*epsilon[at_indx2]);
            en = Vdw_LJ(r[at_indx1],r[at_indx2],f1,f2,sig,scl1*eps);
            energy += SW*en;
            f12 = (SW*f1 - en*dSW);
            f[at_indx1] += f12;
            f[at_indx2] -= f12;

            tp.tensor_product(rij , f12);   at_stress += tscale*tp;
            tp.tensor_product(gij , f12);   fr_stress += tscale*tp;
          }
        }// scl1>0.0
      }//else  - not self-interaction
    }// for exclij
  }// for excli



  return energy;

}


double LJ_Coulomb(VECTOR* r, VECTOR* g, VECTOR* m, VECTOR* f,
                  MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress,
                  int sz,double* epsilon, double* sigma,double* q,int is_cutoff, double R_on, double R_off,
                  int nexcl, int* excl1, int* excl2, double* scale)
{
/** 
  This function computes all electrostatic and vdW interactions in a collection of particles

  \param[in] r The pointer to the array containing the positions of all particles (atoms)
  \param[in] g The pointer to the array containing the center of mass positions of all particles (groups)
  \param[in] m The pointer to the array containing the center of mass positions of all particles (molecules)
  \param[out] f The pointer to the array containing the forces acting on all atoms
  \param[out] at_stress Atomic stress tensor
  \param[out] fr_stress Fragmental stress tensor
  \param[out] ml_stress Molecular stress tensor
  \param[in] sz The number of particles in the system (size of the arrays r and f)
  \param[in] epsilon The array of the vdW epsilon parameters (stength of interaction) for all atoms
  \param[in] sigma The array of the vdW sigma parameters (atomic size) for all atoms
  \param[in] q The array of the atomic charges
  \param[in] is_cutoff The flag telling whether the cutoff for vdW is utilized. Note - no cutoffs to electrostatic potentials 
  will be applied
  \param[in] R_on The cutoff start distance (switching function is 1)
  \param[in] R_off The cutoff stop distance (switching function is 0)
  \param[in] nexcl The number of exclusions
  \param[in] excl1 excl1[i] is the index of the first atom excuded in the pair i
  \param[in] excl2 excl2[i] is the index of the second atom excuded in the pair i
  \param[in] scale scale[i] is the scaling constant for the exclusion i

  --- So far this function only works for cluster, no PBC ---

  --- Atomic units of length, energy, parameters.etc. are assumed -----

*/

  int excl = 0;  // for exclusions
  int i,j;
  double SW,sig,eps,en;
  VECTOR dSW,fmod;
  int na,nb,nc;
  int xshift,yshift,zshift;
  double r2,energy;
  VECTOR dij,rij,gij,mij,f1,f2,f12;
  VECTOR tv1,tv2,tv3,g1,g2,g3,tv;
  MATRIX3x3 tp;
  MATRIX3x3 I;  I.identity();
  double fscl = 1.0;
  double tscale = 1.0;


  //------------------ Initialize forces and stress -----------------
  energy = 0.0;
  for(i=0;i<sz;i++){ f[i] = 0.0; }
  at_stress = 0.0;
  fr_stress = 0.0;
  ml_stress = 0.0;

  double scl_const = 1.0; // use 1 if full i,j loops are used: (i in [1,sz]) and (j in [1,sz])
                          // use 2 if half i,j loops are used: (i in [1,sz]) and (j in [i,sz])
  double scl;

  for(i=0;i<sz;i++){
    for(j=i;j<sz;j++){
     // Get scaling constant
     // Here we use the ordering of the excl1 and excl2 arrays which corresponds to 
     // given way of double loop organization
     if(((i==excl1[excl])&&(j==excl2[excl]))||
        ((j==excl1[excl])&&(i==excl2[excl]))  
       ){ scl = scl_const * scale[excl]; excl++; }
     else{ scl = scl_const; }

 
     if(i==j){ /* skip this - atom interacts with itself*/   }
     else{
       if(scl>0.0){
         
         rij = r[i] - r[j];
         gij = g[i] - g[j];
         mij = m[i] - m[j];
         r2 = rij.length2();

         // Vdw part
         SW = 1.0; dSW = 0.0;
         VECTOR rj = r[j];
         if(is_cutoff){ SWITCH(r[i],rj,R_on,R_off,SW,dSW); }
         if(SW>0.0){
           f1 = f2 = 0.0;
           sig = (sigma[i]*sigma[j]);
           eps = (epsilon[i]*epsilon[j]);
           en = Vdw_LJ(r[i],rj,f1,f2,sig,scl*eps);
           energy += SW*en;
           f12 = (SW*f1 - en*dSW);
           f[i] += f12;
           f[j] -= f12; 

           tp.tensor_product(rij , f12);   at_stress += tscale*tp;
           tp.tensor_product(gij , f12);   fr_stress += tscale*tp;
         }

         // Electrostatic part
         f1 = f2 = 0.0;
         en = Elec_Coulomb(r[i],r[j],f1,f2,q[i],q[j],1.0, 0.0 );
         energy += scl*en;
         f12 = scl*f1;
         f[i] += f12;
         f[j] -= f12; 


       }// scl>0.0
     }//else  - not self-atom

   }// for j
  }// for i


  return energy;
}


}// namespace libpot
}// liblibra

