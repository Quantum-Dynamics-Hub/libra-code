#include "Nuclear.h"

/****************************************************************************
  This file contains following functions:

  int Nuclear::add_atom(std::string elt_nam, std::string at_typ, double z, double zeff, double m,VECTOR r)

  VECTOR Nuclear::get_tv(int i,int j, int x_period, int y_period, int z_period,
                         const VECTOR& t1, const VECTOR& t2, const VECTOR& t3)

  void Nuclear::dipole_moment(VECTOR& mu, double& mu_mod)

  double Energy_nucl(Nuclear& mol,vector<int>& fragment)


****************************************************************************/

int Nuclear::add_atom(std::string elt_nam, std::string at_typ, double z, double zeff, double m,VECTOR r){

  // Atomic properties
  elt_name.push_back(elt_nam);
  at_type.push_back(at_typ);
  Z.push_back(z);                      // Atomic number
  Zeff.push_back(zeff);                // Effective charge
  mass.push_back(m*1836.0);            // convert to atomic units (in which electron mass = 1)

  // Dynamic variables
  R.push_back(r);                      // Coordinates
  P.push_back(VECTOR(0.0,0.0,0.0));    // Momentum
  grad.push_back(VECTOR(0.0,0.0,0.0));
  frcs.push_back(VECTOR(0.0,0.0,0.0));
  Mull_charges_gross.push_back(0.0);
  Mull_charges_net.push_back(0.0);

  Nnucl++;

  return Nnucl;
}



VECTOR Nuclear::get_tv(int i,int j, int x_period, int y_period, int z_period,
                       const VECTOR& t1, const VECTOR& t2, const VECTOR& t3){
// Compute the translation vector between AOs i and j that satisfies minimal image convention
// that is R[i] and R[j] + TV are the closest abong other periodic translations

  VECTOR dIdA,dIdB,TV, Rij,Rij0;
  double dist, dist_min;  
  int opt_x,opt_y,opt_z;


  // Here we compute minimal distance convention
  Rij0 = R[i] - R[j]; 
  dist = dist_min = Rij0.length2();
  opt_x = opt_y = opt_z = 0;

  // Consider 7 cases - i know this is not general route, but it is efficient
  // Periodic along x
  if(x_period && !y_period && !z_period){

    for(int nx=-1;nx<=1;nx++){
      Rij = Rij0 - nx*t1;
      dist = Rij.length2(); // this is also not very wise, but it is only for now

      if(dist<dist_min){  dist_min = dist; opt_x = nx; }
    }// for nx

  }// X-periodic


  else if(!x_period && y_period && !z_period){

    for(int ny=-1;ny<=1;ny++){
      Rij = Rij0 - ny*t2;
      dist = Rij.length2(); // this is also not very wise, but it is only for now

      if(dist<dist_min){  dist_min = dist; opt_y = ny; }
    }// for ny

  }// Y-periodic

  else if(!x_period && !y_period && z_period){

    for(int nz=-1;nz<=1;nz++){
      Rij = Rij0 - nz*t3;
      dist = Rij.length2(); // this is also not very wise, but it is only for now

      if(dist<dist_min){  dist_min = dist; opt_z = nz; }
    }// for nz

  }// Z-periodic

  else if(x_period && y_period && !z_period){

    for(int nx=-1;nx<=1;nx++){
      for(int ny=-1;ny<=1;ny++){
        Rij = Rij0 - nx*t1 - ny*t2;
        dist = Rij.length2(); // this is also not very wise, but it is only for now

        if(dist<dist_min){  dist_min = dist; opt_x = nx; opt_y = ny; }
      }// for ny
    }// for nx

  }// XY-periodic

  else if(x_period && !y_period && z_period){

    for(int nx=-1;nx<=1;nx++){
      for(int nz=-1;nz<=1;nz++){
        Rij = Rij0 - nx*t1 - nz*t3;
        dist = Rij.length2(); // this is also not very wise, but it is only for now

        if(dist<dist_min){  dist_min = dist; opt_x = nx; opt_z = nz; }
      }// for nz
    }// for nx

  }// XZ-periodic

  else if(!x_period && y_period && z_period){

    for(int ny=-1;ny<=1;ny++){
      for(int nz=-1;nz<=1;nz++){
        Rij = Rij0 - ny*t2 - nz*t3;
        dist = Rij.length2(); // this is also not very wise, but it is only for now

        if(dist<dist_min){  dist_min = dist; opt_y = ny; opt_z = nz; }
      }// for nz
    }// for ny

  }// YZ-periodic

  else if(x_period && y_period && z_period){

    for(int nx=-1;nx<=1;nx++){
      for(int ny=-1;ny<=1;ny++){
        for(int nz=-1;nz<=1;nz++){
          Rij = Rij0 - nx*t1 - ny*t2 - nz*t3;
          dist = Rij.length2(); // this is also not very wise, but it is only for now

          if(dist<dist_min){  dist_min = dist; opt_x = nx; opt_y = ny; opt_z = nz; }
        }// for nz
      }// for ny
    }// for nx

  }// XYZ-periodic


  TV = opt_x * t1 + opt_y * t2 + opt_z * t3;

  return TV;
}

void Nuclear::dipole_moment(VECTOR& mu, double& mu_mod){

  mu = 0.0; 

  for(int n=0;n<Nnucl;n++){  mu += R[n] * Mull_charges_gross[n];  }// for n

  mu_mod = mu.length();


}

double Energy_nucl(Nuclear& mol,vector<int>& fragment){
/// Compute nuclear energy of a sub-system

  cout<<"Computing nuclear energy:\n";
  

  int I,J;
  int sz = fragment.size();
  double en = 0.0;

  for(int i=0;i<sz;i++){
    I = fragment[i];
    for(int j=i+1;j<sz;j++){
      J = fragment[j];

      en += mol.Zeff[I] * mol.Zeff[J] / (mol.R[I] - mol.R[J]).length();
      
    }// for j
    cout<<"i= "<<i<<" I = "<<I<<" Zeff[I]= "<<mol.Zeff[I]<<"  R[I]= "<<mol.R[I]<<endl;
  }// for i

  return en;

}

