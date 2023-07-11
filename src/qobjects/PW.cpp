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
  \file PW.cpp
  \brief The file implements basic operations on/with plane-wave objects    
*/

#include "PW.h"
#include "../util/libutil.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;
using namespace libutil;


/// libqobjects namespace
namespace libqobjects{


//------------------ Members of the PW class -------------------------
PW PW::operator-(){
/**
  Negation operator: Returns a new wavefunction with all coefficients c_new[i] = -c_old[i]
*/

  PW res(*this); 
  for(int i=0;i<npw;i++){ res.coeff[i] = -coeff[i]; }
  return res;
}

PW PW::operator+(const PW& m){
/**
  Addition operator: Returns a new wavefunction with all coefficients c_new[i] = c_old[i] + m[i],
  but only if the k-points are the same
*/

  if(!(m.kx==this->kx && m.ky==this->ky && m.kz==this->kz)){  
    cout<<"Error in PW::operator+: can add orbitals only if the k-points are the same!\n";
    exit(0);
  }

  PW res(*this);
  for(int i=0;i<npw;i++){ res.coeff[i] = coeff[i] + m.coeff[i]; }
  return res;

}

PW PW::operator-(const PW& m){
/**
  Subtraction operator: Returns a new wavefunction with all coefficients c_new[i] = c_old[i] - m[i],
  but only if the k-points are the same
*/

  if(!(m.kx==this->kx && m.ky==this->ky && m.kz==this->kz)){  
    cout<<"Error in PW::operator-: can subtract orbitals only if the k-points are the same!\n";
    exit(0);
  }

  PW res(*this);
  for(int i=0;i<npw;i++){ res.coeff[i] = coeff[i] - m.coeff[i]; }
  return res;

}

void PW::operator+=(const PW& m){
/**
  Increment operator: Modifies existing coefficient by: c[i] = c[i] + m[i],
  but only if the k-points are the same
*/

  if(!(m.kx==this->kx && m.ky==this->ky && m.kz==this->kz)){  
    cout<<"Error in PW::operator+=: can increment orbitals only if the k-points are the same!\n";
    exit(0);
  }

  for(int i=0;i<npw;i++){ coeff[i] += m.coeff[i]; }

}

void PW::operator-=(const PW& m){
/**
  Decrement operator: Modifies existing coefficient by: c[i] = c[i] - m[i],
  but only if the k-points are the same
*/

  if(!(m.kx==this->kx && m.ky==this->ky && m.kz==this->kz)){  
    cout<<"Error in PW::operator-=: can increment orbitals only if the k-points are the same!\n";
    exit(0);
  }

  for(int i=0;i<npw;i++){ coeff[i] -= m.coeff[i]; }

}


PW PW::operator/(double num){
/**
  Returns a new PW object with the coefficients divided by a real number: c_new[i] = c_old[i] / num
*/

  PW res(*this);
  for(int i=0;i<npw;i++){ res.coeff[i] /= num;  }
  return res;
}

PW PW::operator/(complex<double> num){
/**
  Returns a new PW object with the coefficients divided by a complex number: c_new[i] = c_old[i] / num
*/

  PW res(*this);
  for(int i=0;i<npw;i++){ res.coeff[i] /= num; }
  return res;
}


PW PW::conj(){
/**
  Returns a PW object which is conjugate to original one
*/

  PW res(*this);
  for(int i=0;i<npw;i++){ res.coeff[i] = std::conj(coeff[i]); }
  return res;
}

void PW::normalize(){
/**
  Normalizes the PW object. Function modifies the coefficients of the original object
*/
  int i;
  double norm = 0.0;
  for(i=0;i<npw;i++){ norm += (std::conj(coeff[i]) * coeff[i] ).real();  }
  norm = sqrt(1.0/norm);
  for(i=0;i<npw;i++){ coeff[i] *= norm; }
}

void PW::complete(){
/** 
  Complete the PW object by adding the complex conjugate part
  This is the storage saving trick, due to c(-G) = conj(G), which is usually applied only to
  gamma-point (k_point==0). We also assume that the very first coefficient coeff[0] corresponds
  to G-point = 0,0,0 This is true for QE wavefunctions
  
  The function changes the original object
*/
  int i;
  if(kx==0 && ky==0 && kz==0){

    // Increase storage
    coeff.resize(2*npw-1);

    // Add new coefficients and compute norm of the 'completed' wavefunction
    double norm = (std::conj(coeff[0]) * coeff[0] ).real();
    for(i=1;i<npw;i++){
      coeff[npw+i] = std::conj(coeff[i]); 
      norm += 2.0*(std::conj(coeff[i]) * coeff[i] ).real();
    }

    // Update the number of the planewaves in new (completed) wfc
    npw =  2*npw - 1;

    // Finally, normalize the completed wfc
    norm = sqrt(1.0/norm);
    for(i=0;i<npw;i++){ coeff[i] *= norm; }

  }
  else{
    cout<<"Error in PW::complete() - can't complete non-Gamma-point orbitals\n";
    exit(0);
  }

    
}


template<class T>
PW multiply(T& f, const PW& m1){
  PW res(m1);
  for(int i=0;i<res.npw;i++){ res.coeff[i] *= f; }
  return res;
}

PW operator*(const double& f, const PW& m1){ return multiply<const double>(f,m1); }
PW operator*(const PW& m1, const double& f){ return multiply<const double>(f,m1); }
PW operator*(const float& f, const PW& m1){ return multiply<const float>(f,m1); }
PW operator*(const PW& m1, const float& f){ return multiply<const float>(f,m1); }
PW operator*(const complex<double>& f, const PW& m1){ return multiply<const complex<double> >(f,m1); }
PW operator*(const PW& m1, const complex<double>& f){ return multiply<const complex<double> >(f,m1); }
PW operator*(const complex<float>& f, const PW& m1){ return multiply<const complex<float> >(f,m1); }
PW operator*(const PW& m1, const complex<float>& f){ return multiply<const complex<float> >(f,m1); }


complex<double> operator*(const PW& m1, const PW& m2){
  complex<double> res(0.0,0.0);
  if(m1.npw!=m2.npw){ cout<<"Error: Can not multiply MOs with different basis sizes\n"; exit(0); }
  
  else{   for(int i=0;i<m1.npw;i++){ res += m1.coeff[i] * m2.coeff[i];}  }
  return res;
}





complex<double> I_1D(double kx, double kxp, double gx, double gxp){
/*
    # All arguments are float
    # kx,  gx  - refer to k-point 1 and its corresponding grid points
    # kxp, gxp - refer to k-point 2 and its corresponding grid points
*/

    double zero = 1e-12;
    complex<double> res(0.0, 0.0);
    complex<double> one(0.0, 1.0);

    double delt = kx + gx - kxp - gxp;

    if(fabs(delt) <= zero){   res = complex<double>(1.0, 0.0); }
    else{  
      double argg = -2.0*M_PI*delt;
      res =  -one *  complex<double>( cos(argg) - 1.0 , sin(argg) ) / argg; 
    }

    return res;
}


complex<double> I_3D(VECTOR& k, VECTOR& kp, VECTOR& g, VECTOR& gp){
/*
    # All arguments are VECTOR (float)
    # k,  g  - refer to k-point 1 and its corresponding grid points
    # kp, gp - refer to k-point 2 and its corresponding grid points
*/

    complex<double> res = I_1D(k.x, kp.x, g.x, gp.x) * I_1D(k.y, kp.y, g.y, gp.y) * I_1D(k.z, kp.z, g.z, gp.z);

    return res;

}


CMATRIX pw_overlap(VECTOR& k1, VECTOR& k2, CMATRIX& coeff1, CMATRIX& coeff2, vector<VECTOR>& grid1, vector<VECTOR>& grid2){
/*
    # all k- and g-points are in units of 2*pi/a
    # k1, k2 - are k-point vectors (VECTOR of float) 
    # coeff1 - is a matrix (complex) of coefficeints for all states for given k-point (1), dimensions: npw1 x nbands1
    # coeff2 - is a matrix (complex) of coefficeints for all states for given k-point (2), dimensions: npw2 x nbands2
    # grid1 - a list of vectors for all G-points for given k-point (1): dimension npw1
    # grid2 - a list of vectors for all G-points for given k-point (2): dimension npw2
*/

    int npw1 = coeff1.n_rows;
    int nbands1 = coeff1.n_cols;

    int npw2 = coeff2.n_rows;
    int nbands2 = coeff2.n_cols;


    CMATRIX* S;
    S = new CMATRIX(nbands1, nbands2);  // all orbitals for given pair of k-points (a block of entire matrix)


    // A double sum over the grid points (may be different for the two k-points)
    for(int g1=0; g1<npw1; g1++){

        for(int g2=0; g2<npw2; g2++){

// Original version - neat, but not very efficient
//            complex<double> s = I_3D(k1, k2, grid1[g1], grid2[g2]);

// More efficient version:
            double zero = 1e-12;
            complex<double> one(0.0, 1.0);
            complex<double> s(0.0, 0.0);  // final result

            //======= Work on I_1D (x) ==========
            complex<double> ix(0.0, 0.0);  
            
            double delt = k1.x + grid1[g1].x - k2.x - grid2[g2].x;
            if(fabs(delt) <= zero){   ix = complex<double>(1.0, 0.0); }
            else{  
                double argg = -2.0*M_PI*delt;
                ix =  -one *  complex<double>( cos(argg) - 1.0 , sin(argg) ) / argg; 
            }

            //======= Work on I_1D (y) ==========
            complex<double> iy(0.0, 0.0);  
            
            delt = k1.y + grid1[g1].y - k2.y - grid2[g2].y;
            if(fabs(delt) <= zero){   iy = complex<double>(1.0, 0.0); }
            else{  
                double argg = -2.0*M_PI*delt;
                iy =  -one *  complex<double>( cos(argg) - 1.0 , sin(argg) ) / argg; 
            }

            //======= Work on I_1D (z) ==========
            complex<double> iz(0.0, 0.0);  
            
            delt = k1.z + grid1[g1].z - k2.z - grid2[g2].z;
            if(fabs(delt) <= zero){   iz = complex<double>(1.0, 0.0); }
            else{  
                double argg = -2.0*M_PI*delt;
                iz =  -one *  complex<double>( cos(argg) - 1.0 , sin(argg) ) / argg; 
            }


            s = ix * iy * iz;


            for(int i1=0; i1<nbands1; i1++){
                for(int i2=0; i2<nbands2; i2++){

                    complex<double> tmp = std::conj(coeff1.get(g1,i1)) * s * coeff2.get(g2,i2);

                    S->set(i1,i2, S->get(i1,i2) + tmp);

                }// for i2
            }// for i1
        }// for g2
    }// for g1

    return *S;

}




CMATRIX QE_read_acsii_wfc(std::string filename, int kpt, vector<int>& act_space, int verbose){
/**
   filename - the name of the file from which we'll read this info
   kpt - the index of the k-point that we want to read (starting from 0)
   act_space - is a list of integers that correspond to the the orbitals to be read (starting from 0)
   verbose - the level of verbosity
   
   Returns: 
   wfc - is a npw x nmo  matrix containing certain orbitals in the PW basis   
*/

  // Read all lines
  vector<std::string> A;
  int filesize = read_file(filename,1,A);

  // Find k-point sections
  vector<int> beg,end;
  int nkpts = 0;
  int start = 0;
  int status = 1;
  while(status!=0){
    int b,e;
    status = find_section(A,"<Kpoint.","</Kpoint.",start,filesize,b,e);
    if(status==1){ nkpts++; start = e; beg.push_back(b); end.push_back(e); }
  }
  if(verbose){  cout<<"Number of K-points = "<<nkpts<<endl; }
 
  
  // Now find the positions of the bands for all k-points
  vector< vector<int> > kbeg( nkpts,vector<int>() );
  vector< vector<int> > kend( nkpts,vector<int>() );

  if(kpt<0){ cout<<"Error: The k-point index ( "<<kpt<<" ) can not be less than 0\n"; exit(0); }
  if(kpt>nkpts-1){ cout<<"Error: The k-point index ("<<kpt<<" ) is not allowed. The maximal\
                   value is = "<<nkpts-1<<"\n"; exit(0); }


  start = beg[kpt]; status = 1; int knbnd = 0;

  while(status!=0){
    int b,e;
    status = find_section(A,"<Wfc.","</Wfc.",start,end[kpt],b,e);
    if(status==1) { start = e+1; knbnd++; kbeg[kpt].push_back(b); kend[kpt].push_back(e);}
  }
  if(verbose){ cout<<"Number of bands for k-point "<<kpt<<" is = "<<knbnd<<endl; }


  // Construct MOs:
  int nbands = kbeg[kpt].size();                   // Number of bands(MOs) in each k-point
  int npw =  kend[kpt][0] - kbeg[kpt][0] - 1;      // Number of plane waves in each band

  // Sanity check:
  int sz = act_space.size();
  for(int i=0; i<sz; i++){
    if(act_space[i]<0){ 
      cout<<"Error: Orbital index can not be less than 0, given = "<<act_space[i]<<endl; 
      cout<<"Exiting...\n";
      exit(0);
    }
    if(act_space[i]>nbands-1){ 
      cout<<"Error: Orbital index can not be larger than "<<nbands-1<<", given = "<<act_space[i]<<endl; 
      cout<<"Exiting...\n";
      exit(0);
    }
  }// for i
  nbands = sz; // if all is fine - we'll just read what we need, as dictated by the act_space


  CMATRIX wfc(npw, nbands);
    

  

  // Now finally get the coefficients 
  for(int iband=0;iband<nbands;iband++){
      int band = act_space[iband];

    for(int pw=0;pw<npw;pw++){
      
      int i = kbeg[kpt][band]+1+pw;

      vector<std::string> At;
      split_line(A[i],At,',');
      double re = atof(At[0].c_str());
      double im = atof(At[1].c_str());

      wfc.set(pw, band, re, im);

    }//for npw
  }//for band

  // Free the memory
  A.clear();

  return wfc;
}


MATRIX QE_read_acsii_grid(std::string filename, int verbose){

  // Read all lines
  vector<std::string> A;
  int filesize = read_file(filename,1,A);

  // Find k-point sections
  int beg,end;
  int start = 0;
  int status = 0;
  status = find_section(A,"<grid","</grid>",start,filesize,beg,end);
  int sz = end-beg-1;
  if(verbose){  cout<<"Size of the grid = "<<sz<<endl; }

  MATRIX grid(sz, 3);

  // Read all grid points now
  for(int i=0;i<sz;i++){
    vector<std::string> At;
    split_line(A[i+beg+1],At);

    grid.set(i,0, atof(At[0].c_str()));
    grid.set(i,1, atof(At[1].c_str()));
    grid.set(i,2, atof(At[2].c_str()));

  }// for i

  // Free memory
  A.clear();

  return grid;  
}




vector<CMATRIX> compute_Hprime(CMATRIX& wfc, MATRIX& grid, MATRIX& reci){
/**
   wfc - is a npw x nmo  matrix containing certain orbitals in the PW basis   
   grid - is a npw x 3 matrix containing the grid point coordinates in terms of reciprocal vectors
   reci - is a 3 x 3 matrix of reciprocal vectors - 1-column = bx, 2-nd column = by, 3-rd column = bz

   Returns: 
   hprime - an 3 x nmo x nmo vector of transition dipoles 

*/

  int g_sz = grid.n_rows; 
  int nmo = wfc.n_cols; //used to be sz
  int npw = wfc.n_rows;

  int is_compl = 0;
  if(g_sz!=npw){ 
    if(npw==(2*g_sz-1)){
      cout<<"Warning: Using reconstructed (completed) wavefunction\n";
      is_compl = 1;
    }
    else{
      cout<<"Warning in compute_Hprime: number of plane waves = "
          <<npw<<" is different from the grid size ="<<g_sz<<endl;
    }
    g_sz = min(g_sz,npw);
  }

  // Compute
  complex<double> scl(1.0,0.0); // scaling factor
  complex<double> scl1(0.0,1.0);

  vector<CMATRIX> hprime = vector<CMATRIX>(3, CMATRIX(nmo, nmo));

  complex<double> tmp;
  double gx,gy,gz;  

  
  for(int i=0;i<nmo;i++){
    for(int j=0;j<nmo;j++){

      for(int g=0;g<g_sz;g++){

        tmp = ( std::conj(wfc.get(g,i)) * wfc.get(g, j) );

        gx = grid.get(g,0)*reci.get(0,0) +  grid.get(g,1)*reci.get(0,1) +  grid.get(g,2)*reci.get(0,2);
        gy = grid.get(g,0)*reci.get(1,0) +  grid.get(g,1)*reci.get(1,1) +  grid.get(g,2)*reci.get(1,2);
        gz = grid.get(g,0)*reci.get(2,0) +  grid.get(g,1)*reci.get(2,1) +  grid.get(g,2)*reci.get(2,2);
 
        if(is_compl==0){
          hprime[0].add(i,j, tmp*gx);
          hprime[1].add(i,j, tmp*gy);
          hprime[2].add(i,j, tmp*gz);
        }
        else if(is_compl==1){
          if(g==0){ 
            // This should give zero!
            hprime[0].add(i,j, tmp*gx);
            hprime[1].add(i,j, tmp*gy);
            hprime[2].add(i,j, tmp*gz);
          }
          else{ 
            // Now the Hprime_ matrices are purely imaginary, for the case of gamma-symmetry.
            hprime[0].add(i,j, complex<double>(0.0, 2.0*tmp.real()*gx) );
            hprime[1].add(i,j, complex<double>(0.0, 2.0*tmp.real()*gy) );
            hprime[2].add(i,j, complex<double>(0.0, 2.0*tmp.real()*gz) );
          }
        }// is_compl==1
        
      }// for g

    }// for j
  }// for i

  return hprime;

}


    


}// namespace libqobjects
}// namespace liblibra

