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
  \file Grid_functions.cpp
  \brief The file implements some auxiliary functions for grid operations
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <sstream>
#endif 

#include "Grid_functions.h"

/// liblibra namespace
namespace liblibra{


using namespace std;

/// libdyn namespace
namespace libdyn{

/// libwfcgrid namespace
namespace libwfcgrid{


//--------------------- General ---------------------

int compute_imapping(vector<int>& inp, vector<int>& npts){
/**

  By giving the index of the point along each dimension, we want to compute 
  the index of the point in the overall scheme

  The formula is like:

  1D: i = i_0 
  2D: i = i_0*n_1 + i_1
  3D: i = i_0*n_1*n_2 + i_1*n_2 + i_2 

  ...

*/

  if(inp.size()!=npts.size()){
    cout<<"ERROR in compute_imapping: inp.size()= "<<inp.size()<<" but npts.size()= "<<npts.size()<<endl;
    cout<<"Exiting now...\n";
    exit(0);
  }

  int ndof = npts.size();

  int i = inp[0];
  for(int dof=1; dof<ndof; dof++){    
    i *= npts[dof];
    i += inp[dof]; 
  }
  
  return i;
}



vector<int> compute_mapping(int indx, vector<int>& npts){
/** 
  Maps an integer to a vector according to the grid dimensions   
*/

  int i;
  int ndof = npts.size();
  vector<int> sizes;

  int sz = 1;
  for(i=0; i<ndof; i++){
    sizes.push_back(sz);
    sz = sz * npts[ndof-1-i];
  }

  vector<int> res;
  int _indx = indx;
  for(i=0; i<ndof; i++){
    int rem = _indx % sizes[ndof-1-i];
    int ni = int((_indx - rem) / sizes[ndof-1-i]);
    _indx = _indx - ni*sizes[ndof-1-i];

    res.push_back(ni);
  }
 
  return res;

}

vector<int> compute_hyperplane(vector<int>& npts, int idim_const, int ipt_const){
/**
  This function computes the set of grid points (their indices) that
  belong to a given hyperplan of a multi-dimensional grid.

  npts - grid dimensions
  idim_const - index of the constrainted dimension
  ipt_const - index of the point on the constrained dimension
*/

  int i;
  int ndof = npts.size();

  // Compute the total number of point on the grid
  int tot_npts = 1;
  for(i=0; i<ndof; i++){
    tot_npts = tot_npts * npts[i];
  }

  // Now compute the indices of the points belonging to the hyperplane
  vector<int> res;
  vector<int> vec;
  for(i=0; i<tot_npts; i++){
    vec = compute_mapping(i, npts);
    if(vec[idim_const] == ipt_const){ res.push_back(i); }
  }

  return res;
}


vector<vector<int> > compute_mapping(vector<vector<int> >& inp, vector<int>& npts){
/**
     The ordering is like this...

     DOF0      DOF1 ...  Index

                x         0
      x         x         1
                x         2
 
                x         3
      x         x         4
                x         5

                x         6
      x         x         7
                x         8

*/


  int sz,i,ipt;
  sz = inp.size();             // how many vectors are in

  if(sz==0){                   // starting with an empty container
      
    vector<vector<int> > res; 

    for(ipt=0; ipt<npts[0]; ipt++){ // all points along a given dimension

      vector<int> res_i(1, ipt);
      res.push_back(res_i);      
    }
    return compute_mapping(res, npts);
  }
  else{
  
    int lvl = inp[0].size();     // level is the same as the length of each vector

    if (lvl<npts.size()){

      vector<vector<int> > res; 
 
      for(i=0; i<sz; i++){  // take each initial input vectors
   
          vector<int> res_i = inp[i];
          res_i.push_back(0);

          for(ipt=0; ipt<npts[lvl]; ipt++){ // all points along a given dimension
            res_i[lvl] = ipt; 
            res.push_back(res_i);
          }// for ipt
      
      }// for i

      return compute_mapping(res, npts);

    }
    else{  return inp; }
  }

}



int find_grid_size(double xmin,double xmax, double dx){
/**
  \brief Compute the minimal number of points that is a power of 2 and
  encloses the selected interval with the given spacing between the points

  \param[in] xmin The minimal (leftmost) boundary of the grid
  \param[in] xmax The maximal (rightmost) boundary of the grid
  \param[in] dx The spacing between the grid points.

*/
  int sz = int((xmax - xmin)/dx) + 1;

  // Automatically expands the grid to the next closest power of 2
  int n1 = 2; do{ n1 *= 2; }while(sz>n1); sz = n1;

  return sz;
}

CMATRIX init_grid(double xmin,double xmax, double dx){
/**
  \brief Initialize 1D grid

  The number of points that is a power of 2 and encloses the selected interval
  with the given spacing between the points is computed, and then these points 
  will be generated as a Nx x 1 complex-valued matrix (elements are actually real-valued)

  \param[in] xmin The minimal (leftmost) boundary of the grid
  \param[in] xmax The maximal (rightmost) boundary of the grid
  \param[in] dx The spacing between the grid points.

*/

  int sz = find_grid_size(xmin, xmax, dx);
                                     
  CMATRIX x(sz,1);
  for(int i=0;i<sz;i++){ x.M[i] = complex<double>(xmin + i*dx,0.0);   }
  return x;
}

//------------------ 1D specific --------------------------

void init_gauss_1D(vector<CMATRIX>& wfc,CMATRIX& X,double x_,double px_,double dx, int nstates, int occ_state, complex<double> scl){
/**
  \brief Initialize a Gaussian wavepacket on a 1D grid
  \param[out] wfc Is a list of complex matrices (vectors), each containing numerical wavefunction for different electronic state
  \param[in] X is the complex matrix(vector) containing the grid points
  \param[in] x_ Position of the center of the Gaussian wavepacket
  \param[in] px_ Momentum of the Gaussian wavepacket
  \param[in] dx Spread (distribution width) of the spatial component of the Gaussian wavepacket
  \param[in] nstates The number of electronic states for which to initialize the wavefunction
  \param[in] occ_state Index of the occupied electronic state on which the wavepacket is initialized

  G(x) = [ (1/(2.0*pi*dx^2))^(1/4) ] * exp(-((x-x_)/(2*dx))^2 + i*(x-x_)*px_)

  P(x) = |G(x)|^2 = [ (1/(2.0*pi*dx^2)^2) ] * exp( -(x-x_)^2/(2*dx^2) )

  That is according to: https://en.wikipedia.org/wiki/Normal_distribution,  
  dx - corresponds to standard deviation (in the classical distribution)
  dx^2 - variance (in the classical distribution)

  That is, this wavepacket would correspond to the probability density (classical coordinates)
  that are generated like this:

  x_i = x_ * dx * rnd.normal()


  Other connections:     
  2*a = 1/2*dx^2 =>  a = 1/(2*dx)^2 , where a is such that:  alpha/2 = a + i*b, where a and b are defined in:

  (1) Heller, E. J. Guided Gaussian Wave Packets. Acc. Chem. Res. 2006, 39, 127–134.  
  (2) Akimov, A. V.; Prezhdo, O. V. Formulation of Quantized Hamiltonian Dynamics in Terms of Natural Variables. J. Chem. Phys. 2012, 137, 224115.

*/

  // Get the size of the 1D grid
  int Nx = X.n_elts; 

  // Allocate memory, if not yet done
  wfc = vector<CMATRIX>(nstates,CMATRIX(Nx,1));


  // Constants
  const double nrm = pow((1.0/(2.0*M_PI*dx*dx)),0.25);
  const complex<double> one(0.0, 1.0);


  // Copute wfc values at the grid points
  for(int nx=0;nx<Nx;nx++){ 

    for(int st=0;st<nstates;st++){
      if(st==occ_state){

        double deltx = 0.5*(X.M[nx].real() - x_)/dx;

        double c1 = -deltx*deltx;
        double c2 = px_*(X.M[nx].real() - x_);

        wfc[st].M[nx*1+0] = scl * nrm * exp(c1) * (cos(c2)+one*sin(c2));
      }
      else{ wfc[st].M[nx*1+0] = 0.0; }

    }// for st
  }// for nx

}// init_gauss_1D


void init_gauss_1D(vector<CMATRIX>& wfc,CMATRIX& X,double x_,double px_,double dx, int nstates, int occ_state){

  complex<double> one(1.0, 0.0);
  init_gauss_1D(wfc, X, x_, px_, dx, nstates, occ_state, one); 

}


void add_gauss_1D(vector<CMATRIX>& wfc,CMATRIX& X,double x_,double px_,double dx, int nstates, int occ_state, complex<double> weight){
/**
  \brief Adds a Gaussian wavepacket with a given weight to a 1D grid

  \param[out] wfc Is a list of complex matrices (vectors), each containing numerical wavefunction for different electronic state
  \param[in] X is the complex matrix(vector) containing the grid points
  \param[in] x_ Position of the center of the Gaussian wavepacket
  \param[in] px_ Momentum of the Gaussian wavepacket
  \param[in] dx Spread (distribution width) of the spatial component of the Gaussian wavepacket
  \param[in] nstates The number of electronic states for which to initialize the wavefunction
  \param[in] occ_state Index of the occupied electronic state on which the wavepacket is initialized

  G(x) = [ (1/(2.0*pi*dx^2))^(1/4) ] * exp(-((x-x_)/(2*dx))^2 + i*(x-x_)*px_)

  P(x) = |G(x)|^2 = [ (1/(2.0*pi*dx^2)^2) ] * exp( -(x-x_)^2/(2*dx^2) )

  That is according to: https://en.wikipedia.org/wiki/Normal_distribution,  
  dx - corresponds to standard deviation (in the classical distribution)
  dx^2 - variance (in the classical distribution)

  That is, this wavepacket would correspond to the probability density (classical coordinates)
  that are generated like this:

  x_i = x_ * dx * rnd.normal()


  Other connections:     
  2*a = 1/2*dx^2 =>  a = 1/(2*dx)^2 , where a is such that:  alpha/2 = a + i*b, where a and b are defined in:

  (1) Heller, E. J. Guided Gaussian Wave Packets. Acc. Chem. Res. 2006, 39, 127–134.  
  (2) Akimov, A. V.; Prezhdo, O. V. Formulation of Quantized Hamiltonian Dynamics in Terms of Natural Variables. J. Chem. Phys. 2012, 137, 224115.

*/

  // Get the size of the 1D grid
  int Nx = X.n_elts; 

  // The memory should be allocated already
  // wfc = vector<CMATRIX>(nstates,CMATRIX(Nx,1));

  // Constants
  const double nrm = pow((1.0/(2.0*M_PI*dx*dx)),0.25);
  const complex<double> one(0.0, 1.0);


  // Copute wfc values at the grid points
  for(int nx=0;nx<Nx;nx++){ 

    for(int st=0;st<nstates;st++){
      if(st==occ_state){

        double deltx = 0.5*(X.M[nx].real() - x_)/dx;

        double c1 = -deltx*deltx;
        double c2 = px_*(X.M[nx].real() - x_);

        wfc[st].M[nx*1+0] += weight * nrm * exp(c1) * (cos(c2)+one*sin(c2));
      }
      else{ wfc[st].M[nx*1+0] += 0.0; }

    }// for st
  }// for nx

}// init_gauss_1D




void add_ho_1D(vector<CMATRIX>& wfc, CMATRIX& X, int nu, double x_, double px_, complex<double> weight, int occ_state, int alpha){
/**
  \brief Adds a moving! Harmonic Oscillator (HO) with a given weight to a 1D grid wavefunction

  \param[out] wfc Is a list of complex matrices (vectors), each containing numerical wavefunction for different electronic state
  \param[in] X is the complex matrix(vector) containing the grid points
  \param[in] nu Quantum number of the HO basis function
  \param[in] x_ Position of the center of the HO basis function
  \param[in] px_ Momentum of the HO basis wavepacket 
  \param[in] weight The weight with which the basis function is added to the grid
  \param[in] occ_state Index of the electronic state to which the wavepacket is added
  \param[in] alpha The parameter related to the reference HO Hamiltonian:  alpha = sqrt(k*mu/hbar^2)
             Where:  H = -hbar^2 /(2*mu) d^2/dx^2  + 1/2 * k * (x-x_)^2 


  This function adds:

  weight * HO_nu(x-x_) * exp(i*px_*(x-x_))

  HO_nu(x) = N_nu * H_nu(sqrt(alpha) * (x)) * exp(-1/2 * alpha * x^2 ) 

  N_nu = 1/sqrt(2^nu * nu!)  * (alpha/pi)^(1/4) - normalization factor

  H_nu(ksi):  H_{nu+1}(ksi) - 2 ksi*H_{nu}(ksi) + 2*nu *H_{nu-1}(ksi) = 0  - Hermite polynomial

*/

  // Get the size of the 1D grid
  int Nx = X.n_elts; 

  // Allocate memory, if not yet done
  // wfc = vector<CMATRIX>(nstates,CMATRIX(Nx,1));

  // Constants
  const double nrm = (1.0/sqrt(pow(2.0, nu) * FACTORIAL(nu)) ) * pow((alpha/M_PI),0.25);
  const complex<double> one(0.0, 1.0);

  double H, dH;

  // Copute wfc values at the grid points
  for(int nx=0;nx<Nx;nx++){ 
    
    double ksi = sqrt(alpha) * (X.M[nx].real() - x_);
    HERMITE(nu, ksi, H, dH);

    double c2 = px_*(X.M[nx].real() - x_);

    wfc[occ_state].M[nx] += weight * nrm *  H * exp(-0.5*ksi*ksi) * (cos(c2)+one*sin(c2));

  }// for nx

}// init_ho_1D






void print_1D(CMATRIX& X,vector<CMATRIX>& PSI,string prefix, int frame){
/**
  \brief Printing populations of all states in different files (one state per file)
  \param[in] X The grid points (Nx x 1 matrix)
  \param[in] PSI The wavefunctions for all states to be printed out. This is a vector of Nx x 1 matrices
  \param[in] prefix The common part of the names of the files to which the wfc will be printed out
  \param[in] frame The integer index to be added into the file name - for instance, when printing snapshots of dynamical simulations
  
  All files have a common prefix, which also includes indexing via frame parameter, but all 
  files have different suffixes - depending on the electronic state
*/

  std::string sframe,filename;
  stringstream ss(stringstream::in | stringstream::out);
  ss << frame;  ss >> sframe;  
  filename = prefix + ".frame."+ sframe;

  print_1D(X,PSI,filename);

}

void print_1D(CMATRIX& X,vector<CMATRIX>& PSI,string filename){
/**
  \brief Printing populations of all states in different files (one state per file)
  \param[in] X The grid points (Nx x 1 matrix)
  \param[in] PSI The wavefunctions for all states to be printed out. This is a vector of Nx x 1 matrices
  \param[in] filename The common part of the names of the files to which the wfc will be printed out
  
  All files have a common prefix, which also includes indexing via frame parameter, but all 
  files have different suffixes - depending on the electronic state
*/

  
  int nstates = PSI.size();
  int Nx = X.n_elts;

  for(int st=0;st<nstates;st++){

    std::string sst, filename1; 
    stringstream ss(stringstream::in | stringstream::out);
    ss << st;  ss >> sst;
    filename1 = filename + ".st."+sst;
    

    std::ofstream out(filename1.c_str(),ios::out);

    for(int nx=0;nx<Nx;nx++){
      out<<real(X.M[nx])<<"  "<<real(std::conj(PSI[st].M[nx])*PSI[st].M[nx])<<endl;
    }// for nx

    out.close();

  }// for st
  
}

void ft_1D(vector<CMATRIX>& psi,vector<CMATRIX>& reci_psi,int opt,
           double xmin,double kxmin,double dx){
/**
  \brief Do 1D cfft (complex fast Fourier transform) for each state independently
  \param[in,out] psi  Reals space wavefunctions: nstates x Nx x 1 - that is nstates complex-valued matrices Nx x 1 each
  \param[in,out] reci_psi  Reciprocal space wavefunctions: nstates x Nx x 1 - that is nstates complex-valued matrices Nx x 1 each
  \param[in] opt option to control the direction of transformation:  opt = 1:  psi ---> reci_psi;  opt = 2:  psi <--- reci_psi
  \param[in] xmin The lower boundary of the real-space grid
  \param[in] kxmin The lower boundary of the reciprocal-space grid
  \param[in] dx The real-space grid point spacing 

  Note: the actual meaning of psi and reci_psi is defined by opt parameter
*/

  int nstates = psi.size();
  int Nx = psi[0].n_rows;

  CMATRIX in(Nx,1);
  CMATRIX out(Nx,1);

  for(int i=0;i<nstates;i++){ 

    in = psi[i];
    if(opt==1){  cfft1(in,out,xmin,kxmin,dx); }
    else if(opt==2){ inv_cfft1(in,out,xmin,kxmin,dx); }

    reci_psi[i] = out;

  }// for i
}


//========================= 2D specific ================================

void init_gauss_2D(vector<CMATRIX>& wfc,
                   CMATRIX& X,double x_,double px_,double dx,
                   CMATRIX& Y,double y_,double py_,double dy,
                   int nstates, int occ_state){
/**
  \brief Initialize a Gaussian wavepacket on a 2D grid
  \param[out] wfc Is a list of complex matrices (vectors), each containing numerical wavefunction for different electronic state
  \param[in] X is the complex matrix(vector) containing the grid points along X direction
  \param[in] x_ Position of the center of the Gaussian wavepacket in X direction
  \param[in] px_ Momentum of the Gaussian wavepacket in X direction
  \param[in] dx Spread (distribution width) of the spatial component of the Gaussian wavepacket along the axis X
  \param[in] Y is the complex matrix(vector) containing the grid points along Y direction
  \param[in] y_ Position of the center of the Gaussian wavepacket in Y direction
  \param[in] py_ Momentum of the Gaussian wavepacket in Y direction
  \param[in] dy Spread (distribution width) of the spatial component of the Gaussian wavepacket along the axis Y
  \param[in] nstates The number of electronic states for which to initialize the wavefunction
  \param[in] occ_state Index of the occupied electronic state on which the wavepacket is initialized

  G(x,y) = G(x)*G(y)
  G(x) = [ (1/(2.0*pi*dx^2))^(1/4) ] * exp(-((x-x_)/(2*dx))^2 + i*(x-x_)*px_)
  G(y) = [ (1/(2.0*pi*dy^2))^(1/4) ] * exp(-((y-y_)/(2*dy))^2 + i*(y-y_)*py_)

*/


  int Nx = X.n_elts; 
  int Ny = Y.n_elts;

  CMATRIX tmp(Nx,Ny);
  vector<CMATRIX> wfcx(nstates,CMATRIX(Nx,1));
  vector<CMATRIX> wfcy(nstates,CMATRIX(Ny,1));

  init_gauss_1D(wfcx,X,x_,px_,dx,nstates,occ_state);
  init_gauss_1D(wfcy,Y,y_,py_,dy,nstates,occ_state);


  for(int nx=0;nx<Nx;nx++){
    for(int ny=0;ny<Ny;ny++){

      for(int st=0;st<nstates;st++){

          wfc[st].M[nx*Ny+ny] = wfcx[st].M[nx] * wfcy[st].M[ny];

      }// for st

    }// for ny
  }// for nx

}


void print_2D(CMATRIX& X,CMATRIX& Y,vector<CMATRIX>& PSI,string prefix, int frame){
/**
  \brief Printing populations of all states in different files (one state per file)
  \param[in] X The grid points (Nx x 1 matrix) along X direction
  \param[in] Y The grid points (Ny x 1 matrix) along Y direction
  \param[in] PSI The wavefunctions for all states to be printed out. This is a vector of Nx x Ny matrices
  \param[in] prefix The common part of the names of the files to which the wfc will be printed out
  \param[in] frame The integer index to be added into the file name - for instance, when printing snapshots of dynamical simulations
  
  All files have a common prefix, which also includes indexing via frame parameter, but all 
  files have different suffixes - depending on the electronic state
*/


  std::string sframe,filename;
  stringstream ss(stringstream::in | stringstream::out);
  ss << frame;  ss >> sframe;  
  filename = prefix + ".frame."+ sframe;

  print_2D(X,Y,PSI,filename);

}


void print_2D(CMATRIX& X,CMATRIX& Y,vector<CMATRIX>& PSI,string filename){
/**
  \brief Printing populations of all states in different files (one state per file)
  \param[in] X The grid points (Nx x 1 matrix) along X direction
  \param[in] Y The grid points (Ny x 1 matrix) along Y direction
  \param[in] PSI The wavefunctions for all states to be printed out. This is a vector of Nx x Ny matrices
  \param[in] filename The common part of the names of the files to which the wfc will be printed out
  
  All files have a common prefix, which also includes indexing via frame parameter, but all 
  files have different suffixes - depending on the electronic state
*/


  int nstates = PSI.size();
  int Nx = X.n_elts;
  int Ny = Y.n_elts;

  for(int st=0;st<nstates;st++){

    std::string sst, filename1; 
    stringstream ss(stringstream::in | stringstream::out);
    ss << st;  ss >> sst;
    filename1 = filename + ".st."+sst;
    

    std::ofstream out(filename1.c_str(),ios::out);


    for(int nx=0;nx<Nx;nx++){
      for(int ny=0;ny<Ny;ny++){
        out<<real(X.M[nx])<<"  "<<real(Y.M[ny])<<"  "<<real(std::conj(PSI[st].M[nx*Ny+ny])*PSI[st].M[nx*Ny+ny])<<endl;
      }// for ny
      out<<"\n";
    }// for nx

    out.close();

  }// for st
  
}


void ft_2D(vector<CMATRIX>& psi,vector<CMATRIX>& reci_psi,int opt,
           double xmin,double ymin,double kxmin,double kymin,double dx,double dy){
/**
  \brief Do 2D cfft (complex fast Fourier transform) for each state independently
  \param[in,out] psi  Reals space wavefunctions: nstates x Nx x Ny - that is nstates complex-valued matrices Nx x Ny each
  \param[in,out] reci_psi  Reciprocal space wavefunctions: nstates x Nx x Ny - that is nstates complex-valued matrices Nx x Ny each
  \param[in] opt option to control the direction of transformation:  opt = 1:  psi ---> reci_psi;  opt = 2:  psi <--- reci_psi
  \param[in] xmin The lower boundary of the real-space grid
  \param[in] kxmin The lower boundary of the reciprocal-space grid
  \param[in] dx The real-space grid point spacing 

  Note: the actual meaning of psi and reci_psi is defined by opt parameter
*/

  int nstates = psi.size();
  int Nx = psi[0].n_rows;
  int Ny = psi[0].n_cols;

  CMATRIX in(Nx,Ny);
  CMATRIX out(Nx,Ny);

  for(int i=0;i<nstates;i++){ 

    in = psi[i];
    if(opt==1){  cfft1_2D(in,out,xmin,ymin,kxmin,kymin,dx,dy); }
    else if(opt==2){ inv_cfft1_2D(in,out,xmin,ymin,kxmin,kymin,dx,dy); }

    reci_psi[i] = out;

  }// for i
}


} // namespace libwfcgrid
} // namespace libdyn
}// liblibra
