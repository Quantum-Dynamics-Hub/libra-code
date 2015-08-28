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

#include <sstream>
#include "Grid_functions.h"

using namespace std;


namespace libdyn{
namespace libwfcgrid{


//--------------------- General ---------------------

int find_grid_size(double xmin,double xmax, double dx){
  int sz = int((xmax - xmin)/dx) + 1;

  // Automatically expands the grid to the next closest power of 2
  int n1 = 2; do{ n1 *= 2; }while(sz>n1); sz = n1;

  return sz;
}

CMATRIX init_grid(double xmin,double xmax, double dx){
  int sz = find_grid_size(xmin, xmax, dx);
                                     
  CMATRIX x(sz,1);
  for(int i=0;i<sz;i++){ x.M[i] = complex<double>(xmin + i*dx,0.0);   }
  return x;
}

//------------------ 1D specific --------------------------
//void init_gauss_1D(vector<CMATRIX>& wfc, CMATRIX& X, double x_, double p_, double dx, int nstates, int occ_state);
//void print_1D(CMATRIX& X,vector<CMATRIX>& OUT,string filename);
//void ft_1D(vector<CMATRIX>& psi,vector<CMATRIX>& reci_psi,int opt,double xmin,double kxmin,double dx);


void init_gauss_1D(vector<CMATRIX>& wfc,CMATRIX& X,double x_,double px_,double dx, int nstates, int occ_state){
// Gaussian wavepacket on 1D grid
// G(x) = [ (1/(2.0*pi*dx^2))^(1/4) ] * exp(-((x-X)/(2*dx))^2 + i*x*px)

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
        double c2 = 2.0* px_* deltx;

        wfc[st].M[nx*1+0] = nrm * exp(c1) * (cos(c2)+one*sin(c2));
      }
      else{ wfc[st].M[nx*1+0] = 0.0; }

    }// for st
  }// for nx

}// init_gauss_1D


void print_1D(CMATRIX& X,vector<CMATRIX>& PSI,string prefix, int frame){

  std::string sframe,filename;
  stringstream ss(stringstream::in | stringstream::out);
  ss << frame;  ss >> sframe;  
  filename = prefix + ".frame."+ sframe;

  print_1D(X,PSI,filename);

}

void print_1D(CMATRIX& X,vector<CMATRIX>& PSI,string filename){
  
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
// Do cfft for each state independently
// psi = nstates x Nx x 1 - that is nstates matrices Nx x Ny each
// opt = 1:  psi ---> reci_psi
// opt = 2:  psi <--- reci_psi
// Note: the actual meaning of psi and reci_psi is defined by opt parameter

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

  std::string sframe,filename;
  stringstream ss(stringstream::in | stringstream::out);
  ss << frame;  ss >> sframe;  
  filename = prefix + ".frame."+ sframe;

  print_2D(X,Y,PSI,filename);

}


void print_2D(CMATRIX& X,CMATRIX& Y,vector<CMATRIX>& PSI,string filename){

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
// Do cfft for each state independently
// psi = nstates x Nx x Ny - that is nstates matrices Nx x Ny each
// opt = 1:  psi ---> reci_psi
// opt = 2:  psi <--- reci_psi
// Note: the actual meaning of psi and reci_psi is defined by opt parameter

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

