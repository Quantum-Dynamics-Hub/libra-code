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

#include "Wfcgrid.h"

namespace libdyn{
namespace libwfcgrid{


void Wfcgrid::init_numbers(double minx_, double maxx_, double dx_, int nstates_){

  nstates = nstates_;

  xmin = minx_;
  xmax = maxx_;
  dx = dx_;

  kxmin = -0.5/dx;

  // Expand the grid to be power of 2
  Nx = find_grid_size(xmin,xmax,dx);
  Ny = 1;

  cout<<"Grid size is calculated\n";
  cout<<"Nx = "<<Nx<<endl;
  cout<<"Ny = "<<Ny<<endl;


}// init_numbers


void Wfcgrid::init_numbers(double minx_, double maxx_, double dx_, double miny_, double maxy_, double dy_, int nstates_){

  nstates = nstates_;

  xmin = minx_;
  xmax = maxx_;
  dx = dx_;

  ymin = miny_;
  ymax = maxy_;
  dy = dy_;

  kxmin = -0.5/dx;
  kymin = -0.5/dy;

  // Expand the grid to be power of 2
  Nx = find_grid_size(xmin,xmax,dx);
  Ny = find_grid_size(ymin,ymax,dy);

  cout<<"Grid size is calculated\n";
  cout<<"Nx = "<<Nx<<endl;
  cout<<"Ny = "<<Ny<<endl;


}// init_numbers



// Memory allocator
void Wfcgrid::allocate_1D(){

  // all numbers are assumed to be defined by this time
  cout<<"Nx = "<<Nx<<endl;
  cout<<"Ny = "<<Ny<<endl;
  cout<<"nstates = "<<nstates<<endl;


  // Allocate arrays
  CMATRIX psi(Nx,Ny);  psi = 0.0; // is a matrix placeholder

  PSI = vector<CMATRIX>(nstates,psi); // wavefunction  nstates x Nx x Ny
  reciPSI = vector<CMATRIX>(nstates,psi);   // same as PSI but in Fourier (reciprocal) space with 1.0*kmin
  DtreciPSI = vector<CMATRIX>(nstates,psi); // d/dt * reciPSI - time derivative of PSI in reciprocal space
  DxPSI = vector<CMATRIX>(nstates,psi);     // d/dx * PSI in real space (up to factor)
  KxreciPSI = vector<CMATRIX>(nstates,psi); // Kx * reciPSI  -  d/dx * PSI in k space (up to factor)

  // Potential-related variables
  H  = vector< vector<CMATRIX> >(nstates, vector<CMATRIX>(nstates, psi) ); // Hamiltonian
  Dx = vector< vector<CMATRIX> >(nstates, vector<CMATRIX>(nstates, psi) ); // for each i,j pair - a separate x projection for all r-points of the grid
  expH = vector< vector<CMATRIX> >(nstates,vector<CMATRIX>(nstates, psi) ); 
  expK = vector<CMATRIX>(nstates,psi);

  
  X = new CMATRIX(1,Nx);        // grid of r-points
  Kx = new CMATRIX(1,Nx);


}// allocate_1D


void Wfcgrid::allocate_2D(){

  // all numbers are assumed to be defined by this time
  cout<<"Nx = "<<Nx<<endl;
  cout<<"Ny = "<<Ny<<endl;
  cout<<"nstates = "<<nstates<<endl;


  // Allocate arrays
  CMATRIX psi(Nx,Ny);  psi = 0.0; // is a matrix placeholder

  PSI = vector<CMATRIX>(nstates,psi); // wavefunction  nstates x Nx x Ny
  reciPSI = vector<CMATRIX>(nstates,psi);   // same as PSI but in Fourier (reciprocal) space with 1.0*kmin
  DtreciPSI = vector<CMATRIX>(nstates,psi); // d/dt * reciPSI - time derivative of PSI in reciprocal space
  DxPSI = vector<CMATRIX>(nstates,psi);     // d/dx * PSI in real space (up to factor)
  DyPSI = vector<CMATRIX>(nstates,psi);     // d/dx * PSI in real space (up to factor)  
  KxreciPSI = vector<CMATRIX>(nstates,psi); // Kx * reciPSI  -  d/dx * PSI in k space (up to factor)
  KyreciPSI = vector<CMATRIX>(nstates,psi); // Ky * reciPSI  -  d/dy * PSI in k space (up to factor)

  // Potential-related variables
  H  = vector< vector<CMATRIX> >(nstates, vector<CMATRIX>(nstates, psi) ); // Hamiltonian
  Dx = vector< vector<CMATRIX> >(nstates, vector<CMATRIX>(nstates, psi) ); // for each i,j pair - a separate x projection for all r-points of the grid
  Dy = vector< vector<CMATRIX> >(nstates, vector<CMATRIX>(nstates, psi) ); // for each i,j pair - a separate y projection for all r-points of the grid
  expH = vector< vector<CMATRIX> >(nstates,vector<CMATRIX>(nstates, psi) ); 
  expK = vector<CMATRIX>(nstates,psi);

  
  X = new CMATRIX(1,Nx);        // grid of r-points
  Y = new CMATRIX(1,Ny);        // grid of r-points
  Kx = new CMATRIX(1,Nx);
  Ky = new CMATRIX(1,Ny);



}// allocate_2D


void Wfcgrid::init_grid_1D(){

  // all numbers are assumed to be defined by this time

  // Initialize grids and wavefunctions
  *X = libdyn::libwfcgrid::init_grid(xmin,xmax,dx);                                        // real space
  for(int nx=0;nx<Nx;nx++){ Kx->M[nx] = kxmin + nx/((double)Nx*dx);}   // reciprocal space
  cout<<"Grids are initialized\n";

}// init_grid_1D


void Wfcgrid::init_grid_2D(){

  // all numbers are assumed to be defined by this time

  // Initialize grids and wavefunctions
  *X = libdyn::libwfcgrid::init_grid(xmin,xmax,dx);                                        // real space
  for(int nx=0;nx<Nx;nx++){ Kx->M[nx] = kxmin + nx/((double)Nx*dx);}   // reciprocal space
  *Y = libdyn::libwfcgrid::init_grid(ymin,ymax,dy);                                        // real space
  for(int ny=0;ny<Ny;ny++){ Ky->M[ny] = kymin + ny/((double)Ny*dy);}   // reciprocal space            
  cout<<"Grids are initialized\n";

}// init_grid_2D



// 1D Constructor

Wfcgrid::Wfcgrid(double minx_, double maxx_, double dx_, int nstates_){

  init_numbers(minx_, maxx_, dx_, nstates_);
  allocate_1D();
  init_grid_1D();

}
// 2D Constructor
Wfcgrid::Wfcgrid(double minx_, double maxx_, double dx_, double miny_, double maxy_, double dy_, int nstates_){

  init_numbers(minx_, maxx_, dx_, miny_, maxy_, dy_, nstates_);
  allocate_2D();
  init_grid_2D();

}


void Wfcgrid::init_wfc_1D(double x0, double px0, double dx0, int init_state){

  init_gauss_1D(PSI, *X, x0, px0, dx0, nstates, init_state);

  cout<<"Wavefunction is initialized\n";


}// init_wavefunction


void Wfcgrid::init_wfc_2D(double x0, double y0, double px0, double py0, double dx0, double dy0, int init_state){

  init_gauss_2D(PSI, *X, x0, px0, dx0,    *Y, y0, py0, dy0,  nstates, init_state);

  cout<<"Wavefunction is initialized\n";


}// init_wavefunction


void Wfcgrid::print_wfc_1D(std::string prefix, int snap, int state){

// for 1D profile on XY plane

  std::string filename, snaps, states;
  stringstream ss(stringstream::in | stringstream::out);
  stringstream ss1(stringstream::in | stringstream::out);

  ss  << snap;   ss  >> snaps;
  ss1 << state;  ss1 >> states;


  filename = prefix+".state"+states+".frame"+snaps;
  ofstream out(filename.c_str(),ios::out);


  for(int nx=0;nx<Nx;nx++){

    out<<real(X->M[nx])<<"  "<<real(std::conj(PSI[state].M[nx])*PSI[state].M[nx])<<endl;

  }// for nx

  out.close();


}// print_wfc_1D


void Wfcgrid::print_wfc_2D(std::string prefix, int snap, int state){

// for 2D projections on XY plane

  std::string filename, snaps, states;
  stringstream ss(stringstream::in | stringstream::out);
  stringstream ss1(stringstream::in | stringstream::out);

  ss  << snap;   ss  >> snaps;
  ss1 << state;  ss1 >> states;


  filename = prefix+".state"+states+".frame"+snaps;
  ofstream out(filename.c_str(),ios::out);


  for(int nx=0;nx<Nx;nx++){
    for(int ny=0;ny<Ny;ny++){

      out<<real(X->M[nx])<<"  "<<real(Y->M[ny])<<"  "<<real(std::conj(PSI[state].M[nx*Ny+ny])*PSI[state].M[nx*Ny+ny])<<endl;

    }// for ny
    out<<"\n";
  }// for nx

  out.close();


}// print_wfc_2D


void Wfcgrid::print_reci_wfc_1D(std::string prefix, int snap, int state){

// for 1D profile on XY plane

  std::string filename, snaps, states;
  stringstream ss(stringstream::in | stringstream::out);
  stringstream ss1(stringstream::in | stringstream::out);

  ss  << snap;   ss  >> snaps;
  ss1 << state;  ss1 >> states;


  filename = prefix+".state"+states+".frame"+snaps;
  ofstream out(filename.c_str(),ios::out);


  for(int nx=0;nx<Nx;nx++){

    out<<real(Kx->M[nx])<<"  "<<real(std::conj(reciPSI[state].M[nx])*reciPSI[state].M[nx])<<endl;

  }// for nx

  out.close();

}// print_wfc_1D

void Wfcgrid::print_complex_matrix_1D(CMATRIX& CM, std::string filename){

  ofstream out(filename.c_str(),ios::out);

  for(int nx=0;nx<Nx;nx++){
    out<<real(X->M[nx])<<"  "<<real(Kx->M[nx])<<"  "<<CM.M[nx].real()<<"   "<<CM.M[nx].imag()<<endl;
  }// for nx

  out.close();

}

void Wfcgrid::print_ham_1D(std::string prefix, int i, int j){
  print_complex_matrix_1D(H[i][j],prefix);
}

void Wfcgrid::print_expH_1D(std::string prefix, int i, int j){
  print_complex_matrix_1D(expH[i][j],prefix);
}

void Wfcgrid::print_expK_1D(std::string prefix, int i){
  print_complex_matrix_1D(expK[i],prefix);
}




double Wfcgrid::print_populations_1D(string filename,int snap){

  vector<double> Pop(nstates,0.0);
  double Pop_tot = 0.0;
  double Pop_tot_active = 0.0;

  for(int nx=0;nx<Nx;nx++){
    for(int nst=0;nst<nstates;nst++){

      double pii = real(std::conj(PSI[nst].M[nx])*PSI[nst].M[nx]);

      Pop[nst] += dx*pii;  // population on state nst
      Pop_tot += dx*pii;   // total population

    }// for nst
  }// for nx

  ofstream out(filename.c_str(),ios::app);

  out<<"nsnap= "<<snap;
  // Probabilities of the wavepackets that still remain in active calculation area
  for(int nst=0;nst<nstates;nst++){
    out<<"  P("<<nst<<")= "<<setprecision(5)<<Pop[nst];
  }
  out<<" P_total(active)= "<<Pop_tot<<" | \n";
  Pop_tot_active = Pop_tot;


/*
  // Now the reflection and transition probabilities on all states
  double sum = 0.0;
  for(nst=0;nst<nstates;nst++){
    out<<"  P_re("<<nst<<")= "<<setprecision(5)<<Pops[nst][0];
    out<<"  P_tr("<<nst<<")= "<<setprecision(5)<<Pops[nst][1];
    sum += Pops[nst][0] + Pops[nst][1];
  }
  out<<" P_total(absorbed)= "<<sum<<" | ";
  out<<" P_total(all)= "<<Pop_tot + sum<<endl;
*/
  
  out.close();

  return Pop_tot_active;
}


void Wfcgrid::flux_1D(double xf,vector<double>& res, double m0){
// xf - the point at which flux is computed

  if(res.size()<nstates){ res = vector<double>(nstates,0.0);  }

  // index of the grid point at which the counting plane is installed  
  int i = (xf - real(X->M[0]))/dx;


  // Compute wfc derivatives:
  // Form Kx * reciPSI and Ky * reciPSI
  for(int nst=0;nst<nstates;nst++){
    KxreciPSI[nst] = reciPSI[nst];

    for(int kx=0;kx<Nx;kx++){  
      KxreciPSI[nst].M[kx] *= Kx->M[kx];
    }// kx
  }// for nst

  // Kxrec(k) -> DxPSI(r)
  ft_1D(KxreciPSI,DxPSI,2,xmin,kxmin,dx);


  // Im(i*z) = Im(i(re+i*im)) = re = Re(z)
  // z - z* = (re + i*im) - (re - i*im) = 2*i*im
  // Im(z - z*) = 2 * im
  for(int nst=0;nst<nstates;nst++){  
/*
    // Original expression
    res[nst] = 2.0*M_PI*(0.5*hbar/m0)*imag(std::conj(PSI[nst].M[i])*one*DxPSI[nst].M[i] - 
                                           std::conj(one*DxPSI[nst].M[i])*PSI[nst].M[i]
                                          );
    // Simplified 1
    res[nst] = 4.0*M_PI*(0.5*hbar/m0)*imag(std::conj(PSI[nst].M[i])*one*DxPSI[nst].M[i] );
*/
    // Simplified 2
    res[nst] = (2.0*M_PI/m0) * real(std::conj(PSI[nst].M[i])*DxPSI[nst].M[i] );

  }

}// flux_1D




}// namespace libwfcgrid
}// namespace libdyn
