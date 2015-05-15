#include "Wfcgrid.h"

namespace libdyn{
namespace libwfcgrid{


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
void Wfcgrid::allocate(){

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
  expH = vector<CMATRIX>(nstates,psi); 
  expK = vector<CMATRIX>(nstates,psi);

  
  X = new CMATRIX(1,Nx);        // grid of r-points
  Y = new CMATRIX(1,Ny);        // grid of r-points
  Kx = new CMATRIX(1,Nx);
  Ky = new CMATRIX(1,Ny);



}// allocate


void Wfcgrid::init_grid(){

  // all numbers are assumed to be defined by this time

  // Initialize grids and wavefunctions
  *X = libdyn::libwfcgrid::init_grid(xmin,xmax,dx);                                        // real space
  for(int nx=0;nx<Nx;nx++){ Kx->M[nx] = kxmin + nx/((double)Nx*dx);}   // reciprocal space
  *Y = libdyn::libwfcgrid::init_grid(ymin,ymax,dy);                                        // real space
  for(int ny=0;ny<Ny;ny++){ Ky->M[ny] = kymin + ny/((double)Ny*dy);}   // reciprocal space            
  cout<<"Grids are initialized\n";

}// init_grid



// Constructor
Wfcgrid::Wfcgrid(double minx_, double maxx_, double dx_, double miny_, double maxy_, double dy_, int nstates_){

  init_numbers(minx_, maxx_, dx_, miny_, maxy_, dy_, nstates_);
  allocate();
  init_grid();

}


void Wfcgrid::init_wfc(double x0, double y0, double px0, double py0, double dx0, double dy0, int init_state){

  init_gauss_2D(PSI, *X, x0, px0, dx0,    *Y, y0, py0, dy0,  nstates, init_state);

  cout<<"Wavefunction is initialized\n";


}// init_wavefunction



void Wfcgrid::print_map(std::string prefix, int snap, int state){

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


}


void Wfcgrid::update_potential(Hamiltonian_Model& ham){
// This function recomputes the Hamiltonian for all points
// working in atomic units: hbar = 1

  vector<double> q(2,0.0);

  // Precompute H, d_ij, ... along grid
  for(int nx=0;nx<Nx;nx++){
    for(int ny=0;ny<Ny;ny++){

      // For all r points of the grid compute local(potential energy) part
      q[0] = real(X->M[nx]);
      q[1] = real(Y->M[ny]);

      ham.set_q(q);
      ham.compute();

      for(int nst=0;nst<nstates;nst++){        
        for(int nst1=0;nst1<nstates;nst1++){               

          H[nst][nst1].M[nx*Ny+ny] = ham.Hvib(nst, nst1);  
                                                        
          Dx[nst][nst1].M[nx*Ny+ny] = ham.D(nst, nst1, 0);
          Dy[nst][nst1].M[nx*Ny+ny] = ham.D(nst, nst1, 1); 

        }// for nst1
      }// for nst

    }// for ny
  }// for nx

}


void Wfcgrid::update_propagator_K(double dt,double m0){
// This function updates real- and recoprocal-space propagators (diagonal matrices)
// working in atomic units: hbar = 1

  // Precompute H, d_ij, ... along grid
  for(int nx=0;nx<Nx;nx++){
    for(int ny=0;ny<Ny;ny++){

      double kx_ = Kx->M[nx].real();
      double ky_ = Ky->M[ny].real();
      double argg = -(2.0*M_PI*M_PI/m0)*(kx_*kx_ + ky_*ky_)*dt;

      complex<double> scl(std::cos(argg),std::sin(argg));

      for(int nst=0;nst<nstates;nst++){    expK[nst].M[nx*Ny+ny] = scl;    }//for nst

    }// for ny
  }// for nx


}// update_propagator


void Wfcgrid::update_propagator(double dt,double m0){
// This function updates real- and recoprocal-space propagators (diagonal matrices)
// working in atomic units: hbar = 1

  // Precompute H, d_ij, ... along grid
  for(int nx=0;nx<Nx;nx++){
    for(int ny=0;ny<Ny;ny++){

      for(int nst=0;nst<nstates;nst++){  

        double argg = -dt*H[nst][nst].M[nx*Ny+ny].real();
        expH[nst].M[nx*Ny+ny] = complex<double>(std::cos(argg), std::sin(argg));  // exp(-i*H*dt)

      }//for nst

    }// for ny
  }// for nx


}// update_propagator




}// namespace libwfcgrid
}// namespace libdyn
