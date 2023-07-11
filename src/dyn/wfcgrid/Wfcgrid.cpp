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
  \file Wfcgrid.cpp
  \brief The file implements some basic methods of the Wfcgrid class: initialization, printing, memory allocation, etc.
    
*/

#include "Wfcgrid.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

/// libwfcgrid namespace
namespace libwfcgrid{


void Wfcgrid::init_numbers(double minx_, double maxx_, double dx_, int nstates_){
/**
  \brief Initialize 1D grid dimensions
  \param[in] minx_ The minimal (leftmost) boundary of the grid
  \param[in] maxx_ The maximal (rightmost) boundary of the grid
  \param[in] dx_ The spacing between the grid points.
  \param[in] nstates_ The number of electronic states to consider

  Computes the number of points in X dimension so that the boundaries are satisfied and the number
  of point is the lowest power of 2 needed to enclose the interval. Also sets other class variables.
*/

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
/**
  \brief Initialize 2D grid dimensions
  \param[in] minx_ The minimal (leftmost) boundary of the grid in X direction
  \param[in] maxx_ The maximal (rightmost) boundary of the grid in X direction
  \param[in] dx_ The spacing between the grid points in X direction
  \param[in] miny_ The minimal (leftmost) boundary of the grid in Y direction
  \param[in] maxy_ The maximal (rightmost) boundary of the grid in Y direction
  \param[in] dy_ The spacing between the grid points in Y direction
  \param[in] nstates_ The number of electronic states to consider

  Computes the number of points in X and Y dimension so that the boundaries are satisfied and the number
  of point is the lowest power of 2 needed to enclose the interval. Also sets other class variables.
*/


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
/**
  \brief Allocates memory for 1D wavefucntion and associted objects

  all numbers are assumed to be defined by this time, so the init_numbers() function must be called first
*/
  
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
/**
  \brief Allocates memory for 2D wavefucntion and associted objects

  all numbers are assumed to be defined by this time, so the init_numbers() function must be called first
*/

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
/**
  \brief Initialize 1D grids (X and Kx points)

  memory is assumed to be allocated by this time. If not, call allocate_1D() function first
*/

  // Initialize grids
  *X = libdyn::libwfcgrid::init_grid(xmin,xmax,dx);                    // real space
  for(int nx=0;nx<Nx;nx++){ Kx->M[nx] = kxmin + nx/((double)Nx*dx);}   // reciprocal space
  cout<<"Grids are initialized\n";

}// init_grid_1D


void Wfcgrid::init_grid_2D(){
/**
  \brief Initialize 2D grids (X, Y and Kx, Ky points)

  memory is assumed to be allocated by this time. If not, call allocate_2D() function first
*/

  // Initialize grids
  *X = libdyn::libwfcgrid::init_grid(xmin,xmax,dx);                    // real space
  for(int nx=0;nx<Nx;nx++){ Kx->M[nx] = kxmin + nx/((double)Nx*dx);}   // reciprocal space
  *Y = libdyn::libwfcgrid::init_grid(ymin,ymax,dy);                    // real space
  for(int ny=0;ny<Ny;ny++){ Ky->M[ny] = kymin + ny/((double)Ny*dy);}   // reciprocal space            
  cout<<"Grids are initialized\n";

}// init_grid_2D



// 1D Constructor
Wfcgrid::Wfcgrid(double minx_, double maxx_, double dx_, int nstates_){
/**
  \brief 1D wavefunction constructors with parameters
  \param[in] minx_ The minimal (leftmost) boundary of the grid in X direction
  \param[in] maxx_ The maximal (rightmost) boundary of the grid in X direction
  \param[in] dx_ The spacing between the grid points in X direction
  \param[in] nstates_ The number of electronic states to consider

  This constructor will: 1) initialize numbers, 2) allocate memory; 3) initialize grids
*/

  init_numbers(minx_, maxx_, dx_, nstates_);
  allocate_1D();
  init_grid_1D();

}

// 2D Constructor
Wfcgrid::Wfcgrid(double minx_, double maxx_, double dx_, double miny_, double maxy_, double dy_, int nstates_){
/**
  \brief 2D wavefunction constructors with parameters
  \param[in] minx_ The minimal (leftmost) boundary of the grid in X direction
  \param[in] maxx_ The maximal (rightmost) boundary of the grid in X direction
  \param[in] dx_ The spacing between the grid points in X direction
  \param[in] miny_ The minimal (leftmost) boundary of the grid in Y direction
  \param[in] maxy_ The maximal (rightmost) boundary of the grid in Y direction
  \param[in] dy_ The spacing between the grid points in Y direction
  \param[in] nstates_ The number of electronic states to consider

  This constructor will: 1) initialize numbers, 2) allocate memory; 3) initialize grids
*/

  init_numbers(minx_, maxx_, dx_, miny_, maxy_, dy_, nstates_);
  allocate_2D();
  init_grid_2D();

}


void Wfcgrid::init_wfc_1D(double x0, double px0, double dx0, int init_state){
/**
  \brief Initialize 1D wavefunction - taken as moving Gaussian wavepacket
  \param[in] x0 Position of the center of the Gaussian wavepacket
  \param[in] px0 Momentum of the Gaussian wavepacket
  \param[in] dx0 Spread (distribution width) of the spatial component of the Gaussian wavepacket
  \param[in] init_state Index of the electronic state on which the wavepacket is initialized

   G(x) = [ (1/(2.0*pi*dx0^2))^(1/4) ] * exp(-((x-x0)/(2*dx0))^2 + i*(x-x0)*px0)

*/

  init_gauss_1D(PSI, *X, x0, px0, dx0, nstates, init_state);

  // PSI(r)->PSI(k)=reciPSI
  ft_1D(PSI,reciPSI,1,xmin,kxmin,dx);


  cout<<"Wavefunction is initialized\n";


}// init_wavefunction


void Wfcgrid::init_wfc_1D(vector<double>& x0, vector<double>& px0, vector<double>& dx0, vector<int>& init_state, vector<complex<double> >& weights){
/**
  \brief Initialize 1D wavefunction - taken as moving Gaussian wavepacket
  \param[in] x0 Position of the center of the Gaussian wavepacket
  \param[in] px0 Momentum of the Gaussian wavepacket
  \param[in] dx0 Spread (distribution width) of the spatial component of the Gaussian wavepacket
  \param[in] init_state Index of the electronic state on which the wavepacket is initialized

   G(x) = [ (1/(2.0*pi*dx0^2))^(1/4) ] * exp(-((x-x0)/(2*dx0))^2 + i*(x-x0)*px0)
*/

  int sz = x0.size();
//  init_gauss_1D(PSI, *X, x0[0], px0[0], dx0[0], nstates, init_state[0], weights[0]);

  // PSI should be already allocated by this point
  for(int i=0; i<sz; i++){
    add_gauss_1D(PSI, *X, x0[i], px0[i], dx0[i], nstates, init_state[i], weights[i]);
  }
   

  // PSI(r)->PSI(k)=reciPSI
  ft_1D(PSI,reciPSI,1,xmin,kxmin,dx);

  cout<<"Wavefunction is initialized\n";

}


void Wfcgrid::init_wfc_1D_HO(vector<int>& init_state, vector<int>& nu, 
                             vector<complex<double> >& weights,
                             vector<double>& x0, vector<double>& px0, vector<double>& alpha ){
/**
  \brief Initialize 1D wavefunction as a superposition of moving Harmonic oscillator eigenstates

  \param[in] init_state Index of the electronic state to which the wavepacket is added
  \param[in] nu Vibrational quantum number of the HO basis function
  \param[in] weight The weight with which the basis function is added to the grid
  \param[in] x0 Position of the center of the HO basis function
  \param[in] px0 Momentum of the HO basis wavepacket 
  \param[in] alpha The parameter related to the reference HO Hamiltonian:  alpha = sqrt(k*mu/hbar^2)
             Where:  H = -hbar^2 /(2*mu) d^2/dx^2  + 1/2 * k * (x-x_)^2 

*/

  int sz = x0.size();
  cout<<"Initialize the wfc as a superposition of "<<sz<<" Harmonic Oscillators\n";

  // PSI should be already allocated by this point
  for(int i=0; i<sz; i++){
    cout<<"===== Oscillator "<<i<<" ======\n";
    cout<<"add it to the electronic state = "<<init_state[i]<<endl;
    cout<<"vibrational quantum number, nu = "<<nu[i]<<endl;
    cout<<"add it with the weight, w = "<<weights[i]<<endl;
    cout<<"coordinate, x0 = "<<x0[i]<<endl;
    cout<<"momentum, px0 = "<<px0[i]<<endl;
    cout<<"exponent, alpha = "<<alpha[i]<<endl;

    add_ho_1D(PSI, *X, nu[i], x0[i], px0[i], weights[i], init_state[i], alpha[i]);
  }

  normalize_wfc_1D();
   
  // PSI(r)->PSI(k)=reciPSI
  ft_1D(PSI,reciPSI,1,xmin,kxmin,dx);

  cout<<"Wavefunction is initialized\n";

}





void Wfcgrid::init_wfc_1D_ARB(bp::object py_funct, bp::object params){
/**
  \brief Initialize a 1D wavefunction according the external Python function

  \param[in] py_funct - the name of the Python-defined function. Expectations is that the 
  function returns a complex-values result for different states
  \param[in] params - parameters of that function

*/
  
  CMATRIX res(nstates, 1);

  for(int nx=0; nx<Nx; nx++){
    res = bp::extract< CMATRIX >( py_funct(real(X->M[nx]), params) );
    for(int i=0; i<nstates; i++){  PSI[i].M[nx] = res.get(nx);  }
  }
   
  // PSI(r)->PSI(k)=reciPSI
  ft_1D(PSI,reciPSI,1,xmin,kxmin,dx);

  cout<<"Wavefunction is initialized\n";

}



void Wfcgrid::init_wfc_2D(double x0, double y0, double px0, double py0, double dx0, double dy0, int init_state){
/**
  \brief Initialize 2D wavefunction - taken as moving Gaussian wavepacket
  \param[in] x0 Position of the center of the Gaussian wavepacket in X dimension
  \param[in] y0 Position of the center of the Gaussian wavepacket in Y dimension
  \param[in] px0 Momentum of the Gaussian wavepacket in X dimension
  \param[in] py0 Momentum of the Gaussian wavepacket in Y dimension
  \param[in] dx0 Spread (distribution width) of the spatial component of the Gaussian wavepacket in X dimension
  \param[in] dy0 Spread (distribution width) of the spatial component of the Gaussian wavepacket in Y dimension
  \param[in] init_state Index of the electronic state on which the wavepacket is initialized

  G(x,y) = G(x)*G(y), 
  where
  G(x) = [ (1/(2.0*pi*dx0^2))^(1/4) ] * exp(-((x-x0)/(2*dx0))^2 + i*(x-x0)*px0)
  G(y) = [ (1/(2.0*pi*dy0^2))^(1/4) ] * exp(-((y-y0)/(2*dy0))^2 + i*(y-y0)*py0)

*/

  init_gauss_2D(PSI, *X, x0, px0, dx0,    *Y, y0, py0, dy0,  nstates, init_state);

  // PSI(r)->PSI(k)=reciPSI
  ft_2D(PSI,reciPSI,1,xmin,ymin,kxmin,kymin,dx,dy);


  cout<<"Wavefunction is initialized\n";


}// init_wavefunction


void Wfcgrid::print_wfc_1D(std::string prefix, int snap, int state){
/**
  \brief Print 1D wavefunction into a file
  \param[in] prefix The prefix of the filenames to which the wfc will be printed out
  \param[in] snap This is the integer indexing the printed out wavefunction - e.g. convenient for printing wfc at different times
  \param[in] state The wavefunction projection on this state will be printed out

  for 1D profile on XY plane
*/

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
/**
  \brief Print 2D wavefunction into a file
  \param[in] prefix The prefix of the filenames to which the wfc will be printed out
  \param[in] snap This is the integer indexing the printed out wavefunction - e.g. convenient for printing wfc at different times
  \param[in] state The wavefunction projection on this state will be printed out

  for 2D profile on XY plane
*/

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
/**
  \biref Print 1D wavefunction (in k-representation, reciprocal space) into a file
  \param[in] prefix The prefix of the filenames to which the wfc will be printed out
  \param[in] snap This is the integer indexing the printed out wavefunction - e.g. convenient for printing wfc at different times
  \param[in] state The wavefunction projection on this state will be printed out

  for 1D profile on XY plane
*/

  std::string filename, snaps, states;
  stringstream ss(stringstream::in | stringstream::out);
  stringstream ss1(stringstream::in | stringstream::out);

  ss  << snap;   ss  >> snaps;
  ss1 << state;  ss1 >> states;


  filename = prefix+".state"+states+".frame"+snaps;
  ofstream out(filename.c_str(),ios::out);


  for(int nx=0;nx<Nx;nx++){

    out<<real(2.0*M_PI*Kx->M[nx])<<"  "<<real(std::conj(reciPSI[state].M[nx])*reciPSI[state].M[nx])<<endl;

  }// for nx

  out.close();

}// print_wfc_1D


void Wfcgrid::print_reci_wfc_2D(std::string prefix, int snap, int state){
/**
  \brief Prints 2D wavefunction (in k-representation, reciprocal space) into a file
  \param[in] prefix The prefix of the filenames to which the wfc will be printed out
  \param[in] snap This is the integer indexing the printed out wavefunction - e.g. convenient for printing wfc at different times
  \param[in] state The wavefunction projection on this state will be printed out
  
  for 2D profile on XY plane 
*/

  std::string filename, snaps, states;
  stringstream ss(stringstream::in | stringstream::out);
  stringstream ss1(stringstream::in | stringstream::out);

  ss  << snap;   ss  >> snaps;
  ss1 << state;  ss1 >> states;


  filename = prefix+".state"+states+".frame"+snaps;
  ofstream out(filename.c_str(),ios::out);


  for(int nx=0;nx<Nx;nx++){
    for(int ny=0;ny<Ny;ny++){

      out<<real(2.0*M_PI*Kx->M[nx])<<"  "<<real(2.0*M_PI*Ky->M[ny])<<"  "<<real(std::conj(reciPSI[state].M[nx*Ny+ny])*reciPSI[state].M[nx*Ny+ny])<<endl;

    }// for ny
    out<<"\n";
  }// for nx

  out.close();
}


void Wfcgrid::print_complex_matrix_1D(CMATRIX& CM, std::string filename){
/**
  \brief An auxiliary function for printing a complex matrix into a file
  \param[in] CM The complex matrix (1 by Nx or Nx by 1 to be printed out)
  \param[in] filename The name of the file to which the matrix will be printed out
*/

  ofstream out(filename.c_str(),ios::out);

  for(int nx=0;nx<Nx;nx++){
    out<<real(X->M[nx])<<"  "<<real(2.0*M_PI*Kx->M[nx])<<"  "<<CM.M[nx].real()<<"   "<<CM.M[nx].imag()<<endl;
  }// for nx

  out.close();

}

void Wfcgrid::print_ham_1D(std::string prefix, int i, int j){
/**
  \brief Printing of the Hamiltonian of 1D wavefunctions
  \param[in] prefix The prefix of the file to which the info will be printed out
  \param[in] i The index of one of the states
  \param[in] j The index of the other state

  Printing <i|H|j>
*/

  print_complex_matrix_1D(H[i][j],prefix);
}

void Wfcgrid::print_expH_1D(std::string prefix, int i, int j){
/**
  \brief Printing of the exponential of the Hamiltonian of 1D wavefunctions (component of propagator)
  \param[in] prefix The prefix of the file to which the info will be printed out
  \param[in] i The index of one of the states
  \param[in] j The index of the other state

  Printing exp(-i*dt/hbar * <i|H|j>)
*/

  print_complex_matrix_1D(expH[i][j],prefix);
}

void Wfcgrid::print_expK_1D(std::string prefix, int i){
/**
  \brief Printing of the exponential of the kinetic energy operator (in reciprocal space) for 1D wavefunctions
  \param[in] prefix The prefix of the file to which the info will be printed out
  \param[in] i The index of the state - the operator is diagonal, so only need this one index

  Printing exp(-i*dt/hbar * <i|T|i>)
*/

  print_complex_matrix_1D(expK[i],prefix);
}




double Wfcgrid::print_populations_1D(string filename,int snap){
/**
  \brief Print population of 1D wavefunction into a file
  \param[in] filename The name of the file to which the wavefunction will be printed out
  \param[in] snap This is the integer indexing the printed out wavefunction - e.g. convenient for printing wfc at different times
  
  Populations of all states will be printed out
  for 1D profile on XY plane
*/

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


double Wfcgrid::print_populations_2D(string filename,int snap){
/**
  \brief Print population of 1D wavefunction into a file
  \param[in] filename The name of the file to which the wavefunction will be printed out
  \param[in] snap This is the integer indexing the printed out wavefunction - e.g. convenient for printing wfc at different times

  Populations of all states will be printed out
  for 2D profile on XY plane
*/


  vector<double> Pop(nstates,0.0);
  double Pop_tot = 0.0;
  double Pop_tot_active = 0.0;

  for(int nx=0;nx<Nx;nx++){
    for(int ny=0;ny<Ny;ny++){
      for(int nst=0;nst<nstates;nst++){

        double pii = real(std::conj(PSI[nst].M[nx*Ny+ny])*PSI[nst].M[nx*Ny+ny]);

        Pop[nst] += dx*dy*pii;  // population on state nst
        Pop_tot += dx*dy*pii;   // total population

      }// for nst
    }// for ny
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
/**
  \brief Compute the population flux in 1D case
  \param[in] xf The point at which flux is computed
  \param[out] res The collector of the fluxes (for each electronic state projection of the wfc)
  \param[in] m0 The effective mass of the quantum particle (DOF)
*/

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



double Wfcgrid::get_x_1D(){
/**
  Compute expectation value of coordinate of a 1D wavefunction: <psi|x|psi> / <psi|psi>
*/

  double res; res = 0.0;
  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){        
    for(int nst1=0;nst1<nstates;nst1++){        

      for(int nx=0;nx<Nx;nx++){
   
        double d = real(std::conj(PSI[nst].M[nx]) * PSI[nst1].M[nx] );

        norm += d; 
        res += d * real(X->M[nx]);

      }// for nx

    }// for nst1
  }// for nst

  res = res / norm;

  return res;

}


double Wfcgrid::get_pow_x_1D(double n){
/**
  Compute expectation value of a power of coordinate x^n a 1D wavefunction: <psi|x^n|psi> / <psi|psi>
*/

  double res; res = 0.0;
  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){        
    for(int nst1=0;nst1<nstates;nst1++){        

      for(int nx=0;nx<Nx;nx++){
   
        double d = real(std::conj(PSI[nst].M[nx]) * PSI[nst1].M[nx] );

        norm += d; 
        res += d * pow(real(X->M[nx]), n);

      }// for nx

    }// for nst1
  }// for nst

  res = res / norm;

  return res;

}




double Wfcgrid::get_px_1D(){
/**
  Compute the expectation value of momentum of a 1D wavefunction: <psi|-i*hbar*d/dx |psi> / <psi|psi>
*/

  double res; res = 0.0;
  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){        
    for(int nst1=0;nst1<nstates;nst1++){        

      for(int nx=0;nx<Nx;nx++){

        norm += real(std::conj(PSI[nst].M[nx]) * PSI[nst1].M[nx] );  
        res += real(Kx->M[nx]) * real( std::conj(reciPSI[nst].M[nx]) * reciPSI[nst1].M[nx] );

      }// for nx

    }// for nst1
  }// for nst
 
  norm *= dx;
  res *= (2.0*M_PI/((double)Nx*dx));
  res = res / norm;

  return res;

}



double Wfcgrid::get_pow_px_1D(int n){
/**
  Compute the expectation value of power of momentum of a 1D wavefunction: <psi|(-i*hbar*d/dx)^n|psi> / <psi|psi>
*/

  double res; res = 0.0;
  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){        
    for(int nst1=0;nst1<nstates;nst1++){        

      for(int nx=0;nx<Nx;nx++){

        norm += real(std::conj(PSI[nst].M[nx]) * PSI[nst1].M[nx] );  
        res += pow(real(Kx->M[nx]), n) * real( std::conj(reciPSI[nst].M[nx]) * reciPSI[nst1].M[nx] );

      }// for nx

    }// for nst1
  }// for nst
 
  norm *= dx;
  res *= (pow(2.0*M_PI, n)/((double)Nx*dx));
  res = res / norm;

  return res;

}






double Wfcgrid::get_x_2D(){
/**
  Compute expectation value of x coordinate of a 2D wavefunction: <psi|x|psi> / <psi|psi>
*/

  double res; res = 0.0;
  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){        
    for(int nst1=0;nst1<nstates;nst1++){        

      for(int nx=0;nx<Nx;nx++){
        for(int ny=0;ny<Ny;ny++){
   
          double d = real(std::conj(PSI[nst].M[nx*Ny+ny]) * PSI[nst1].M[nx*Ny+ny] );

          norm += d; 
          res += d * real(X->M[nx]);
      
        }// for ny
      }// for nx

    }// for nst1
  }// for nst

  res = res / norm;

  return res;

}


double Wfcgrid::get_y_2D(){
/**
  Compute expectation value of y coordinate of a 2D wavefunction: <psi|y|psi> / <psi|psi>
*/

  double res; res = 0.0;
  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){        
    for(int nst1=0;nst1<nstates;nst1++){        

      for(int nx=0;nx<Nx;nx++){
        for(int ny=0;ny<Ny;ny++){
   
          double d = real(std::conj(PSI[nst].M[nx*Ny+ny]) * PSI[nst1].M[nx*Ny+ny] );

          norm += d; 
          res += d * real(Y->M[ny]);
      
        }// for ny
      }// for nx

    }// for nst1
  }// for nst

  res = res / norm;

  return res;

}



double Wfcgrid::get_px_2D(){
/**
  Compute the expectation value of x component of momentum of a 2D wavefunction: <psi|-i*hbar*d/dx |psi> / <psi|psi>
*/

  double res; res = 0.0;
  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){        
    for(int nst1=0;nst1<nstates;nst1++){        

      for(int nx=0;nx<Nx;nx++){
        for(int ny=0;nx<Ny;ny++){

          norm += real(std::conj(PSI[nst].M[nx*Ny+ny]) * PSI[nst1].M[nx*Ny+ny] );  
          res += real(Kx->M[nx]) * real( std::conj(reciPSI[nst].M[nx*Ny+ny]) * reciPSI[nst1].M[nx*Ny+ny] );

        }// for ny
      }// for nx

    }// for nst1
  }// for nst
 
  norm *= (dx * dy);
  res *= (2.0*M_PI/((double)Nx*dx));
  res *= (2.0*M_PI/((double)Ny*dy));

  res = res / norm;

  return res;


}



double Wfcgrid::get_py_2D(){
/**
  Compute the expectation value of y component of momentum of a 2D wavefunction: <psi|-i*hbar*d/dy |psi> / <psi|psi>
*/

  double res; res = 0.0;
  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){        
    for(int nst1=0;nst1<nstates;nst1++){        

      for(int nx=0;nx<Nx;nx++){
        for(int ny=0;nx<Ny;ny++){

          norm += real(std::conj(PSI[nst].M[nx*Ny+ny]) * PSI[nst1].M[nx*Ny+ny] );  
          res += real(Ky->M[ny]) * real( std::conj(reciPSI[nst].M[nx*Ny+ny]) * reciPSI[nst1].M[nx*Ny+ny] );

        }// for ny
      }// for nx

    }// for nst1
  }// for nst
 
  norm *= (dx * dy);
  res *= (2.0*M_PI/((double)Nx*dx));
  res *= (2.0*M_PI/((double)Ny*dy));

  res = res / norm;

  return res;


}


double Wfcgrid::norm_1D(){
/**
  Compute the norm of the 1D wavefunction: <psi|psi>
*/

  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){        
    for(int nst1=0;nst1<nstates;nst1++){        

      for(int nx=0;nx<Nx;nx++){

        norm += real(std::conj(PSI[nst].M[nx]) * PSI[nst1].M[nx] );  

      }// for nx

    }// for nst1
  }// for nst

  norm *= dx;

  return norm;


}


void Wfcgrid::normalize_wfc_1D(){
/**
 *   Normalize the wavefunction
 *   */

  double sqrt_norm = sqrt( norm_1D() );

  for(int nst=0;nst<nstates;nst++){
    for(int nx=0;nx<Nx;nx++){

      PSI[nst].M[nx] /= sqrt_norm;

      }// for nx
  }// for nst

}


double Wfcgrid::e_kin_1D(double m0){
/**
  Compute kinetic energy for 1D wavefunction: <psi|T|psi> / <psi|psi>
*/
  double res; res = 0.0;
  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){        
    for(int nst1=0;nst1<nstates;nst1++){        

      for(int nx=0;nx<Nx;nx++){

        norm += real(std::conj(PSI[nst].M[nx]) * PSI[nst1].M[nx] );  
        res += real(Kx->M[nx]) * real(Kx->M[nx]) * real(std::conj(reciPSI[nst].M[nx]) * reciPSI[nst1].M[nx] );

      }// for nx

    }// for nst1
  }// for nst

  norm *= dx;
  res *= (1.0/((double)Nx*dx));
  res *= (2.0*M_PI*M_PI/m0);

  res = res/norm;

  return res;


}// e_kin_1D


double Wfcgrid::e_pot_1D(){
/**
  Compute potential energy for 1D wavefunction: <psi|V|psi> / <psi|psi>
*/
  double res; res = 0.0;
  double norm; norm = 0.0;


  for(int nst=0;nst<nstates;nst++){        
    for(int nst1=0;nst1<nstates;nst1++){               

      for(int nx=0;nx<Nx;nx++){

        norm += real(std::conj(PSI[nst].M[nx]) * PSI[nst1].M[nx] ); 
        res += real(std::conj(PSI[nst].M[nx]) * H[nst][nst1].M[nx] * PSI[nst1].M[nx] ); 

      }// for nx

    }// for nst1
  }// for nst


  res = res / norm;

  return res;


}// e_pot_1D



double Wfcgrid::e_tot_1D(double m0){
/**
  Compute total energy for 1D wavefunction: <psi|T+V|psi> / <psi|psi>
*/
  double res = e_kin_1D(m0) + e_pot_1D();
  
  return res;

}// e_tot_1D




double Wfcgrid::norm_2D(){
/**
  Compute the norm for 2D wavefunction: <psi|psi> 
  
*/

  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){
    for(int nst1=0;nst1<nstates;nst1++){     

      for(int nx=0;nx<Nx;nx++){
        for(int ny=0;ny<Ny;ny++){
          
          norm += real(std::conj(PSI[nst].M[nx*Ny+ny])*PSI[nst1].M[nx*Ny+ny]);

        }// for ny
      }// for nx

    }// for nst1
  }// for nst 

  norm *= (dx * dy);

  return norm;


}// e_norm_2D


double Wfcgrid::e_kin_2D(double m1, double m2){
/**
  Compute kinetic energy for 2D wavefunction: <psi|T|psi> / <psi|psi>
  
*/


  double res_x; res_x = 0.0;
  double res_y; res_y = 0.0;
  double norm; norm = 0.0;
  double tmp;

  for(int nst=0;nst<nstates;nst++){
    for(int nst1=0;nst1<nstates;nst1++){     

      for(int nx=0;nx<Nx;nx++){
        for(int ny=0;ny<Ny;ny++){
          
          norm += real(std::conj(PSI[nst].M[nx*Ny+ny])*PSI[nst1].M[nx*Ny+ny]);

          double e1 = real(Kx->M[nx]) * real(Kx->M[nx]);
          double e2 = real(Ky->M[ny]) * real(Ky->M[ny]);
          tmp = real( std::conj(reciPSI[nst].M[nx*Ny+ny]) * reciPSI[nst1].M[nx*Ny+ny]);

          res_x += e1 * tmp; 
          res_y += e2 * tmp; 


        }// for ny
      }// for nx

    }// for nst1
  }// for nst 

  double res = res_x / m1 + res_y / m2;

  norm *= (dx * dy);
  res *= (1.0/((double)Nx*dx));
  res *= (1.0/((double)Ny*dy));
  res *= (2.0*M_PI*M_PI);

  res = res / norm;

  return res;


}// e_kin_2D



double Wfcgrid::e_pot_2D(){
/**
  Compute potential energy for 2D wavefunction: <psi|V|psi> / <psi|psi>
  
*/


  double res; res = 0.0;
  double norm; norm = 0.0;

  for(int nst=0;nst<nstates;nst++){
    for(int nst1=0;nst1<nstates;nst1++){     

      for(int nx=0;nx<Nx;nx++){
        for(int ny=0;ny<Ny;ny++){
          
          norm += real(std::conj(PSI[nst].M[nx*Ny+ny])*PSI[nst1].M[nx*Ny+ny]);
          res += real(std::conj(PSI[nst].M[nx*Ny+ny]) * H[nst][nst1].M[nx*Ny+ny] * PSI[nst1].M[nx*Ny+ny] ); 

        }// for ny
      }// for nx

    }// for nst1
  }// for nst 


  res = res / norm;

  return res;


}// e_pot_2D


double Wfcgrid::e_tot_2D(double m1, double m2){
/**
  Compute total energy for 2D wavefunction: <psi|T+V|psi> / <psi|psi>
*/
  double res = e_kin_2D(m1, m2) + e_pot_2D();
  
  return res;

}// e_tot_2D






}// namespace libwfcgrid
}// namespace libdyn
}// liblibra

