#ifndef WFCGRID_H
#define WFCGRID_H


#include "../../mmath/libmmath.h"
#include "../../hamiltonian/libhamiltonian.h"
#include "Grid_functions.h"
using namespace libmmath;
using namespace libhamiltonian;

namespace libdyn{
namespace libwfcgrid{

// 1D and 2D grid wavefunction
class Wfcgrid{

  void init_numbers(double minx_, double maxx_, double dx_, int nstates_);
  void init_numbers(double minx_, double maxx_, double dx_, double miny_, double maxy_, double dy_, int nstates_);
  void allocate_1D();
  void allocate_2D();
  void init_grid_1D();
  void init_grid_2D();

  void print_complex_matrix_1D(CMATRIX& CM, std::string filename);

public:

  // 2-D Grid
  int nstates;         // number of electronic states
  int Nx, Ny;          // grid size along X and Y coordinates
  double xmin, ymin;   // grid boundaries in x and y directions in real space 
  double xmax, ymax;   // grid boundaries in x and y directions in real space
  double dx, dy;       // grid point spacing in x and y dimensions
  double kxmin, kymin; // 
  CMATRIX* X;          // grid of r-points 1 x Nx
  CMATRIX* Y;          // grid of r-points 1 x Ny
  CMATRIX* Kx;         // grid of k-points  1 x Nx
  CMATRIX* Ky;         // grid of k-points  1 x Ny




  // True quantum properties
  vector<CMATRIX> PSI;       // wavefunction:  nstates x Nx x Ny
  vector<CMATRIX> reciPSI;   // same as PSI but in Fourier (reciprocal) space with 1.0*kmin - nstates x Nx x Ny
  vector<CMATRIX> DtreciPSI; // d/dt * reciPSI - time derivative of PSI in reciprocal space  - nstates x Nx x Ny
  vector<CMATRIX> DxPSI;     // d/dx * PSI in real space (up to a factor) - nstates x Nx x Ny
  vector<CMATRIX> DyPSI;     // d/dx * PSI in real space (up to a factor)  - nstates x Nx x Ny
  vector<CMATRIX> KxreciPSI; // Kx * reciPSI  -  d/dx * PSI in k space (up to a factor) - nstates x Nx x Ny
  vector<CMATRIX> KyreciPSI; // Ky * reciPSI  -  d/dy * PSI in k space (up to a factor) - nstates x Nx x Ny

  // Hamiltonian and Propagator
  vector< vector<CMATRIX> > H;  // PES - nstates x nstates x Nx x Ny - full Hamiltonian 
  vector< vector<CMATRIX> > Dx; // for each i,j pair - a separate x projection for all r-points of the grid - nstates x nstates x Nx x Ny
  vector< vector<CMATRIX> > Dy; // for each i,j pair - a separate y projection for all r-points of the grid - nstates x nstates x Nx x Ny  
  vector< vector<CMATRIX> > expH; // nstates x nstates x Nx x Ny
  vector<CMATRIX> expK;         // nstates x Nx x Ny



  // Constructor : 1D and 2D
  Wfcgrid(double minx_, double maxx_, double dx_, int nstates_); // 1D
  Wfcgrid(double minx_, double maxx_, double dx_, double miny_, double maxy_, double dy_, int nstates_); // 2D

  // Populate wfc
  void init_wfc_1D(double x0, double px0, double dx, int init_state); // 1D
  void init_wfc_2D(double x0, double y0, double px0, double py0, double dx, double dy, int init_state); // 2D

  // Print 1D and 2D wavefunctions to file
  void print_wfc_1D(std::string prefix, int snap, int state);
  void print_wfc_2D(std::string prefix, int snap, int state);
  void print_reci_wfc_1D(std::string prefix, int snap, int state);
  void print_ham_1D(std::string prefix, int i, int j);
  void print_expH_1D(std::string prefix, int i, int j);
  void print_expK_1D(std::string prefix, int i);

 
  // Print state-resolved populations 
  double print_populations_1D(string filename,int snap);

  // Flux
  void flux_1D(double xf,vector<double>& res, double m0);


  //--------------- in Wfcgrid_Dynamics1 ------------------

  void update_potential_1D(Hamiltonian& ham);
  void update_potential_2D(Hamiltonian& ham);

  void update_propagator_1D(double dt,double m0);
  void update_propagator_2D(double dt,double m0);

  void update_propagator_K_1D(double dt,double m0);
  void update_propagator_K_2D(double dt,double m0);


  void propagate_exact_1D(int Nmts);
  void propagate_exact_2D(int Nmts);

  void absorb_1D(double dL,vector<double>& Pops_l,vector<double>& Pops_r);
  boost::python::list absorb_1D(double dL);


}; //  class Wfcgrid



}// namespace libwfcgrid
}// namespace libdyn

#endif  // WFCGRID_H
