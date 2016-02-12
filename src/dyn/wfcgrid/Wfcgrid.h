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
/**
  \file Wfcgrid.h
  \brief The file describes a Wfcgrid class for exact numerical solution of 
  time-dependent Schrodinger equation on the grid. 1D and 2D grids are available
    
*/

#ifndef WFCGRID_H
#define WFCGRID_H


#include "../../mmath/libmmath.h"
#include "../../hamiltonian/libhamiltonian.h"
#include "Grid_functions.h"
using namespace libmmath;
using namespace libhamiltonian;

/// libdyn namespace
namespace libdyn{

/// libwfcgrid namespace
namespace libwfcgrid{

// 1D and 2D grid wavefunction
class Wfcgrid{
/**
  \brief The Wfcgrid class

  The class that represents wavefunction numerically on the 1D or 2D grids

*/


  void init_numbers(double minx_, double maxx_, double dx_, int nstates_);
  void init_numbers(double minx_, double maxx_, double dx_, double miny_, double maxy_, double dy_, int nstates_);
  void allocate_1D();
  void allocate_2D();
  void init_grid_1D();
  void init_grid_2D();

  void print_complex_matrix_1D(CMATRIX& CM, std::string filename);

public:

  // 2-D Grid
  int nstates;         ///< number of electronic states
  int Nx, Ny;          ///< grid size along X and Y coordinates
  double xmin, ymin;   ///< grid boundaries in x and y directions in real space 
  double xmax, ymax;   ///< grid boundaries in x and y directions in real space
  double dx, dy;       ///< grid point spacing in x and y dimensions
  double kxmin, kymin; ///< minimal values of grid in the reciprocal (momentum)  space
  CMATRIX* X;          ///< grid of r-points 1 x Nx
  CMATRIX* Y;          ///< grid of r-points 1 x Ny
  CMATRIX* Kx;         ///< grid of k-points 1 x Nx
  CMATRIX* Ky;         ///< grid of k-points 1 x Ny




  // True quantum properties
  vector<CMATRIX> PSI;       ///< wavefunction:  nstates x Nx x Ny
  vector<CMATRIX> reciPSI;   ///< same as PSI but in Fourier (reciprocal) space with 1.0*kmin - nstates x Nx x Ny
  vector<CMATRIX> DtreciPSI; ///< d/dt * reciPSI - time derivative of PSI in reciprocal space  - nstates x Nx x Ny
  vector<CMATRIX> DxPSI;     ///< d/dx * PSI in real space (up to a factor) - nstates x Nx x Ny
  vector<CMATRIX> DyPSI;     ///< d/dx * PSI in real space (up to a factor)  - nstates x Nx x Ny
  vector<CMATRIX> KxreciPSI; ///< Kx * reciPSI  -  d/dx * PSI in k space (up to a factor) - nstates x Nx x Ny
  vector<CMATRIX> KyreciPSI; ///< Ky * reciPSI  -  d/dy * PSI in k space (up to a factor) - nstates x Nx x Ny

  // Hamiltonian and Propagator
  vector< vector<CMATRIX> > H;  ///< PES - nstates x nstates x Nx x Ny - full Hamiltonian 
  vector< vector<CMATRIX> > Dx; ///< for each i,j pair - a separate x projection for all r-points of the grid - nstates x nstates x Nx x Ny
  vector< vector<CMATRIX> > Dy; ///< for each i,j pair - a separate y projection for all r-points of the grid - nstates x nstates x Nx x Ny  
  vector< vector<CMATRIX> > expH; ///< nstates x nstates x Nx x Ny
  vector<CMATRIX> expK;         ///< nstates x Nx x Ny



  // Constructor : 1D and 2D
  Wfcgrid(double minx_, double maxx_, double dx_, int nstates_); ///< constructor for 1D wavefunction
  Wfcgrid(double minx_, double maxx_, double dx_, double miny_, double maxy_, double dy_, int nstates_); ///< constructor for 2D wavefunction

  // Populate wfc
  void init_wfc_1D(double x0, double px0, double dx, int init_state); ///< initialization of 1D wavefunction
  void init_wfc_2D(double x0, double y0, double px0, double py0, double dx, double dy, int init_state); ///< initialization of 2D wavefunction

  // Print 1D and 2D wavefunctions to file
  void print_wfc_1D(std::string prefix, int snap, int state);
  void print_wfc_2D(std::string prefix, int snap, int state);
  void print_reci_wfc_1D(std::string prefix, int snap, int state);
  void print_ham_1D(std::string prefix, int i, int j);
  void print_expH_1D(std::string prefix, int i, int j);
  void print_expK_1D(std::string prefix, int i);

 
  // Print state-resolved populations 
  double print_populations_1D(string filename,int snap);
  double print_populations_2D(string filename,int snap);

  // Flux
  void flux_1D(double xf,vector<double>& res, double m0);

  // Energy
  double e_pot_1D();
  double e_kin_1D(double m0);
  double e_tot_1D(double m0);


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
