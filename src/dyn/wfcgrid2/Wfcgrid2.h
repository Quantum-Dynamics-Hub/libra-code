/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Wfcgrid2.h
  \brief The file describes a Wfcgrid2 class for exact numerical solution of 
  time-dependent Schrodinger equation on the arbitrary-dimensional grids.
    
*/

#ifndef WFCGRID2_H
#define WFCGRID2_H

#include "../wfcgrid/libwfcgrid.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace librandom;
using namespace libspecialfunctions;
//using namespace libhamiltonian;

/// libdyn namespace
namespace libdyn{

using namespace libwfcgrid;

/// libwfcgrid namespace
namespace libwfcgrid2{


///=============== In the Wfcgrid2_initialize.cpp ====================
CMATRIX Gaussian(vector<double>& x, vector<double>& x0, vector<double>& px0, vector<double>& dx, int init_state, int nstates, complex<double> scl);
CMATRIX HO(vector<double>& x, vector<double>& x0, vector<double>& px0, 
           vector<double>& alpha, int init_state, int nstates, vector<int>& nu, complex<double> scl);


///=============== In the Wfcgrid2_ColbertMiller.cpp ==================
int points_on_same_line(int idof, vector<int>& pt1, vector<int>& pt2);

void dvr0(MATRIX& T, int npts, double rmin, double rmax, double mass);
void dvr1(MATRIX& T, int npts, double dr, double mass);
void dvr2(MATRIX& T, int npts, double dr, double mass);

//vector<CMATRIX> expV(vector<CMATRIX>& V, complex<double> dt);
//vector<CMATRIX> expV_diag(vector<CMATRIX>& V, complex<double> dt);

// N-D grid wavefunction
class Wfcgrid2{
/**
  \brief The Wfcgrid2 class

  The class that represents wavefunction numerically on the arbitrary-dimensional grids

*/

  ///=============== In the Wfcgrid2.cpp ====================
  void init_numbers(vector<double>& rmin_, vector<double>& rmax_, vector<double>& dr_, int nstates_);
  void allocate();
  void init_grids();
  void compute_mapping();

  //void print_complex_matrix_1D(CMATRIX& CM, std::string filename);

  /// For direct integrator:
  vector<CMATRIX> PSI_dia_past;      ///< wavefunction at time t-dt:    list of Npts matrices  nstates x 1
  vector<CMATRIX> PSI_adi_past;      ///< wavefunction at time t-dt:    list of Npts matrices  nstates x 1


public:

  // ndof-D Grid
  int nstates;            ///< number of electronic states
  int ndof;               ///< number of nuclear DOFs


  ///< Grids
  int Npts;               ///< the total number of grid points
  vector<int> npts;       ///< number of grid points in each DOF dimension
  vector<double> rmin;    ///< minimal grid boundary in each DOF dimension in real space 
  vector<double> rmax;    ///< maximal grid boundary in each DOF dimension in real space 
  vector<double> dr;      ///< real space grid points spacing for each DOF dimension 
  vector<double> kmin;    ///< minimal grid boundary in each DOF dimension in reciprocal space 
  vector<double> dk;      ///< reciprocal space grid points spacing for each DOF dimension 

  vector< MATRIX* > rgrid;  ///< rgrid[dof].get(ipt) - coordinate of the real-space point ipt along the dof dimension
  vector< MATRIX* > kgrid;  ///< kgrid[dof].get(ipt) - coordinate of the reciptocal-space point ipt along the dof dimension
  /**
  Note the relationship:  p[dof][ipt] = 2*M_PI*kgrid[dof][ipt] , the kgrid stored consists of the wavevectors rather than momenta
  */

  ///< Grid mapping : the wavefunctions below are stored in a consecutive order (see below)
  /// to get from a single integer (which is just an order of the point in a real or reciprocal space) from 
  /// the indices of the point on the 1D grid in each dimensions, we use the mapping below:
  /// e.g. igmap[1] = [0, 1, 0, 0] means that the second (index 1) entry in the PSI array below corresponds to
  /// a grid point that is first (lower boundary) in dimensions 0, 2, and 3, but is second (index 1) in the 
  /// dimension 1. Same for the reciprocal space
  vector< vector<int> > gmap;  


  ///=============== In the Wfcgrid2.cpp ====================
  ///< Inverse grid points mapping     
  int imap(vector<int>& point);

 
  /// Quantum properties
  /// Wavefunction is the diabatic basis
  vector<CMATRIX> PSI_dia;       ///< wavefunction:  list of Npts matrices  nstates x 1
  vector<CMATRIX> reciPSI_dia;   ///< same as PSI but in Fourier (reciprocal) space with 1.0*kmin 

  /// Wavefunction is the adiabatic basis
  vector<CMATRIX> PSI_adi;       ///< wavefunction:  list of Npts matrices  nstates x 1
  vector<CMATRIX> reciPSI_adi;   ///< same as PSI but in Fourier (reciprocal) space with 1.0*kmin 

  /// Wavefunction spatial derivatives
  vector< vector<CMATRIX> > nabla_PSI_dia;  /// d/dR PSI_dia in real space; nablaPSI_dia[dof][ipt].get(istate, 0)   dPSI_dia/dR_dof ( r[ipt]) - nstates x 1 matrix
  vector< vector<CMATRIX> > nabla_PSI_adi;  /// d/dR PSI_adi in real space; nablaPSI_adi[dof][ipt].get(istate, 0)   dPSI_adi/dR_dof ( r[ipt]) - nstates x 1 matrix
  vector< vector<CMATRIX> > nabla_reciPSI_dia; /// d/dR PSI_dia in k-space - same structure as nabla_PSI_dia;
  vector< vector<CMATRIX> > nabla_reciPSI_adi; /// d/dR PSI_adi in k-space - same structure as nabla_PSI_adi;

/*
  vector<CMATRIX> DtreciPSI; ///< d/dt * reciPSI - time derivative of PSI in reciprocal space 
*/

  /// Hamiltonian and Propagator
  vector<CMATRIX> Hdia;      ///<  diabatic Hamiltoninans for all the Npts points
  vector<CMATRIX> Hadi;      ///<  adiabatic Hamiltoninans for all the Npts points
  vector<CMATRIX> Vcomplex;  ///<  complex absorbing potential that will be added to the real-space propagator
  vector< vector<CMATRIX> > NAC1;  ///<  1-st order NACS: NAC1[ipt][alpha].get(i,j) = <psi_i| nabla_alpha | psi_j> 
  vector< vector<CMATRIX> > NAC2;  ///<  2-nd order NACS: NAC2[ipt][alpha].get(i,j) = <psi_i| nabla_alpha^2 | psi_j> 
  vector<CMATRIX> U;         ///<  |adi> = |dia> * U : diabatic-to-adiabatic transformation for all the Npts points
  vector<CMATRIX> expH;      ///<  exponent of the diabatic Hamiltoninans for all the Npts points
  vector<CMATRIX> expK;      ///<  exponent of the kinetik energy propagator for all the Npts points


  ///=============== In the Wfcgrid2.cpp ====================
  ///< Grid constructor
  Wfcgrid2(vector<double>& rmin_, vector<double>& rmax_, vector<double>& dr_, int nstates_); ///< constructor for n-D wavefunction


  ///=============== In the Wfcgrid2_ColbertMiller.cpp ===============
  vector<CMATRIX> T_PSI(vector<CMATRIX>& inp_psi, vector<int>& bc_type, vector<double>& mass, complex<double> scaling);

  vector<CMATRIX> T_PSI_adi(vector<int>& bc_type, vector<double>& mass, complex<double> scaling);
  vector<CMATRIX> T_PSI_dia(vector<int>& bc_type, vector<double>& mass, complex<double> scaling);

  vector<CMATRIX> expT_PSI(double dt, vector<double>& mass, int nterms);

  CMATRIX operator_T(vector<int>& bc_type, vector<double>& mass, complex<double> scaling);

  void Colbert_Miller_SOFT(CMATRIX& expT, vector<CMATRIX>& expV, int opt);

  vector<CMATRIX> nubla_PSI(int idof, vector<CMATRIX>& inp_psi, complex<double> scaling);
  vector<CMATRIX> nubla_PSI_adi(int idof, complex<double> scaling);
  vector<CMATRIX> nubla_PSI_dia(int idof, complex<double> scaling);

  vector<CMATRIX> V_PSI(vector<CMATRIX>& V, vector<CMATRIX>& inp_psi, complex<double> scaling);
  vector<CMATRIX> V_PSI_adi(complex<double> scaling);
  vector<CMATRIX> V_PSI_dia(complex<double> scaling);

  vector<CMATRIX> H_PSI_adi(vector<double>& mass);
  vector<CMATRIX> H_PSI_dia(vector<double>& mass);

  void Colbert_Miller_propagate_adi1(double dt, vector<double>& mass);
  void Colbert_Miller_propagate_adi2(double dt, vector<double>& mass);
  void Colbert_Miller_propagate_dia1(double dt, vector<double>& mass);
  void Colbert_Miller_propagate_dia2(double dt, vector<double>& mass);

  ///=============== In the Wfcgrid2_direct.cpp ====================
  void direct_allocate_tmp_vars(double rep);
  void direct_propagate_adi2(double dt, vector<double>& mass);
  void direct_propagate_adi1(double dt, vector<double>& mass);
  void direct_propagate_dia2(double dt, vector<double>& mass);
  void direct_propagate_dia1(double dt, vector<double>& mass);



  ///=============== In the Wfcgrid2_initialize.cpp ====================
  ///< Wavefunction initialization on the grid
  void add_wfc_Gau(vector<double>& x0, vector<double>& px0, vector<double>& dx0, int init_state, complex<double> weight, int rep);
  void add_wfc_HO(vector<double>& x0, vector<double>& px0, vector<double>& alpha, int init_state, vector<int>& nu, complex<double> weight, int rep);
  void add_wfc_ARB(bp::object py_funct, bp::object params, int rep);


  ///=============== In the Wfcgrid2_properties.cpp ====================  
  double norm(int rep);
  double e_kin(vector<double>& mass, int rep);
  double e_pot(int rep);
  double e_tot(vector<double>& mass, int rep);

  CMATRIX get_pow_q(int rep, int n);
  CMATRIX get_pow_p(int rep, int n);

  CMATRIX get_den_mat(int rep);

  MATRIX get_pops(int rep);
  MATRIX get_pops(int rep, vector<double>& bmin, vector<double>& bmax);

  void compute_wfc_gradients(int rep, int idof, double mass);


  ///=============== In the Wfcgrid2_SOFT.cpp ====================  
  void update_propagator_H(double dt);
  void update_propagator_K(double dt, vector<double>& mass);
  void SOFT_propagate();


  ///=============== In the Wfcgrid2_transforms.cpp ====================
  ///< Update reciprocal or real functions 
  void update_reciprocal(int rep);
  void update_real(int rep);

  ///< Normalize the wfc
  void normalize(int rep);

  ///< Reshape wfc from the internal format to something more suitable for the 1D or 2D FFTs
  void reshape_wfc_1D(int _rep, int _r_or_k, int _dir, vector<CMATRIX>& _tmp); // reshape wfc into/from the nstates x CMATRIX(Nx, 1) format
  void reshape_wfc_2D(int _rep, int _r_or_k, int _dir, vector<CMATRIX>& _tmp); // reshape wfc into/from the nstates x CMATRIX(Nx, Ny) format



  ///=============== In the Wfcgrid2_updates.cpp ====================  
  ///< Precompute Hamiltonian on the grid
  void update_Hamiltonian(bp::object py_funct, bp::object params, int rep);


  ///< Convert diabatic and adiabatic wfc into one another
  void update_adiabatic();
  void update_diabatic();


  // Print 1D and 2D wavefunctions to file
  void print_wfc_1D(std::string prefix, int rep, vector<int>& states, int do_real, int do_imag, int do_dens);
  void print_reci_wfc_1D(std::string prefix, int rep, vector<int>& states, int do_real, int do_imag, int do_dens);

  void print_wfc_2D(std::string prefix, int rep, int state, int do_real, int do_imag, int do_dens);
  void print_reci_wfc_2D(std::string prefix, int rep, int state, int do_real, int do_imag, int do_dens);


/* 
  // Print state-resolved populations 
  double print_populations_1D(string filename,int snap);
  double print_populations_2D(string filename,int snap);

  // Flux
  void flux_1D(double xf,vector<double>& res, double m0);

  void absorb_1D(double dL,vector<double>& Pops_l,vector<double>& Pops_r);
  boost::python::list absorb_1D(double dL);

*/

}; //  class Wfcgrid2



}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

#endif  // WFCGRID2_H
