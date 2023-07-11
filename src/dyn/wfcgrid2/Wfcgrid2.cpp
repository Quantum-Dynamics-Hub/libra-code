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
  \file Wfcgrid2.cpp
  \brief The file implements some basic methods of the Wfcgrid2 class: initialization, printing, memory allocation, etc.
    
*/

#include "Wfcgrid2.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

using namespace libwfcgrid;

/// libwfcgrid namespace
namespace libwfcgrid2{


void Wfcgrid2::init_numbers(vector<double>& rmin_, vector<double>& rmax_, vector<double>& dr_, int nstates_){
/**
  \brief Initialize n-D grid dimensions
  \param[in] rmin_ The minimal (leftmost) boundaries of the grid in each dimension
  \param[in] rmax_ The maximal (rightmost) boundaries of the grid in each dimension
  \param[in] dr_ The spacing between the grid points in each dimension
  \param[in] nstates_ The number of electronic states to consider

  Computes the number of points in each dimension so that the boundaries are satisfied and the number
  of point is the lowest power of 2 needed to enclose the interval. Also sets other class variables.
*/

  int dof;

  nstates = nstates_;
  ndof = rmin_.size();

  rmin = rmin_;
  rmax = rmax_;
  dr = dr_;


  if(npts.size()>0){  npts.clear(); }

  Npts = 1;
  for(dof=0; dof<ndof; dof++){

    // Expand the grid to be power of 2 along each dimension
    npts.push_back( libdyn::libwfcgrid::find_grid_size(rmin[dof],rmax[dof],dr[dof]) );
    kmin.push_back( -0.5/dr[dof] );
    dk.push_back( 1.0/((double)npts[dof]*dr[dof]));

    Npts *= npts[dof];

    cout<<"Dimension "<<dof<<" has "<<npts[dof]<<" grid points"<<endl;  
  }
  cout<<"Grid size is calculated\n";



  if(kmin.size()>0){  kmin.clear(); }
  for(dof=0; dof<ndof; dof++){
    kmin.push_back( -0.5/dr[dof] );
    cout<<"Minimal wavevector along dimension "<<dof<<" is "<<kmin[dof]<<endl;  
  }

  cout<<"Lower wavevectors are computed \n";


}// init_numbers




// Memory allocator
void Wfcgrid2::allocate(){
/**
  \brief Allocates memory for n-D wavefucntion and associted objects

  all numbers are assumed to be defined by this time, so the init_numbers() function must have been called first
*/

  int dof;
  
  cout<<"nstates = "<<nstates<<endl;
  cout<<"ndof = "<<ndof<<endl;
  for(dof=0; dof<ndof; dof++){
    cout<<"Dimension "<<dof<<" has "<<npts[dof]<<" grid points"<<endl;  
  }



  // Allocate arrays
  PSI_dia = vector<CMATRIX>(Npts, CMATRIX(nstates, 1));          ///< wavefunction  Npts x nstates x 1
  reciPSI_dia = vector<CMATRIX>(Npts, CMATRIX(nstates,1));      ///< same as PSI but in Fourier (reciprocal) space with 1.0*kmin
  
  PSI_adi = vector<CMATRIX>(Npts, CMATRIX(nstates, 1));          ///< wavefunction  Npts x nstates x 1
  reciPSI_adi = vector<CMATRIX>(Npts, CMATRIX(nstates,1));      ///< same as PSI but in Fourier (reciprocal) space with 1.0*kmin

  nabla_PSI_dia = vector< vector<CMATRIX> >(ndof, vector<CMATRIX>(Npts, CMATRIX(nstates, 1) ) );  ///   dPSI_dia/dR_dof ( r[ipt]) - nstates x 1 matrix
  nabla_PSI_adi = vector< vector<CMATRIX> >(ndof, vector<CMATRIX>(Npts, CMATRIX(nstates, 1) ) );  ///   dPSI_adi/dR_dof ( r[ipt]) - nstates x 1 matrix
  nabla_reciPSI_dia = vector< vector<CMATRIX> >(ndof, vector<CMATRIX>(Npts, CMATRIX(nstates, 1) ) );  ///   dPSI_dia/dR_dof ( r[ipt]) - nstates x 1 matrix in k-space
  nabla_reciPSI_adi = vector< vector<CMATRIX> >(ndof, vector<CMATRIX>(Npts, CMATRIX(nstates, 1) ) );  ///   dPSI_adi/dR_dof ( r[ipt]) - nstates x 1 matrix in k-space

  Hdia = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));   ///<  diabatic Hamiltoninans for all the Npts points
  Hadi = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));   ///<  adiabatic Hamiltoninans for all the Npts points
  Vcomplex = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));   ///<  complex absorbing potential
  U = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));      ///<  |adi> = |dia> * U : diabatic-to-adiabatic transformation for all the Npts points
  expH = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));   ///<  exponent of the diabatic Hamiltoninans for all the Npts points
  expK = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));   ///<  exponent of the kinetic energy propagator for all the Npts points

  
  rgrid = vector<MATRIX*>(ndof);
  kgrid = vector<MATRIX*>(ndof);
  for(dof=0; dof<ndof; dof++){
      rgrid[dof] = new MATRIX(npts[dof], 1);
      kgrid[dof] = new MATRIX(npts[dof], 1);
  }


}// allocate


void Wfcgrid2::init_grids(){
/**
  \brief Initialize n-D grid (real and reciprocal spaces)

  memory is assumed to be allocated by this time. If not, call allocate_1D() function first
*/
  int dof;

  // Initialize grids
  for(dof=0; dof<ndof; dof++){

      // Because init_grid returns data as CMATRIX(1 , npts[dof]) object
      *rgrid[dof]  = libdyn::libwfcgrid::init_grid(rmin[dof], rmax[dof] ,dr[dof]).real().T(); // real space

      // reciprocal space
      for(int nx=0;nx<npts[dof];nx++){  
        kgrid[dof]->M[nx] = kmin[dof] + nx * dk[dof];
      }   

  }
  cout<<"Grids are initialized\n";

}// init_grid



void Wfcgrid2::compute_mapping(){

  if(gmap.size()>0){  gmap.clear(); }

  gmap = libdyn::libwfcgrid::compute_mapping(gmap, npts);

}

int Wfcgrid2::imap(vector<int>& inp){

  return libdyn::libwfcgrid::compute_imapping(inp, npts);

}


// Constructor
Wfcgrid2::Wfcgrid2(vector<double>& rmin_, vector<double>& rmax_, vector<double>& dr_, int nstates_){
/**
  \brief n-D wavefunction constructors with parameters
  \param[in] rmin_ The minimal (leftmost) boundary of the grids in all dimensions
  \param[in] rmax_ The maximal (rightmost) boundary of the grids in all dimensions
  \param[in] dr_ The spacing between the grid points in all dimensions
  \param[in] nstates_ The number of electronic states

  This constructor will: 
  1) initialize numbers;
  2) allocate memory; 
  3) initialize grids
  4) setup grid mappings (direct and inverse)
*/

  init_numbers(rmin_, rmax_, dr_, nstates_);
  allocate();
  init_grids();
  compute_mapping();

}






}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

