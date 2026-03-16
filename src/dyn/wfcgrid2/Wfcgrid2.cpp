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


void Wfcgrid2::init_numbers(const vector<double>& rmin_, const vector<double>& rmax_, const vector<double>& dr_, int nstates_){
/**
  \brief Initialize n-D grid dimensions
  \param[in] rmin_ The minimal (leftmost) boundaries of the grid in each dimension
  \param[in] rmax_ The maximal (rightmost) boundaries of the grid in each dimension
  \param[in] dr_ The spacing between the grid points in each dimension
  \param[in] nstates_ The number of electronic states to consider

  Computes the number of points in each dimension so that the boundaries are satisfied and the number
  of point is the lowest power of 2 needed to enclose the interval. Also sets other class variables.
*/
  cout<<" == in init_numbers ==\n";
  int dof;

  nstates = nstates_;
  ndof = rmin_.size();

  rmin = vector<double>(rmin_);
  rmax = vector<double>(rmax_);
  dr = vector<double>(dr_);

  if(npts.size()>0){  npts.clear(); }

  Npts = 1;
  for(dof=0; dof<ndof; dof++){

    // Expand the grid to be power of 2 along each dimension
    npts.push_back( libdyn::libwfcgrid::find_grid_size(rmin[dof],rmax[dof],dr[dof]) );
    kmin.push_back( -0.5/dr[dof] );
    dk.push_back( 1.0/(((double)npts[dof])*dr[dof]));

    Npts *= npts[dof];

    cout<<"Dimension "<<dof<<" has "<<npts[dof]<<" grid points"<<endl;  
  }
  cout<<"Grid size is calculated\n";


  if(kmin.size()>0){  kmin.clear(); }
  for(dof=0; dof<ndof; dof++){
    kmin.push_back( -0.5/dr[dof] );
    cout<<"Minimal wavevector along dimension "<<dof<<" is "<<kmin[dof]<<endl;  
  }

  cout<<" == done with init_numbers ==\n";


}// init_numbers




// Memory allocator
void Wfcgrid2::allocate(){
/**
  \brief Allocates memory for n-D wavefucntion and associted objects

  all numbers are assumed to be defined by this time, so the init_numbers() function must have been called first
*/

  int dof;
  cout<<" == in allocate ==\n";
  cout<<"nstates = "<<nstates<<endl;
  cout<<"ndof = "<<ndof<<endl;
  for(dof=0; dof<ndof; dof++){  cout<<"Dimension "<<dof<<" has "<<npts[dof]<<" grid points"<<endl;   }

  // Allocate arrays
  PSI_dia = vector<CMATRIX>(Npts, CMATRIX(nstates, 1));          ///< wavefunction  Npts x nstates x 1
  reciPSI_dia = vector<CMATRIX>(Npts, CMATRIX(nstates,1));      ///< same as PSI but in Fourier (reciprocal) space with 1.0*kmin
//  lin_PSI_dia = new CMATRIX(dim, 1);
  
  PSI_adi = vector<CMATRIX>(Npts, CMATRIX(nstates, 1));          ///< wavefunction  Npts x nstates x 1
  reciPSI_adi = vector<CMATRIX>(Npts, CMATRIX(nstates,1));      ///< same as PSI but in Fourier (reciprocal) space with 1.0*kmin
//  lin_PSI_adi = new CMATRIX(dim, 1);

  nabla_PSI_dia = vector< vector<CMATRIX> >(ndof, vector<CMATRIX>(Npts, CMATRIX(nstates, 1) ) );  ///   dPSI_dia/dR_dof ( r[ipt]) - nstates x 1 matrix
  nabla_PSI_adi = vector< vector<CMATRIX> >(ndof, vector<CMATRIX>(Npts, CMATRIX(nstates, 1) ) );  ///   dPSI_adi/dR_dof ( r[ipt]) - nstates x 1 matrix
  nabla_reciPSI_dia = vector< vector<CMATRIX> >(ndof, vector<CMATRIX>(Npts, CMATRIX(nstates, 1) ) );  ///   dPSI_dia/dR_dof ( r[ipt]) - nstates x 1 matrix in k-space
  nabla_reciPSI_adi = vector< vector<CMATRIX> >(ndof, vector<CMATRIX>(Npts, CMATRIX(nstates, 1) ) );  ///   dPSI_adi/dR_dof ( r[ipt]) - nstates x 1 matrix in k-space

  Hdia = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));   ///<  diabatic Hamiltoninans for all the Npts points
//  lin_Hdia = new CMATRIX(dim, dim);
  Hadi = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));   ///<  adiabatic Hamiltoninans for all the Npts points
//  lin_Hadi = new CMATRIX(dim, dim);
  Vcomplex = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));   ///<  complex absorbing potential
  U = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));      ///<  |adi> = |dia> * U : diabatic-to-adiabatic transformation for all the Npts points
//  lin_U = new CMATRIX(dim, dim);
  expH = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));   ///<  exponent of the diabatic Hamiltoninans for all the Npts points
//  lin_expH = new CMATRIX(dim, dim);
  expK = vector<CMATRIX>(Npts, CMATRIX(nstates, nstates));   ///<  exponent of the kinetic energy propagator for all the Npts points

  
  rgrid = vector<MATRIX*>(ndof);
  kgrid = vector<MATRIX*>(ndof);
  for(dof=0; dof<ndof; dof++){
      rgrid[dof] = new MATRIX(npts[dof], 1);
      kgrid[dof] = new MATRIX(npts[dof], 1);
  }

  cout<<" == done with allocate ==\n";

}// allocate


void Wfcgrid2::init_grids(){
/**
  \brief Initialize n-D grid (real and reciprocal spaces)

  memory is assumed to be allocated by this time. If not, call allocate_1D() function first
*/
  int dof;

  cout<<"In init_grids, ndof = "<<ndof<<endl;

  // Initialize grids
  for(dof=0; dof<ndof; dof++){

      // Because init_grid returns data as CMATRIX(1 , npts[dof]) object
      *rgrid[dof]  = libdyn::libwfcgrid::init_grid(rmin[dof], rmax[dof] ,dr[dof]).real().T(); // real space

      // reciprocal space
      for(int nx=0;nx<npts[dof];nx++){  
        kgrid[dof]->M[nx] = kmin[dof] + nx * dk[dof];
      }   

      cout<<"For idof = "<<dof<<"\n";
      cout<<"Limits of the r-grid: ["<<rgrid[dof]->M[0]<<" , "<<rgrid[dof]->M[ npts[dof]-1 ]<<" ]\n";
      cout<<"Limits of the k-grid: ["<<kgrid[dof]->M[0]<<" , "<<kgrid[dof]->M[ npts[dof]-1 ]<<" ]\n";
  }
  cout<<"Grids are initialized, done with init_grids\n";

}// init_grid



void Wfcgrid2::compute_mapping(){

  if(gmap.size()>0){  gmap.clear(); }

  gmap = libdyn::libwfcgrid::compute_mapping(gmap, npts);

}

int Wfcgrid2::imap(vector<int>& inp){

  return libdyn::libwfcgrid::compute_imapping(inp, npts);

}


// Constructor
Wfcgrid2::Wfcgrid2(const vector<double>& rmin_, const vector<double>& rmax_, const vector<double>& dr_, int nstates_){
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
  cout<<" === In constructor ===\n";
  init_numbers(rmin_, rmax_, dr_, nstates_);
  allocate();
  init_grids();
  compute_mapping();
  cout<<" === done with constructor ===\n";
}


Wfcgrid2::Wfcgrid2(const Wfcgrid2& obj){

  cout<<"In copy constructor\n";
  vector<double> rmin_(obj.rmin);
  vector<double> rmax_(obj.rmax);
  vector<double> dr_(obj.dr);

  init_numbers( rmin_, rmax_, dr_, obj.nstates);
  allocate();
  init_grids();
  compute_mapping();

  // Now copy the content:
  PSI_dia = obj.PSI_dia;
  reciPSI_dia = obj.reciPSI_dia;
  
  PSI_adi = obj.PSI_adi;
  reciPSI_adi = obj.reciPSI_adi;  

  nabla_PSI_dia = obj.nabla_PSI_dia;
  nabla_PSI_adi = obj.nabla_PSI_adi;

  nabla_reciPSI_dia = obj.nabla_reciPSI_dia;
  nabla_reciPSI_adi = obj.nabla_reciPSI_adi;

  Hdia = obj.Hdia;
  Hadi = obj.Hadi;
  Vcomplex = obj.Vcomplex;
  NAC1 = obj.NAC1;
  NAC2 = obj.NAC2;
  U = obj.U;
  expH = obj.expH;
  expK = obj.expK; 

  cout<<"done with copy constructor\n";
}




void Wfcgrid2::convert_PSI(int _rep, int _dir){
/** 
  converts  PSI_dia (PSI_adi) <-> lin_PSI_dia (lin_PSI_adi)
  _rep - representation: 0 - diabatic, 1 - adiabatic
  _dir - direction: 1 - forward: PSI -> lin_PSI;  -1 - backward: lin_PSI -> PSI
*/
  int i, ipt;
  
  if(_rep==0 && _dir==1){ 
    for(i=0; i<nstates; i++){ 
      for(ipt=0;ipt<Npts;ipt++){ lin_PSI_dia->set(i*Npts + ipt, 0,  PSI_dia[ipt].get(i, 0) );  }
    }
  }// diabtic, direct-> lin

  if(_rep==1 && _dir==1){
    for(i=0; i<nstates; i++){
      for(ipt=0;ipt<Npts;ipt++){ lin_PSI_adi->set(i*Npts + ipt, 0,  PSI_adi[ipt].get(i, 0) );  }
    }
  }// adiabtic, direct-> lin

  if(_rep==0 && _dir==-1){
    for(i=0; i<nstates; i++){
      for(ipt=0;ipt<Npts;ipt++){ PSI_dia[ipt].set(i, 0, lin_PSI_dia->get(i*Npts + ipt, 0 ) );  }
    }
  }// diabtic, lin -> direct

  if(_rep==1 && _dir==-1){
    for(i=0; i<nstates; i++){
      for(ipt=0;ipt<Npts;ipt++){ PSI_adi[ipt].set(i, 0, lin_PSI_adi->get(i*Npts + ipt, 0 ) );  }
    }
  }// adiabtic, lin -> direct



/**
  for(int i=0; i<nstates; i++){  
    for(int ipt=0;ipt<Npts;ipt++){      

      if(_rep==0){
        if(_dir==1){           lin_PSI_dia->set(i*Npts + ipt, 0,  PSI_dia[ipt].get(i, 0) );     }
        else if(_dir==-1){     PSI_dia[ipt].set(i, 0, lin_PSI_dia->get(i*Npts + ipt, 0 ) );     }
      }

      else if(_rep==1){
        if(_dir==1){           lin_PSI_adi->set(i*Npts + ipt, 0,  PSI_adi[ipt].get(i, 0) );     }
        else if(_dir==-1){     PSI_adi[ipt].set(i, 0, lin_PSI_adi->get(i*Npts + ipt, 0 ) );     }
      }

    }// for ipt - points  
  }// for i - states
*/

}

void Wfcgrid2::convert_Ham(int _rep, int _dir){
/** 
  converts  Hdia (Hadi) <-> lin_Hdia (lin_Hadi)
  _rep - representation: 0 - diabatic, 1 - adiabatic
  _dir - direction: 1 - forward: Ham -> lin_Ham;  -1 - backward: lin_Ham -> Ham
*/

  for(int i=0; i<nstates; i++){
    for(int j=0; j<nstates; j++){
      for(int ipt=0;ipt<Npts;ipt++){

          if(_rep==0){
            if(_dir==1){           lin_Hdia->set(i*Npts + ipt,  j*Npts + ipt,   Hdia[ipt].get(i, j) );  }
            else if(_dir==-1){     Hdia[ipt].set(i, j, lin_Hdia->get(i*Npts + ipt, j*Npts + ipt ) );    }
          }

          else if(_rep==1){
            if(_dir==1){           lin_Hadi->set(i*Npts + ipt,  j*Npts + ipt,   Hadi[ipt].get(i, j) );  }
            else if(_dir==-1){     Hadi[ipt].set(i, j, lin_Hadi->get(i*Npts + ipt, j*Npts + ipt ) );    }
          }

      }// for ipt - points
    }// for j - states
  }// for i - states

}

Wfcgrid2::~Wfcgrid2(){

  PSI_dia.clear();
  reciPSI_dia.clear();  
//  delete lin_PSI_dia;

  PSI_adi.clear();
  reciPSI_adi.clear();
//  delete lin_PSI_adi;

  nabla_PSI_dia.clear();
  nabla_PSI_adi.clear();
  nabla_reciPSI_dia.clear();
  nabla_reciPSI_adi.clear();

  Hdia.clear();
//  delete lin_Hdia;

  Hadi.clear();
//  delete lin_Hadi;
  Vcomplex.clear();
  U.clear();
//  delete lin_U;
  expH.clear();
//  delete lin_expH;  
  expK.clear();

  for(int dof=0; dof<ndof; dof++){
      delete rgrid[dof];
      delete kgrid[dof]; 
  }
  rgrid.clear();
  kgrid.clear();



}



}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

