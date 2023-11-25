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
  \file Wfcgrid2_SOFT.cpp
  \brief The file implements the TD-SE integrators based on SOFT (Split-Operator-Fourier-Transform)
  method of Kosloff and Kosloff
    
*/

#include "Wfcgrid2.h"
#include "../../math_meigen/libmeigen.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

using namespace libwfcgrid;

/// libwfcgrid namespace
namespace libwfcgrid2{

using namespace libmeigen;




void Wfcgrid2::update_propagator_H(double dt){
/**
  \brief Update the exp(-i*dt*Hamiltonian) for nd-D grid

  Should have called ```update_Hamiltonian``` prior to this

  H * U = U * E  =>  H = U * E * U^+

  exp(-i*dt*H) = U * exp(-i*dt*E) * U^+

  exp(-i*dt*E) = exp(-i * dt * [Re(E)+i*Im(E)]) = exp(-i*dt*Re(E)) * exp(dt*Im(E))
                                                  ----- ex -------   ---- scl ----

  ex = cos(dt*Re(E)) - i * sin(dt*Re(E))
*/

  int nst, nst1;

  CMATRIX S(nstates, nstates); S.Init_Unit_Matrix(1.0);
//  CMATRIX si(nstates, nstates);  // sin(dt*E), where E is adiabatic Ham.
//  CMATRIX cs(nstates, nstates);  // cos(dt*E), where E is adiabatic Ham.
//  CMATRIX ex(nstates, nstates);  // exp()

  complex<double> one(1.0, 0.0);
  complex<double> eye(0.0, 1.0);
  //complex<double> ex(0.0, 0.0);
  double scl = 1.0;
  double si, cs;

  // For all grid points
  for(int npt1=0; npt1<Npts; npt1++){

    // Transformation to adiabatic basis
    // Hdia * U = S * U * Hadi
    solve_eigen(Hdia[npt1], S, Hadi[npt1], U[npt1], 0); 
  }

  //============= Add phase correction ================

  if(ndof==1){
    for(int npt1=0; npt1<Npts-1; npt1++){
      CMATRIX ovlp(U[npt1].H() * U[npt1+1]);
      for(int st=0; st<nstates; st++){
        if(ovlp.get(st, st).real()<0.0){
          U[npt1+1].scale(-1, st, complex<double>(-1.0, 0.0));
        }
      }// for st
    }// for npt1
  }// if ndof==1
  else{
    cout<<"WARNING: Beware of the phase problem when initializing the dynamics!\n";
    cout<<"for problems with ndof > 1, no phase consistency correction applied yet,\n";
    cout<<"which may result in wrong kinetic energy and overall dynamics problems\n";
  }


  // For all grid points 
  for(int npt1=0; npt1<Npts; npt1++){
    
    for(int nst=0;nst<nstates;nst++){  
      cs = std::cos( dt*( Hadi[npt1].get(nst, nst).real() + Vcomplex[npt1].get(nst, nst).real() ) );
      si = std::sin( dt*( Hadi[npt1].get(nst, nst).real() + Vcomplex[npt1].get(nst, nst).real() ) );
      scl = std::exp(dt*( Hadi[npt1].get(nst, nst).imag() + Vcomplex[npt1].get(nst, nst).imag() ) );
       
      complex<double> ex(cs, -si);
      expH[npt1].set(nst, nst,  scl * ex);
    }

    // Transform cs and si according to matrix U:
    //cs = U[npt1] * cs * U[npt1].H(); 
    //si = U[npt1] * si * U[npt1].H(); 
    expH[npt1] = U[npt1] * expH[npt1] * U[npt1].H();

    // Finally construct complex exp(-i*dt*H) matrix from the cs and si matrices
    //expH[npt1] = cs - eye*si;


  }// for npt1  

}// update_propagator_H



void Wfcgrid2::update_propagator_H_lin(double dt){
/**
  \brief Update the exp(-i*dt*Hamiltonian) for nd-D grid

  Should have called ```update_Hamiltonian``` prior to this

  H * U = U * E  =>  H = U * E * U^+

  exp(-i*dt*H) = U * exp(-i*dt*E) * U^+

  exp(-i*dt*E) = exp(-i * dt * [Re(E)+i*Im(E)]) = exp(-i*dt*Re(E)) * exp(dt*Im(E))
                                                  ----- ex -------   ---- scl ----

  ex = cos(dt*Re(E)) - i * sin(dt*Re(E))

  Instead of the `update_propagator_H`, this function work in the full space,
  not point-by-point
*/

  int i,j, a;
  int nst, nst1;

  int dim = nstates * Npts;

  CMATRIX S(dim, dim); S.Init_Unit_Matrix(1.0);
  CMATRIX tmp(dim, dim);
  CMATRIX U(dim, dim);
  CMATRIX Udag(dim, dim);

  complex<double> one(1.0, 0.0);
  complex<double> eye(0.0, 1.0);
  double scl = 1.0;
  double si, cs;

  // Transformation to adiabatic basis
  // lin_Hdia * lin_U = S * lin_U * lin_Hadi
  solve_eigen(*lin_Hdia, S, tmp, U, 0);


  int cnt = 0;
  for(i=0; i<dim; i++){
    for(j=0; j<dim;j++){
      if(std::abs(U.get(i,j)) >1e-5){  cnt += 1; }
    }
  }
  cout<<"number of non-zero elements = "<<cnt<<" out of "<<(dim*dim)<<endl;

 
  *lin_Hadi = tmp;
  *lin_U = U;
  Udag = U.H();

  vector< complex<double> > expE(dim, complex<double>(0.0, 0.0));

  for(i=0; i<dim; i++){

    cs = std::cos( dt*( lin_Hadi->get(i, i).real() ) );
    si = std::sin( dt*( lin_Hadi->get(i, i).real() ) );
    scl = std::exp(dt*( lin_Hadi->get(i, i).imag() ) );

    complex<double> ex(cs, -si);
    //lin_expH->set(i, i, scl * ex);
    expE[i] = scl * ex;
  }


  // Transform cs and si according to matrix U:
  //cs = lin_U * cs * lin_U.H();
  //si = lin_U * si * lin_U.H();

// Instead of this...
//  *lin_expH = (*lin_U) * (*lin_expH) * lin_U->H();

// We do this...
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      complex<double> sum; sum = 0.0;
      for(a=0; a<dim; a++){  sum += U.get(i,a) * expE[a] * Udag.get(a, j); }  
      lin_expH->set(i, j, sum);
    }
  }// for i

  // Finally construct complex exp(-i*dt*H) matrix from the cs and si matrices
  //lin_expH = cs - eye*si;


}// update_propagator_H_lin





void Wfcgrid2::update_propagator_K(double dt, vector<double>& mass){
/**
  \brief Update reciprocal-space propagators for nd-D grid
  \param[in] dt Integration time [a.u.]
  \param[in] mass Masses of the particle in all dimensions (effective DOF) [a.u.]

  working in atomic units: hbar = 1
*/

  int idof, ipt;
  double k, kfactor;

  for(int npt1=0; npt1<Npts; npt1++){

    kfactor = 0.0;
    for(idof=0; idof<ndof; idof++){
      ipt = gmap[npt1][idof];

      k = kgrid[idof]->get(ipt);
      kfactor += k*k/mass[idof];
    } 

    kfactor *= -(2.0*M_PI*M_PI*dt);

    for(int nst=0;nst<nstates;nst++){  expK[npt1].set(nst, nst, complex<double>(std::cos(kfactor),std::sin(kfactor)) );   }//for nst

  }// for npt1


}// update_propagator_K



void Wfcgrid2::SOFT_propagate(){
/**
  \brief Propagator for nd-D grid wavefunction

  All the propagation is done in the diabatic representation

*/

  int npt1;

  //=================== Wavefunction propagation part ===================
  //--------------------- exp(-0.5*dt*i/hbar*H_loc) ---------------------
  // For all grid points 
  for(npt1=0; npt1<Npts; npt1++){ PSI_dia[npt1] = expH[npt1] * PSI_dia[npt1];  }
/*
  convert_PSI(0, 1); // dia; direct -> lin
  *lin_PSI_dia = (*lin_expH) * (*lin_PSI_dia);
  convert_PSI(0, -1); // dia; lin-> direct
*/   
  //--------------------- exp(-dt*i/hbar*H_non-loc) ----------------------
  // PSI(r)->PSI(k)=reciPSI
  update_reciprocal(0);

  // Propagate in reciprocal space, for all grid points
  for(npt1=0; npt1<Npts; npt1++){ reciPSI_dia[npt1] = expK[npt1] * reciPSI_dia[npt1];  }
  
  // PSI(k)=reciPSI -> PSI(r)
  update_real(0);

  //--------------------- exp(-0.5*dt*i/hbar*H_loc) ---------------------
  for(npt1=0; npt1<Npts; npt1++){  PSI_dia[npt1] = expH[npt1] * PSI_dia[npt1];  }
/*
  convert_PSI(0, 1); // dia; direct -> lin
  *lin_PSI_dia = (*lin_expH) * (*lin_PSI_dia);
  convert_PSI(0, -1); // dia; lin-> direct
*/

}// void Wfcgrid::SOFT_propagate()





}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

