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
    
//    cs = 0.0;   si = 0.0;
    for(int nst=0;nst<nstates;nst++){  
      //cs.set(nst, nst, std::cos(dt*Hadi[npt1].get(nst, nst)) * one );
      //si.set(nst, nst, std::sin(dt*Hadi[npt1].get(nst, nst)) * one );
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
   
  //--------------------- exp(-dt*i/hbar*H_non-loc) ----------------------
  // PSI(r)->PSI(k)=reciPSI

  update_reciprocal(0);

  // Propagate in reciprocal space, for all grid points
  for(npt1=0; npt1<Npts; npt1++){ reciPSI_dia[npt1] = expK[npt1] * reciPSI_dia[npt1];  }
  
  // PSI(k)=reciPSI -> PSI(r)
  update_real(0);

  //--------------------- exp(-0.5*dt*i/hbar*H_loc) ---------------------
  for(npt1=0; npt1<Npts; npt1++){  PSI_dia[npt1] = expH[npt1] * PSI_dia[npt1];  }

}// void Wfcgrid::SOFT_propagate()





}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

