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
  \file Wfcgrid2_direct.cpp
  \brief The file implements the TD-SE integrators based on the simplest finite difference 
  method     
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


void Wfcgrid2::direct_allocate_tmp_vars(double rep){
  /**
  \param[in] rep - representation for which we want to allocate memory: 0 - diabatic, 1 - adiabatic
  */

  // Allocate arrays
  if(rep==0){
    PSI_dia_past = vector<CMATRIX>(Npts, CMATRIX(nstates, 1));        
    PSI_dia_past = vector<CMATRIX>(Npts, CMATRIX(nstates, 1));        
  }

  if(rep==1){
    PSI_adi_past = vector<CMATRIX>(Npts, CMATRIX(nstates, 1));        
    PSI_adi_past = vector<CMATRIX>(Npts, CMATRIX(nstates, 1));        
  }

  // Allocate first- and second-order NACs

  NAC1 = vector< vector<CMATRIX> >(Npts, vector<CMATRIX>(ndof, CMATRIX(nstates, nstates) ) );        
  NAC2 = vector< vector<CMATRIX> >(Npts, vector<CMATRIX>(ndof, CMATRIX(nstates, nstates) ) );        

}


void Wfcgrid2::direct_propagate_adi2(double dt, vector<double>& mass){
/**
  \brief Propagator for nd-D grid wavefunction

  Second order in time
  For adiabatic representation

*/

  int idof, ipt, npt1, npt2;
  vector<int> pt1(ndof, 0);
  vector<int> pt2(ndof, 0);
  complex<double> one(1.0, 0.0);
  complex<double> eye(0.0, 1.0);

  vector< complex<double> > A(ndof, complex<double>(0.0, 0.0));
  for(idof=0; idof<ndof; idof++){  
    A[idof] = 0.25 * eye * dt / (mass[idof] * dr[idof]*dr[idof]);
  }


  CMATRIX xi_n_p2(nstates,1);  /// chi_{n+2}
  CMATRIX xi_n_p1(nstates,1);  /// chi_{n+1}
  CMATRIX xi_n(nstates,1);     /// chi_{n}
  CMATRIX xi_n_m1(nstates,1);  /// chi_{n-1}
  CMATRIX xi_n_m2(nstates,1);  /// chi_{n-2}

  /** =============== First terms ================ */
  vector<CMATRIX> PSI_tmp(PSI_adi_past);


  for(npt1=0; npt1<Npts; npt1++){


    pt1 = gmap[npt1];


    for(idof=0; idof<ndof; idof++){
      ipt = pt1[idof];

      /** =============== Second terms ================ */
      if(ipt+2 > npts[idof]-1){  xi_n_p2 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] += 2;
        npt2 = imap(pt2);

        xi_n_p2 = PSI_adi[npt2]; 
      }

      xi_n = PSI_adi[npt1];

      if(ipt-2 < 0){  xi_n_m2 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] -= 2;
        npt2 = imap(pt2);

        xi_n_m2 = PSI_adi[npt2]; 
      }

      
      PSI_tmp[npt1] += A[idof] * (xi_n_p2 - 2.0 * xi_n + xi_n_m2);


      /** =============== Third terms ================ */
      if(ipt+1 > npts[idof]-1){  xi_n_p1 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] += 1;
        npt2 = imap(pt2);

        xi_n_p1 = PSI_adi[npt2]; 
      }

      if(ipt-1 < 0){  xi_n_m1 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] -= 1;
        npt2 = imap(pt2);

        xi_n_m1 = PSI_adi[npt2]; 
      }

      PSI_tmp[npt1] += 4.0 * dr[idof] * A[idof] * NAC1[npt1][idof] * ( xi_n_p1 - xi_n_m1);


      /** =============== Fifth terms ================ */

      PSI_tmp[npt1] += (dt/mass[idof]) * eye * NAC2[npt1][idof] * xi_n;


    } // for idof

    /** =============== Forth terms ================ */

    PSI_tmp[npt1] -= 2.0 * eye * dt * Hadi[npt1] * PSI_adi[npt1];

  }// for npt1


  /// Copy the current wavefunction into the previous one
  for(npt1=0; npt1<Npts; npt1++){ PSI_adi_past[npt1] = PSI_adi[npt1];  }

  /// Copy the temporary results into the current wavefunction
  for(npt1=0; npt1<Npts; npt1++){ PSI_adi[npt1] = PSI_tmp[npt1];  }

}

void Wfcgrid2::direct_propagate_adi1(double dt, vector<double>& mass){
/**
  \brief Propagator for nd-D grid wavefunction

  First order in time
  For adiabatic representation

*/

  int idof, ipt, npt1, npt2;
  vector<int> pt1(ndof, 0);
  vector<int> pt2(ndof, 0);
  complex<double> one(1.0, 0.0);
  complex<double> eye(0.0, 1.0);

  vector< complex<double> > A(ndof, complex<double>(0.0, 0.0));
  for(idof=0; idof<ndof; idof++){  
    A[idof] = 0.25 * eye * dt / (mass[idof] * dr[idof]*dr[idof]);
  }


  CMATRIX xi_n_p2(nstates,1);  /// chi_{n+2}
  CMATRIX xi_n_p1(nstates,1);  /// chi_{n+1}
  CMATRIX xi_n(nstates,1);     /// chi_{n}
  CMATRIX xi_n_m1(nstates,1);  /// chi_{n-1}
  CMATRIX xi_n_m2(nstates,1);  /// chi_{n-2}

  /** =============== First terms ================ */
  vector<CMATRIX> PSI_tmp(PSI_adi);


  for(npt1=0; npt1<Npts; npt1++){


    pt1 = gmap[npt1];


    for(idof=0; idof<ndof; idof++){
      ipt = pt1[idof];

      /** =============== Second terms ================ */
      if(ipt+2 > npts[idof]-1){  xi_n_p2 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] += 2;
        npt2 = imap(pt2);

        xi_n_p2 = PSI_adi[npt2]; 
      }

      xi_n = PSI_adi[npt1];

      if(ipt-2 < 0){  xi_n_m2 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] -= 2;
        npt2 = imap(pt2);

        xi_n_m2 = PSI_adi[npt2]; 
      }

      
      PSI_tmp[npt1] += 0.5 * A[idof] * (xi_n_p2 - 2.0 * xi_n + xi_n_m2);


      /** =============== Third terms ================ */
      if(ipt+1 > npts[idof]-1 ){  xi_n_p1 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] += 1;
        npt2 = imap(pt2);

        xi_n_p1 = PSI_adi[npt2]; 
      }

      if(ipt-1 < 0){  xi_n_m1 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] -= 1;
        npt2 = imap(pt2);

        xi_n_m1 = PSI_adi[npt2]; 
      }

      PSI_tmp[npt1] += 2.0 * dr[idof] * A[idof] * NAC1[npt1][idof] * ( xi_n_p1 - xi_n_m1);


      /** =============== Fifth terms ================ */

      PSI_tmp[npt1] += (0.5*dt/mass[idof]) * eye * NAC2[npt1][idof] * xi_n;


    } // for idof

    /** =============== Forth terms ================ */

    PSI_tmp[npt1] -= eye * dt * Hadi[npt1] * PSI_adi[npt1];

  }// for npt1


  /// Copy the current wavefunction into the previous one
  for(npt1=0; npt1<Npts; npt1++){ PSI_adi_past[npt1] = PSI_adi[npt1];  }

  /// Copy the temporary results into the current wavefunction
  for(npt1=0; npt1<Npts; npt1++){ PSI_adi[npt1] = PSI_tmp[npt1];  }

}







void Wfcgrid2::direct_propagate_dia2(double dt, vector<double>& mass){
/**
  \brief Propagator for nd-D grid wavefunction

  Second order in time
  For diabatic representation

*/

  int idof, ipt, npt1, npt2;
  vector<int> pt1(ndof, 0);
  vector<int> pt2(ndof, 0);
  complex<double> one(1.0, 0.0);
  complex<double> eye(0.0, 1.0);

  vector< complex<double> > A(ndof, complex<double>(0.0, 0.0));
  for(idof=0; idof<ndof; idof++){  
    A[idof] = 0.25 * eye * dt / (mass[idof] * dr[idof]*dr[idof]);
  }


  CMATRIX xi_n_p2(nstates,1);  /// chi_{n+2}
  CMATRIX xi_n_p1(nstates,1);  /// chi_{n+1}
  CMATRIX xi_n(nstates,1);     /// chi_{n}
  CMATRIX xi_n_m1(nstates,1);  /// chi_{n-1}
  CMATRIX xi_n_m2(nstates,1);  /// chi_{n-2}

  /** =============== First terms ================ */
  vector<CMATRIX> PSI_tmp(PSI_dia_past);


  for(npt1=0; npt1<Npts; npt1++){


    pt1 = gmap[npt1];


    for(idof=0; idof<ndof; idof++){
      ipt = pt1[idof];

      /** =============== Second terms ================ */
      if(ipt+2 > npts[idof]-1){  xi_n_p2 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] += 2;
        npt2 = imap(pt2);

        xi_n_p2 = PSI_dia[npt2]; 
      }

      xi_n = PSI_dia[npt1];

      if(ipt-2 < 0){  xi_n_m2 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] -= 2;
        npt2 = imap(pt2);

        xi_n_m2 = PSI_dia[npt2]; 
      }

      
      PSI_tmp[npt1] += A[idof] * (xi_n_p2 - 2.0 * xi_n + xi_n_m2);


    } // for idof

    /** =============== Forth terms ================ */

    PSI_tmp[npt1] -= 2.0 * eye * dt * Hdia[npt1] * PSI_dia[npt1];

  }// for npt1


  /// Copy the current wavefunction into the previous one
  for(npt1=0; npt1<Npts; npt1++){ PSI_dia_past[npt1] = PSI_dia[npt1];  }

  /// Copy the temporary results into the current wavefunction
  for(npt1=0; npt1<Npts; npt1++){ PSI_dia[npt1] = PSI_tmp[npt1];  }

}

void Wfcgrid2::direct_propagate_dia1(double dt, vector<double>& mass){
/**
  \brief Propagator for nd-D grid wavefunction

  First order in time
  For diabatic representation

*/

  int idof, ipt, npt1, npt2;
  vector<int> pt1(ndof, 0);
  vector<int> pt2(ndof, 0);
  complex<double> one(1.0, 0.0);
  complex<double> eye(0.0, 1.0);

  vector< complex<double> > A(ndof, complex<double>(0.0, 0.0));
  for(idof=0; idof<ndof; idof++){  
    A[idof] = 0.25 * eye * dt / (mass[idof] * dr[idof]*dr[idof]);
  }


  CMATRIX xi_n_p2(nstates,1);  /// chi_{n+2}
  CMATRIX xi_n_p1(nstates,1);  /// chi_{n+1}
  CMATRIX xi_n(nstates,1);     /// chi_{n}
  CMATRIX xi_n_m1(nstates,1);  /// chi_{n-1}
  CMATRIX xi_n_m2(nstates,1);  /// chi_{n-2}

  /** =============== First terms ================ */
  vector<CMATRIX> PSI_tmp(PSI_dia);


  for(npt1=0; npt1<Npts; npt1++){


    pt1 = gmap[npt1];


    for(idof=0; idof<ndof; idof++){
      ipt = pt1[idof];

      /** =============== Second terms ================ */
      if(ipt+2 > npts[idof]-1){  xi_n_p2 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] += 2;
        npt2 = imap(pt2);

        xi_n_p2 = PSI_dia[npt2]; 
      }

      xi_n = PSI_dia[npt1];

      if(ipt-2 < 0){  xi_n_m2 *= 0.0; }
      else{ 
        pt2 = pt1;
        pt2[idof] -= 2;
        npt2 = imap(pt2);

        xi_n_m2 = PSI_dia[npt2]; 
      }

      
      PSI_tmp[npt1] += 0.5 * A[idof] * (xi_n_p2 - 2.0 * xi_n + xi_n_m2);


    } // for idof

    /** =============== Forth terms ================ */

    PSI_tmp[npt1] -= eye * dt * Hdia[npt1] * PSI_dia[npt1];

  }// for npt1


  /// Copy the current wavefunction into the previous one
  for(npt1=0; npt1<Npts; npt1++){ PSI_dia_past[npt1] = PSI_dia[npt1];  }

  /// Copy the temporary results into the current wavefunction
  for(npt1=0; npt1<Npts; npt1++){ PSI_dia[npt1] = PSI_tmp[npt1];  }

}




}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

