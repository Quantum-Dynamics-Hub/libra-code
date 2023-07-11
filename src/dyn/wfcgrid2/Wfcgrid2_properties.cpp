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
  \file Wfcgrid2_properties.cpp
  \brief The file implements functions for computing various properties of the wavefunctions
    
*/

#include "Wfcgrid2.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

using namespace libwfcgrid;

/// libwfcgrid namespace
namespace libwfcgrid2{




double Wfcgrid2::norm(int rep){
/**
  Compute the norm for nd-D wavefunction: <psi|psi>   
*/

  double nrm; nrm = 0.0;

  for(int npt1=0; npt1<Npts; npt1++){

    if(rep==0){
      nrm += (PSI_dia[npt1].H() * PSI_dia[npt1]).get(0,0).real();
    }
    else if(rep==1){
      nrm += (PSI_adi[npt1].H() * PSI_adi[npt1]).get(0,0).real();
    }

  }// for npt1


  for(int idof; idof<ndof; idof++){
    nrm *= dr[idof];
  }

  return nrm;

}




double Wfcgrid2::e_kin(vector<double>& mass, int rep){
/**
  Compute kinetic energy for nd-D wavefunction: <psi|T|psi> / <psi|psi>
  
*/
  int idof, ipt;
  double k, kfactor, kfactor2;
  double res; res = 0.0;
  double nrm; nrm = 0.0;

  for(int npt1=0; npt1<Npts; npt1++){

    kfactor = 0.0;

    for(idof=0; idof<ndof; idof++){
      ipt = gmap[npt1][idof];
      k = kgrid[idof]->get(ipt);
      kfactor += k*k/mass[idof];
    } 

    if(rep==0){

      res += kfactor * (reciPSI_dia[npt1].H() * reciPSI_dia[npt1]).get(0,0).real();
      nrm += (PSI_dia[npt1].H() * PSI_dia[npt1]).get(0,0).real();
    }
    else if(rep==1){
      res += kfactor * (reciPSI_adi[npt1].H() * reciPSI_adi[npt1]).get(0,0).real();
      nrm += (PSI_adi[npt1].H() * PSI_adi[npt1]).get(0,0).real();

    }

  }// for npt1


/*
      kfactor2 = 0.0;  

      for(idof=0; idof<ndof; idof++){
        ipt = gmap[npt1][idof];
        k = kgrid[idof]->get(ipt);
        kfactor2 += k*NAC1[npt1][idof]/mass[idof];
      } 
*/



  for(idof=0; idof<ndof; idof++){
    res *= dk[idof];
    nrm *= dr[idof];
  }
  res *= (2.0*M_PI*M_PI);


  res = res / nrm;

  return res;

}// e_kin




double Wfcgrid2::e_pot(int rep){
/**
  Compute potential energy for nd-D wavefunction: <psi|V|psi> / <psi|psi>  
*/

  double res; res = 0.0;
  double nrm; nrm = 0.0;

  for(int npt1=0; npt1<Npts; npt1++){

    if(rep==0){
      res += (PSI_dia[npt1].H() * Hdia[npt1] * PSI_dia[npt1]).get(0,0).real();
      nrm += (PSI_dia[npt1].H() * PSI_dia[npt1]).get(0,0).real();
    }
    else if(rep==1){
      res += (PSI_adi[npt1].H() * Hadi[npt1] * PSI_adi[npt1]).get(0,0).real();
      nrm += (PSI_adi[npt1].H() * PSI_adi[npt1]).get(0,0).real();
    }

  }// for npt1

  res = res / nrm;

  return res;

}// e_pot




double Wfcgrid2::e_tot(vector<double>& mass, int rep){
/**
  Compute total energy for nd-D wavefunction: <psi|T+V|psi> / <psi|psi>
*/
  double res = e_kin(mass, rep) + e_pot(rep);
  
  return res;

}// e_tot



CMATRIX Wfcgrid2::get_pow_q(int rep, int n){
/**
  Compute the expectation value of coordinates of a nd-D wavefunction: <psi|r^n |psi> / <psi|psi>

  Out:  CMATRIX(ndof, 1)
  
*/

  int idof, ipt;
  double q;
  double nrm; nrm = 0.0;

  CMATRIX res(ndof, 1);
  

  for(int npt1=0; npt1<Npts; npt1++){

    complex<double> ampl(0.0, 0.0);

    if(rep==0){

      ampl = (PSI_dia[npt1].H() * PSI_dia[npt1]).get(0,0);
      nrm += (PSI_dia[npt1].H() * PSI_dia[npt1]).get(0,0).real();
    }
    else if(rep==1){
      ampl = (PSI_adi[npt1].H() * PSI_adi[npt1]).get(0,0);
      nrm += (PSI_adi[npt1].H() * PSI_adi[npt1]).get(0,0).real();
    }

    for(idof=0; idof<ndof; idof++){

      ipt = gmap[npt1][idof];
      q = rgrid[idof]->get(ipt);

      res.add(idof, 0, pow(q, n)*ampl );
    } 

  }// for npt1

  res = res / nrm;

  return res;

}




CMATRIX Wfcgrid2::get_pow_p(int rep, int n){
/**
  Compute the expectation value of momentum of a nd-D wavefunction: <psi|(-i*hbar*d/dr)^n |psi> / <psi|psi>

  Out:  CMATRIX(ndof, 1)
  
*/

  int idof, ipt;
  double k;
  double nrm; nrm = 0.0;

  CMATRIX res(ndof, 1);
  

  for(int npt1=0; npt1<Npts; npt1++){

    complex<double> ampl(0.0, 0.0);

    if(rep==0){

      ampl = (reciPSI_dia[npt1].H() * reciPSI_dia[npt1]).get(0,0);
      nrm += (PSI_dia[npt1].H() * PSI_dia[npt1]).get(0,0).real();
    }
    else if(rep==1){
      ampl = (reciPSI_adi[npt1].H() * reciPSI_adi[npt1]).get(0,0);
      nrm += (PSI_adi[npt1].H() * PSI_adi[npt1]).get(0,0).real();
    }

    for(idof=0; idof<ndof; idof++){
      ipt = gmap[npt1][idof];
      k = kgrid[idof]->get(ipt);

      res.add(idof, 0, pow(k, n)*ampl );
    } 

  }// for npt1

  for(idof=0; idof<ndof; idof++){
    res *= dk[idof];
    nrm *= dr[idof];
  }
  res *= pow((2.0*M_PI), n);

  res = res / nrm;

  return res;

}



CMATRIX Wfcgrid2::get_den_mat(int rep){
/**
  Compute the wavefunction density matrix in a given representation

  Out:  CMATRIX(nstates, nstates)
  
*/

  CMATRIX res(nstates, nstates);  

  double nrm; nrm = 0.0;

  for(int npt1=0; npt1<Npts; npt1++){

    if(rep==0){
      res += (PSI_dia[npt1] * PSI_dia[npt1].H());
      nrm += (PSI_dia[npt1].H() * PSI_dia[npt1]).get(0,0).real();
    }
    else if(rep==1){
      res += (PSI_adi[npt1] * PSI_adi[npt1].H());
      nrm += (PSI_adi[npt1].H() * PSI_adi[npt1]).get(0,0).real();
    }

  }// for npt1

  res = res / nrm;

  return res;

}



MATRIX Wfcgrid2::get_pops(int rep){
/**
  Compute the populations of all states in a given representation

  Out:  CMATRIX(nstates, 1)
  
*/

  int istate;
  double pop_i;

  MATRIX res(nstates, 1);  

  for(int npt1=0; npt1<Npts; npt1++){

    if(rep==0){

      for(istate=0; istate<nstates; istate++){
          pop_i = (PSI_dia[npt1].get(istate, 0) * std::conj(PSI_dia[npt1].get(istate, 0))).real(); 
          res.add(istate, 0,  pop_i);
      }
    }
    else if(rep==1){

      for(istate=0; istate<nstates; istate++){
          pop_i = (PSI_adi[npt1].get(istate, 0) * std::conj(PSI_adi[npt1].get(istate, 0))).real(); 
          res.add(istate, 0,  pop_i);
      }
    }

  }// for npt1

  double dV = 1.0;
  for(int idof=0; idof<ndof; idof++){  dV *= dr[idof];  }

  res *= dV;

  return res;

}




MATRIX Wfcgrid2::get_pops(int rep, vector<double>& bmin, vector<double>& bmax){
/**
  Compute the populations of all states in a given representation

  Only the points that belong to a multidimensional box given by bmin and bmax parameters are
  included in the corresponding populations

  Out:  CMATRIX(nstates, 1)
  
*/

  int istate, idof;
  double pop_i;
  vector<int> pt;

  MATRIX res(nstates, 1);  
   
  int maxdim = ndof;
  if(bmin.size() < maxdim) { maxdim = bmin.size(); }

  for(int npt1=0; npt1<Npts; npt1++){

    pt = gmap[npt1];


    int is_included = 1;

    for(idof=0; idof<maxdim && is_included; idof++){

      double coord = rgrid[idof]->get(pt[idof]);

      if(coord < bmin[idof]  || coord > bmax[idof] ){  is_included = 0; }    

    }


    if(is_included){

      if(rep==0){    
     
        for(istate=0; istate<nstates; istate++){
            pop_i = (PSI_dia[npt1].get(istate, 0) * std::conj(PSI_dia[npt1].get(istate, 0))).real(); 
            res.add(istate, 0,  pop_i);
        }
      }
      else if(rep==1){
     
        for(istate=0; istate<nstates; istate++){
            pop_i = (PSI_adi[npt1].get(istate, 0) * std::conj(PSI_adi[npt1].get(istate, 0))).real(); 
            res.add(istate, 0,  pop_i);
        }
      }

    }// if is_included

  }// for npt1

  double dV = 1.0;
  for(int idof=0; idof<ndof; idof++){  dV *= dr[idof];  }

  res *= dV;


  return res;

}


void Wfcgrid2::compute_wfc_gradients(int rep, int idof, double mass){
// Compute wfc derivatives: first compute them in the k-space, then 
// FT to the real space
// Form Kx * reciPSI and Ky * reciPSI, etc.

  for(int ipt=0; ipt<Npts; ipt++){

//    nabla_reciPSI_dia[idof][ipt]  =  reciPSI_dia[ipt];

  }// for ipt

/*
  for(int nst=0;nst<nstates;nst++){
    KxreciPSI[nst] = reciPSI[nst];

    for(int kx=0;kx<Nx;kx++){  
      KxreciPSI[nst].M[kx] *= Kx->M[kx];
    }// kx
  }// for nst

  // Kxrec(k) -> DxPSI(r)
  ft_1D(KxreciPSI,DxPSI,2,xmin,kxmin,dx);

*/

}



//MATRIX void Wfcgrid2::flux(double xf, int fdim){
/**
  \brief Compute the population flux in 1D case
  \param[in] xf The point at which flux is computed
  \param[in] fdim index of the dimension for which the xf is set 
  \param[out] res The collector of the fluxes (for each electronic state projection of the wfc)
  \param[in] m0 The effective mass of the quantum particle (DOF)


  MATRIX res(nstates, 1);
  if(xf < rmin[fdim]) { 
    cout<<"The flux point is outside the grid boundaries xf = "<<xf<<" rmin = "<<rmin[fdim]<<"  for fdim = "<< fdim <<endl; 
    cout<<"Exiting\n";
    exit(0);
  }
  if(xf > rmax[fdim]) { 
    cout<<"The flux point is outside the grid boundaries xf = "<<xf<<" rmax = "<<rmax[fdim]<<"  for fdim = "<< fdim <<endl; 
    cout<<"Exiting\n";
    exit(0);
  }


  // index of the grid point at which the counting plane is installed  
  int ipt_const = (xf - rmin[fdim])/dr[fdim];

  vector<int> hplane;
  hplane = libwfcgrid::compute_hyperplane(npts, int fdim,ipt_const);


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

    // Original expression
//    res[nst] = 2.0*M_PI*(0.5*hbar/m0)*imag(std::conj(PSI[nst].M[i])*one*DxPSI[nst].M[i] - 
//                                           std::conj(one*DxPSI[nst].M[i])*PSI[nst].M[i] );
    // Simplified 1
//    res[nst] = 4.0*M_PI*(0.5*hbar/m0)*imag(std::conj(PSI[nst].M[i])*one*DxPSI[nst].M[i] );

    // Simplified 2
    res[nst] = (2.0*M_PI/m0) * real(std::conj(PSI[nst].M[i])*DxPSI[nst].M[i] );

  }

}// flux_1D

*/





}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

