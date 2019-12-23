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

  double nrm; nrm = 0.0;

  CMATRIX res(nstates, nstates);  

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





}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

