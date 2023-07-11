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
  \file Wfcgrid2_ColbertMiller.cpp
  \brief The file implements the TD-SE integrators based on the Colbert-Miller DVR
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


int points_on_same_line(int idof, vector<int>& pt1, vector<int>& pt2){
/**

 This function checks if the two multidimensional points are on
 the same line, which is the axis along the specified direction idof, 

 Returns:
 0 (false) - if the projections along any but `idof` axis are distinct

 1 (true) - if the projections along all but `idof` axes are the same

*/

  int dof;
  int res = 1;

  int ndof = pt1.size();

  for(dof=0; dof<ndof && res==1; dof++){

    if(dof!=idof){

      if(pt1[dof] != pt2[dof]){  

        res = 0;

      }
    }// exclude a specified dof axiss 

  }// for all dofs

  return res;
}


void dvr0(MATRIX& T, int npts, double rmin, double rmax, double mass){
/**
  J. Chem. Phys. 1992, 96, 1982,  Eq. A6a, A6b
*/

  int i,j;

  double L = rmax - rmin;
  double A = 0.25 * M_PI * M_PI / (mass * L * L);


  T = 0.0;

  double T_ij;

  for(i=0; i<npts; i++){

    for(j=0; j<npts; j++){

      if(i==j){  

        if(i>0){
          double si = sin(M_PI*i/npts);
          T.set(i,j, A * ( (2.0 * npts * npts + 1.0)/3.0 - 1.0/(si*si) ) );

        }
        else{  T_ij = 0.0; }

      }
      else{   
        if(i>0 && j>0){

          double si = sin(0.5*M_PI*(i-j)/npts);
          T_ij = 1.0/(si*si);

          si = sin(0.5*M_PI*(i+j)/npts);
          T_ij -= 1.0/(si*si);

          T.set(i,j, A * T_ij * pow(-1.0, i-j));

        }
        else{ T_ij = 0.0; }
      }

    }// for j
  }// for i

}



void dvr1(MATRIX& T, int npts, double dr, double mass){
/**
  J. Chem. Phys. 1992, 96, 1982,  Eq. A7
*/

  int i,j;

  double A = 0.5 / (mass * dr * dr);
  double diag = M_PI * M_PI/3.0;

  T = 0.0;

  double T_ij;

  for(i=0; i<npts; i++){
    for(j=0; j<npts; j++){

      if(i==j){  T_ij = diag;  }
      else{   T_ij = 2.0/pow((i-j),2);    }

      T_ij *= (A * pow(-1.0, i-j));

      T.set(i,j, T_ij);

    }// for j
  }// for i

}


void dvr2(MATRIX& T, int npts, double dr, double mass){
/**
  J. Chem. Phys. 1992, 96, 1982,  Eq. A8
*/

  int i,j;

  double A = 0.5 / (mass * dr * dr);

  T = 0.0;

  double T_ij;

  for(i=0; i<npts; i++){
    for(j=0; j<npts; j++){

      if(i==j){  
        if(i>0){
          T_ij = M_PI * M_PI/3.0 - 0.5/double(i*i); 
        }
      }
      else{   T_ij = 2.0/double((i-j)*(i-j)) - 2.0/double((i+j)*(i+j)) ;   }

      T_ij *= (A * pow(-1.0, i-j) );

      T.set(i,j, T_ij);

    }// for j
  }// for i

}



vector<CMATRIX> Wfcgrid2::T_PSI(vector<CMATRIX>& inp_psi, vector<int>& bc_type, vector<double>& mass, complex<double> scaling){
/**
  \brief Application of the kinetic energy operator to the wavefunction, according to Colbert-Miller formula

  Optimized version
  
  T_psi = scaling * T * inp_psi

  bc_type - the boundary condition along each dof:

  0 - finite (the most general case) - (a, b) - use rmin and rmax for that dof, in this case
  1 - (-inf, +inf)
  2 - (0, +inf)

*/

  int idof, jdof, ipt1, ipt2, istate, i, j;
  vector<int> pt1(ndof, 0);
  vector<int> pt2(ndof, 0);  

  vector<CMATRIX> Tpsi(Npts, CMATRIX(nstates, 1) );

  // Precompute 1D matrices for all dimensions 
  vector<MATRIX> reT;

  for(idof=0; idof<ndof; idof++){

    reT.push_back( MATRIX(npts[idof], npts[idof]) );

    if(bc_type[idof]==0){  dvr0(reT[idof], npts[idof], rmin[idof], rmax[idof], mass[idof]);  }
    else if(bc_type[idof]==1){  dvr1(reT[idof], npts[idof], dr[idof], mass[idof]);   }
    else if(bc_type[idof]==2){  dvr2(reT[idof], npts[idof], dr[idof], mass[idof]);   }

  }


  // Matrix multiplication  
  for(ipt1=0; ipt1<Npts; ipt1++){

    pt1 = gmap[ipt1];

    for(idof=0; idof<ndof; idof++){

      pt2 = pt1;
      i = pt1[idof];

      for(j = 0; j< npts[idof]; j++){

        pt2[idof] = j; 
        ipt2 = imap(pt2);

        for(istate=0; istate<nstates; istate++){
          Tpsi[ipt1].add(istate, 0,  reT[idof].get(i,j) * inp_psi[ipt2].get(istate, 0) );
        }// for istate

      }// for j
    }// for idof

  }// for ipt1


  // Scale
  for(ipt1=0; ipt1<Npts; ipt1++){  Tpsi[ipt1] *= scaling;  }

  return Tpsi;

}




vector<CMATRIX> Wfcgrid2::T_PSI_adi(vector<int>& bc_type, vector<double>& mass, complex<double> scaling){

  return T_PSI(PSI_adi, bc_type, mass, scaling);

}

vector<CMATRIX> Wfcgrid2::T_PSI_dia(vector<int>& bc_type, vector<double>& mass, complex<double> scaling){

  return T_PSI(PSI_dia, bc_type, mass, scaling);

}



CMATRIX Wfcgrid2::operator_T(vector<int>& bc_type, vector<double>& mass, complex<double> scaling){
/**
  \brief Compute real-space operator T according to Colbert-Miller scheme

  bc_type - the boundary condition along each dof:

  0 - finite (the most general case) - (a, b) - use rmin and rmax for that dof, in this case
  1 - (-inf, +inf)
  2 - (0, +inf)

*/

  int idof, jdof, ipt1, ipt2, istate, i, j;
  vector<int> pt1(ndof, 0);
  vector<int> pt2(ndof, 0);  

  CMATRIX T_mat(Npts, Npts);

  // Precompute 1D matrices for all dimensions 
  vector<MATRIX> reT;

  for(idof=0; idof<ndof; idof++){

    reT.push_back( MATRIX(npts[idof], npts[idof]) );

    if(bc_type[idof]==0){  dvr0(reT[idof], npts[idof], rmin[idof], rmax[idof], mass[idof]);  }
    else if(bc_type[idof]==1){  dvr1(reT[idof], npts[idof], dr[idof], mass[idof]);   }
    else if(bc_type[idof]==2){  dvr2(reT[idof], npts[idof], dr[idof], mass[idof]);   }

  }


  // Matrix multiplication  
  for(ipt1=0; ipt1<Npts; ipt1++){

    pt1 = gmap[ipt1];

    for(idof=0; idof<ndof; idof++){

      pt2 = pt1;
      i = pt1[idof];

      for(j = 0; j<npts[idof]; j++){

        pt2[idof] = j; 
        ipt2 = imap(pt2);

        T_mat.add(ipt1, ipt2,  scaling * reT[idof].get(i,j) );

      }// for j
    }// for idof

  }// for ipt1


  return T_mat;
}




vector<CMATRIX> Wfcgrid2::expT_PSI(double dt, vector<double>& mass, int nterms){
//  Compute exp(-i*dt*T) operator in the real space

  int i, j, term;

  complex<double> unit(1.0, 1.0);
  complex<double> scl(0.0, -dt);

  vector<int> bc_type(ndof, 1); // Colbert-Miller -inf to +inf case
  vector<CMATRIX> tpsi(Npts, CMATRIX(nstates, 1));  
  vector<CMATRIX> res(Npts, CMATRIX(nstates, 1));  

  // Initialize it as the probe psi
  for(i=0; i<Npts; i++){  
    for(j=0; j<nstates; j++){  tpsi[i].set(j, 0, unit);  }
  }

  double fact = 1.0;
  for(term=0; term < nterms; term++){

    for(i=0; i<Npts; i++){   res[i] += fact * tpsi[i];     }

    tpsi = T_PSI(tpsi, bc_type, mass, scl);
    
    fact = fact / double(term+1);
 
  }


  // Normalize  - the norm should be 2*Npts
  double nrm = 0.0;
  for(i=0; i<Npts; i++){  nrm += (res[i].H() * res[i]).get(0, 0).real();  }

  nrm = sqrt(2.0*Npts / nrm); 

  for(i=0; i<Npts; i++){  res[i] *= complex<double>(nrm, 0.0);  }

  return res;
}



void Wfcgrid2::Colbert_Miller_SOFT(CMATRIX& expT, vector<CMATRIX>& expV, int opt){
/**
   expT  - CMATRIX(Npts, Npts)
   expV  - Npts x CMATRIX(nstates, nstates)
   
*/
  int i,j, istate;
  CMATRIX psi_tmp(Npts, 1); 
  

  if(opt==0){
    for(i=0; i<Npts; i++){  PSI_dia[i] = expV[i] * PSI_dia[i];   }


    for(istate=0; istate<nstates; istate++){

      for(i=0; i<Npts; i++){   psi_tmp.set(i, 0, PSI_dia[i].get(istate, 0) );  }
      psi_tmp = expT * psi_tmp;
      for(i=0; i<Npts; i++){   PSI_dia[i].set(istate, 0, psi_tmp.get(i, 0));  }

    }// for istate


    for(i=0; i<Npts; i++){  PSI_dia[i] = expV[i] * PSI_dia[i];   }

  }

  else if(opt==1){

    for(istate=0; istate<nstates; istate++){

      for(i=0; i<Npts; i++){   psi_tmp.set(i, 0, PSI_dia[i].get(istate, 0) );  }
      psi_tmp = expT * psi_tmp;
      for(i=0; i<Npts; i++){   PSI_dia[i].set(istate, 0, psi_tmp.get(i, 0));  }

    }// for istate

    for(i=0; i<Npts; i++){  PSI_dia[i] = expV[i] * PSI_dia[i];   }

    for(istate=0; istate<nstates; istate){

      for(i=0; i<Npts; i++){   psi_tmp.set(i, 0, PSI_dia[i].get(istate, 0) );  }
      psi_tmp = expT * psi_tmp;
      for(i=0; i<Npts; i++){   PSI_dia[i].set(istate, 0, psi_tmp.get(i, 0));  }

    }// for istate

  }

}


vector<CMATRIX> Wfcgrid2::nubla_PSI(int idof, vector<CMATRIX>& inp_psi, complex<double> scaling){


  int ipt1, ipt2, i, j, istate;
  vector<int> pt1(ndof, 0);
  vector<int> pt2(ndof, 0);  


  vector<CMATRIX> res(Npts, CMATRIX(nstates, 1) );
  CMATRIX T(nstates, nstates);

  double T_ij = 0.0;

  for(ipt1=0; ipt1<Npts; ipt1++){

    pt1 = gmap[ipt1];

    for(ipt2=0; ipt2<Npts; ipt2++){

      pt2 = gmap[ipt2];


      if( points_on_same_line(idof, pt1, pt2)==1 ){

          i = pt1[idof];
          j = pt2[idof];

          /// ==========(-inf, +inf) ==========          

          if(i==j){  T_ij = 0.0;   }
          else{   T_ij = pow(-1.0, i-j)/(dr[idof]*(i-j));    }


          for(istate=0; istate<nstates; istate++){
            T.set(istate, istate, complex<double>(T_ij, 0.0));
          }
        
          res[ipt1] += T * inp_psi[ipt2];

        }

    }// for ipt2
  }// for ipt1


  // Scale
  for(ipt1=0; ipt1<Npts; ipt1++){
    res[ipt1] *= scaling;
  }

  return res;

}

vector<CMATRIX> Wfcgrid2::nubla_PSI_adi(int idof, complex<double> scaling){

  return nubla_PSI(idof, PSI_adi, scaling);

}

vector<CMATRIX> Wfcgrid2::nubla_PSI_dia(int idof, complex<double> scaling){

  return nubla_PSI(idof, PSI_dia, scaling);

}



vector<CMATRIX> Wfcgrid2::V_PSI(vector<CMATRIX>& V, vector<CMATRIX>& inp_psi, complex<double> scaling){
/**
  \brief Application of the potential energy operator to the wavefunction
  
  V_psi = scaling * V * inp_psi

*/

  int ipt1;
  vector<int> pt1(ndof, 0);
  vector<int> pt2(ndof, 0);  

  vector<CMATRIX> res(Npts, CMATRIX(nstates, 1) );


  for(ipt1=0; ipt1<Npts; ipt1++){
    res[ipt1] = scaling * V[ipt1] * inp_psi[ipt1];
  }

  return res;

}


vector<CMATRIX> Wfcgrid2::V_PSI_adi(complex<double> scaling){

  return V_PSI(Hadi, PSI_adi, scaling);

}

vector<CMATRIX> Wfcgrid2::V_PSI_dia(complex<double> scaling){

  return V_PSI(Hdia, PSI_dia, scaling);

}



vector<CMATRIX> Wfcgrid2::H_PSI_adi(vector<double>& mass){
/*  Compute the action of the H(effective) on the adiabatic wfc

  i * hbar * d ( xi * psi )/dt = ( T + V ) (xi * psi)

  i * hbar * psi_a * d(xi_a)/dt  = (T xi * psi - sum_{alp} {  nubla_{alp} xi / m_{alp} } * nubla_{alp} psi +  psi * T + (V - nac2/2*m) * xi * psi )
*/

  int ipt, idof;
  complex<double> one(1.0, 0.0);
  vector<int> bc_type(ndof, 1);

  vector<CMATRIX> res(Npts, CMATRIX(nstates, 1) );
  vector<CMATRIX> tmp_psi(Npts, CMATRIX(nstates, 1) );

  // Kinetic energy terms 
  tmp_psi = T_PSI_adi(bc_type, mass, one);  
  for(ipt=0; ipt<Npts; ipt++){  res[ipt] += tmp_psi[ipt];  }

  // Nonadiabatic coupling terms 
  for(idof=0; idof<ndof; idof++){

    tmp_psi = nubla_PSI_adi(idof, complex<double>(-1.0/mass[idof], 0.0)); 
    for(ipt=0; ipt<Npts; ipt++){   res[ipt] += NAC1[ipt][idof] * tmp_psi[ipt];   }
  }

  // Potential and second-order NACs energy
  for(ipt=0; ipt<Npts; ipt++){  
    CMATRIX tmp(nstates, nstates);

    for(idof=0; idof<ndof; idof++){  tmp += 0.5 * NAC2[ipt][idof] / mass[idof] ;  }    
    res[ipt] += (Hadi[ipt] - tmp ) * PSI_adi[ipt]; 
  }

  return res;

}



vector<CMATRIX> Wfcgrid2::H_PSI_dia(vector<double>& mass){
/*  Compute the action of the H(effective) on the adiabatic wfc

  i * hbar * d ( xi * psi )/dt = ( T + V ) (xi * psi)

  i * hbar * psi_a * d(xi_a)/dt  = (T xi * psi +  V * xi * psi )
*/

  int ipt, idof;
  complex<double> one(1.0, 0.0);
  vector<int> bc_type(ndof, 1);

  vector<CMATRIX> res(Npts, CMATRIX(nstates, 1) );

  // Kinetic energy terms 
  res = T_PSI_dia(bc_type, mass, one);  

  // Potential energy
  for(ipt=0; ipt<Npts; ipt++){   res[ipt] += Hdia[ipt] * PSI_dia[ipt];  }

  return res;

}




void Wfcgrid2::Colbert_Miller_propagate_adi1(double dt, vector<double>& mass){
/**
  \brief Propagator for nd-D grid wavefunction

  First order in time
  For adiabatic representation:

  i * h [ psi(t+dt,x) - psi(t,x)] / dt = [ H psi ](x)

*/

  int npt1;
  complex<double> eye(0.0, 1.0);
  vector<CMATRIX> Hpsi(Npts, CMATRIX(nstates, 1));
  vector<CMATRIX> PSI_tmp(PSI_adi);
 

  Hpsi = H_PSI_adi(mass);
  
  for(npt1=0; npt1<Npts; npt1++){  PSI_tmp[npt1] -= eye * dt * Hpsi[npt1];  }

  /// Copy the current wavefunction into the previous one
  for(npt1=0; npt1<Npts; npt1++){ PSI_adi_past[npt1] = PSI_adi[npt1];  }

  /// Copy the temporary results into the current wavefunction
  for(npt1=0; npt1<Npts; npt1++){ PSI_adi[npt1] = PSI_tmp[npt1];  }

}


void Wfcgrid2::Colbert_Miller_propagate_adi2(double dt, vector<double>& mass){
/**
  \brief Propagator for nd-D grid wavefunction

  Second order in time
  For adiabatic representation:

  i * h [ psi(t+dt,x) - psi(t-dt,x)]/(2*dt) = [ H psi ](x)

*/

  int npt1;
  complex<double> eye(0.0, 1.0);
  vector<CMATRIX> Hpsi(Npts, CMATRIX(nstates, 1));
  vector<CMATRIX> PSI_tmp(PSI_adi_past);
 
  Hpsi = H_PSI_adi(mass);

  
  for(npt1=0; npt1<Npts; npt1++){  PSI_tmp[npt1] -= 2.0 * eye * dt * Hpsi[npt1];  }

  /// Copy the current wavefunction into the previous one
  for(npt1=0; npt1<Npts; npt1++){ PSI_adi_past[npt1] = PSI_adi[npt1];  }

  /// Copy the temporary results into the current wavefunction
  for(npt1=0; npt1<Npts; npt1++){ PSI_adi[npt1] = PSI_tmp[npt1];  }

}





void Wfcgrid2::Colbert_Miller_propagate_dia1(double dt, vector<double>& mass){
/**
  \brief Propagator for nd-D grid wavefunction

  First order in time
  For diabatic representation:

  i * h [ psi(t+dt,x) - psi(t,x)] / dt = [ H psi ](x)

*/

  int npt1;
  complex<double> eye(0.0, 1.0);
  vector<CMATRIX> Hpsi(Npts, CMATRIX(nstates, 1));
  vector<CMATRIX> PSI_tmp(PSI_dia);
 

  Hpsi = H_PSI_dia(mass);
  
  for(npt1=0; npt1<Npts; npt1++){  PSI_tmp[npt1] -= eye * dt * Hpsi[npt1];  }

  /// Copy the current wavefunction into the previous one
  for(npt1=0; npt1<Npts; npt1++){ PSI_dia_past[npt1] = PSI_dia[npt1];  }

  /// Copy the temporary results into the current wavefunction
  for(npt1=0; npt1<Npts; npt1++){ PSI_dia[npt1] = PSI_tmp[npt1];  }

}


void Wfcgrid2::Colbert_Miller_propagate_dia2(double dt, vector<double>& mass){
/**
  \brief Propagator for nd-D grid wavefunction

  Second order in time
  For diabatic representation:

  i * h [ psi(t+dt,x) - psi(t-dt,x)]/(2*dt) = [ H psi ](x)

*/

  int npt1;
  complex<double> eye(0.0, 1.0);
  vector<CMATRIX> Hpsi(Npts, CMATRIX(nstates, 1));
  vector<CMATRIX> PSI_tmp(PSI_dia_past);
 
  Hpsi = H_PSI_dia(mass);

  
  for(npt1=0; npt1<Npts; npt1++){  PSI_tmp[npt1] -= 2.0 * eye * dt * Hpsi[npt1];  }

  /// Copy the current wavefunction into the previous one
  for(npt1=0; npt1<Npts; npt1++){ PSI_dia_past[npt1] = PSI_dia[npt1];  }

  /// Copy the temporary results into the current wavefunction
  for(npt1=0; npt1<Npts; npt1++){ PSI_dia[npt1] = PSI_tmp[npt1];  }

}







}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

