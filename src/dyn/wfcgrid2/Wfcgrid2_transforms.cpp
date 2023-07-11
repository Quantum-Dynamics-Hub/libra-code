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
  \file Wfcgrid2_transforms.cpp
  \brief The file implements wavefunction transformation functions 
    
*/

#include "Wfcgrid2.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

using namespace libwfcgrid;

/// libwfcgrid namespace
namespace libwfcgrid2{



void Wfcgrid2::update_reciprocal(int rep){
  // PSI(r)->PSI(k)=reciPSI

  int ipt, ipt2, istate, indx;

  vector<int> point(ndof,0); 

  if(ndof==1){

    CMATRIX in(npts[0],1);
    CMATRIX out(npts[0],1);

    for(istate=0;istate<nstates;istate++){ 

      ///< PSI to internal input
      for(ipt=0; ipt<npts[0]; ipt++){
        point[0] = ipt;
        indx = imap(point);
        if(rep==0){
          in.set(ipt, PSI_dia[indx].get(istate) );
        }
        else if(rep==1){
          in.set(ipt, PSI_adi[indx].get(istate) );
        }

      }

      ///< Do 1D FFT 
      cfft1(in,out,rmin[0],kmin[0],dr[0]);

      ///< Internal output to reciPSI
      for(ipt=0; ipt<npts[0]; ipt++){
        point[0] = ipt;
        indx = imap(point);

        if(rep==0){
          reciPSI_dia[indx].set(istate, out.get(ipt) );
        }
        else if(rep==1){
          reciPSI_adi[indx].set(istate, out.get(ipt) );
        }
      }

    }// for istate

  }// 1D 

  else if(ndof==2){

    CMATRIX in(npts[0], npts[1]);
    CMATRIX out(npts[0], npts[1]);

    for(istate=0;istate<nstates;istate++){ 

      ///< PSI to internal input
      for(ipt=0; ipt<npts[0]; ipt++){
        for(ipt2=0; ipt2<npts[1]; ipt2++){

          point[0] = ipt;
          point[1] = ipt2;
          indx = imap(point);

          if(rep==0){
            in.set(ipt, ipt2, PSI_dia[indx].get(istate) );
          }
          else if(rep==1){
            in.set(ipt, ipt2, PSI_adi[indx].get(istate) );
          }

        }// ipt2
      }// ipt

      ///< Do 1D FFT 
      cfft1_2D(in,out,rmin[0],rmin[1], kmin[0], kmin[1], dr[0],dr[1]);

      ///< Internal output to reciPSI
      for(ipt=0; ipt<npts[0]; ipt++){
        for(ipt2=0; ipt2<npts[1]; ipt2++){

          point[0] = ipt;
          point[1] = ipt2;
          indx = imap(point);

          if(rep==0){
            reciPSI_dia[indx].set(istate, out.get(ipt, ipt2) );
          }
          else if(rep==1){
            reciPSI_adi[indx].set(istate, out.get(ipt, ipt2) );
          }
     
        }// ipt2
      }//ipt

    }// for istate

  }// 2D 
 
  else{
    cout<<"ERROR in the Wfcgrid2: the FFT for dimensions larger than 2 is not yet implemented\n";
    cout<<"Exiting...\n";
    exit(0);
  }

}


void Wfcgrid2::update_real(int rep){
  // reciPSI = PSI(k) -> PSI(r)

  int ipt, ipt2, istate, indx;

  vector<int> point(ndof,0); 

  if(ndof==1){

    CMATRIX in(npts[0],1);
    CMATRIX out(npts[0],1);

    for(istate=0;istate<nstates;istate++){ 

      ///< PSI to internal input
      for(ipt=0; ipt<npts[0]; ipt++){
        point[0] = ipt;
        indx = imap(point);

        if(rep==0){
          in.set(ipt, reciPSI_dia[indx].get(istate) );
        }
        else if(rep==1){
          in.set(ipt, reciPSI_adi[indx].get(istate) );
        }
      }

      ///< Do 1D FFT 
      inv_cfft1(in,out,rmin[0],kmin[0],dr[0]);

      ///< Internal output to reciPSI
      for(ipt=0; ipt<npts[0]; ipt++){
        point[0] = ipt;
        indx = imap(point);
        if(rep==0){
          PSI_dia[indx].set(istate, out.get(ipt) );
        }
        else if(rep==1){
          PSI_adi[indx].set(istate, out.get(ipt) );
        }
      }

    }// for istate

  }// 1D 

  else if(ndof==2){

    CMATRIX in(npts[0], npts[1]);
    CMATRIX out(npts[0], npts[1]);

    for(istate=0;istate<nstates;istate++){ 

      ///< PSI to internal input
      for(ipt=0; ipt<npts[0]; ipt++){
        for(ipt2=0; ipt2<npts[1]; ipt2++){

          point[0] = ipt;
          point[1] = ipt2;
          indx = imap(point);


          if(rep==0){
            in.set(ipt, ipt2, reciPSI_dia[indx].get(istate) );
          }
          else if(rep==1){
            in.set(ipt, ipt2, reciPSI_adi[indx].get(istate) );
          }

        }// ipt2
      }// ipt

      ///< Do 1D FFT 
      inv_cfft1_2D(in,out,rmin[0],rmin[1], kmin[0], kmin[1], dr[0],dr[1]);

      ///< Internal output to reciPSI
      for(ipt=0; ipt<npts[0]; ipt++){
        for(ipt2=0; ipt2<npts[1]; ipt2++){

          point[0] = ipt;
          point[1] = ipt2;
          indx = imap(point);

          if(rep==0){
            PSI_dia[indx].set(istate, out.get(ipt, ipt2) );
          }
          else if(rep==1){
            PSI_adi[indx].set(istate, out.get(ipt, ipt2) );
          }
     
        }// ipt2
      }//ipt

    }// for istate

  }// 2D 
 
  else{
    cout<<"ERROR in the Wfcgrid2: the FFT for dimensions larger than 2 is not yet implemented\n";
    cout<<"Exiting...\n";
    exit(0);
  }

}



void Wfcgrid2::normalize(int rep){
/**
  Normalize the nd-D wavefunction: |psi> -> |psi> * 1/sqrt(norm)
*/

  double nrm = 1.0/sqrt(norm(rep));
  

  for(int npt=0; npt<Npts; npt++){
    if(rep==0){        PSI_dia[npt] *= nrm;   }
    else if(rep==1){   PSI_adi[npt] *= nrm;   }
  }

}


void Wfcgrid2::reshape_wfc_1D(int _rep, int _r_or_k, int _dir, vector<CMATRIX>& _tmp){
// reshape wfc into/from the nstates x CMATRIX(Nx, 1) format

  for(int istate=0;istate<nstates;istate++){
    for(int ipt=0; ipt<npts[0]; ipt++){


      if(_dir==1){ // from internal to external      
        if(_r_or_k==0){ // r-case      
          if(_rep == 0){ // diabatic      
            _tmp[istate].set(ipt, 0, PSI_dia[ipt].get(istate,0) ); 
          }
          else if(_rep==1){ // adiabatic
            _tmp[istate].set(ipt, 0, PSI_adi[ipt].get(istate,0) ); 
          }
        }// r-case
        else if(_r_or_k==1){ // k-case
          if(_rep == 0){ // diabatic      
            _tmp[istate].set(ipt, 0, reciPSI_dia[ipt].get(istate,0) ); 
          }
          else if(_rep==1){ // adiabatic
            _tmp[istate].set(ipt, 0, reciPSI_adi[ipt].get(istate,0) ); 
          }      
        }// k-case
      
      }// internal -> external

      else if(_dir==-1){  // from external to internal
        if(_r_or_k==0){ // r-case      
          if(_rep == 0){ // diabatic      
            PSI_dia[ipt].set(istate, 0, _tmp[istate].get(ipt,0) );
          }
          else if(_rep==1){ // adiabatic
            PSI_adi[ipt].set(istate, 0, _tmp[istate].get(ipt,0) );
          }
        }// r-case
        else if(_r_or_k==1){ // k-case
          if(_rep == 0){ // diabatic      
            reciPSI_dia[ipt].set(istate, 0, _tmp[istate].get(ipt,0) );
          }
          else if(_rep==1){ // adiabatic
            reciPSI_adi[ipt].set(istate, 0, _tmp[istate].get(ipt,0) );
          }      
        }// k-case
      
      }// external -> internal

    }// for ipt - points
  }// for i - states

}

void Wfcgrid2::reshape_wfc_2D(int _rep, int _r_or_k, int _dir, vector<CMATRIX>& _tmp){
// reshape wfc into/from the nstates x CMATRIX(Nx, Ny) format

  for(int istate=0;istate<nstates;istate++){
    for(int ipt1=0; ipt1<npts[0]; ipt1++){
      for(int ipt2=0; ipt2<npts[1]; ipt2++){

        int ipt = ipt1 * npts[1] + ipt2;
  
        if(_dir==1){ // from internal to external      

            
          if(_r_or_k==0){ // r-case      
            if(_rep == 0){ // diabatic      
              _tmp[istate].set(ipt1, ipt2, PSI_dia[ipt].get(istate,0) ); 
            }
            else if(_rep==1){ // adiabatic
              _tmp[istate].set(ipt1, ipt2, PSI_adi[ipt].get(istate,0) ); 
            }
          }// r-case
          else if(_r_or_k==1){ // k-case
            if(_rep == 0){ // diabatic      
              _tmp[istate].set(ipt1, ipt2, reciPSI_dia[ipt].get(istate,0) ); 
            }
            else if(_rep==1){ // adiabatic
              _tmp[istate].set(ipt1, ipt2, reciPSI_adi[ipt].get(istate,0) ); 
            }      
          }// k-case
        
        }// internal -> external
  
        else if(_dir==-1){  // from external to internal          

          if(_r_or_k==0){ // r-case      
            if(_rep == 0){ // diabatic      
              PSI_dia[ipt].set(istate, 0, _tmp[istate].get(ipt1, ipt2) );
            }
            else if(_rep==1){ // adiabatic
              PSI_adi[ipt].set(istate, 0, _tmp[istate].get(ipt1, ipt2) );
            }
          }// r-case
          else if(_r_or_k==1){ // k-case
            if(_rep == 0){ // diabatic      
              reciPSI_dia[ipt].set(istate, 0, _tmp[istate].get(ipt1, ipt2) );
            }
            else if(_rep==1){ // adiabatic
              reciPSI_adi[ipt].set(istate, 0, _tmp[istate].get(ipt1, ipt2) );
            }      
          }// k-case
        
        }// external -> internal
  
      }// for ipt2 - points
    }// for ipt1 - points
  }// for i - states


}




}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

