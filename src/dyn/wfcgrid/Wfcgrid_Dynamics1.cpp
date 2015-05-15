/*********************************************************************************
* Copyright (C) 2012 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Wfcgrid.h"

namespace libdyn{
namespace libwfcgrid{


void Wfcgrid::propagate_exact_2D(double dt, int Nmts){

  int nst,nst1,kx,ky;

  // Auxiliary object
  CMATRIX psi(Nx,Ny);   psi  = 0.0; // is a matrix placeholder
  CMATRIX psi1(Nx,Ny);  psi1 = 0.0; // is a matrix placeholder
  CMATRIX psi2(Nx,Ny);  psi2 = 0.0; // is a matrix placeholder


//  wfc->update_potential();
//  wfc->update_propagator(0.5*dt,wfc->mol->mass[0]);


/*  
  cout<<"in propagate_exact_2D\n";
  cout<<"X = \n"<<*wfc->X<<endl;
  cout<<"Y = \n"<<*wfc->Y<<endl;
  cout<<"H = \n";
  cout<<wfc->H[0][0]<<endl;
  cout<<"expH = \n";
  cout<<wfc->expH[0]<<endl;
*/

  //===================== Wavefunction propagation part ==============================
  //--------------------- exp(-0.5*dt*i/hbar*H_loc) ---------------------
  for(nst=0;nst<nstates;nst++){   psi.dot(PSI[nst],expH[nst]);  PSI[nst] = psi;   }// for nst
    
   
  //--------------------- exp(-dt*i/hbar*H_non-loc) ----------------------
  // PSI(r)->PSI(k)=reciPSI
  ft_2D(PSI,reciPSI,1,xmin,ymin,kxmin,kymin,dx,dy);


  // Diagonal part - this is all we need for adiabatic MD (no nonadiabatic couplings)
  for(nst=0;nst<nstates;nst++){   psi.dot(reciPSI[nst],expK[nst]);  reciPSI[nst] = psi;   }// for nst

/*

    //------------------ Couplings part --------------------
    for(int nmts=0;nmts<Nmts;nmts++){

      // Form Kx * reciPSI and Ky * reciPSI
      for(nst=0;nst<nstates;nst++){

        wfc->KxreciPSI[nst] = wfc->reciPSI[nst];
        wfc->KyreciPSI[nst] = wfc->reciPSI[nst];

        for(kx=0;kx<wfc->Nx;kx++){  
          for(ky=0;ky<wfc->Ny;ky++){
            wfc->KxreciPSI[nst].M[kx*Ny+ky] *= wfc->Kx.M[kx];
            wfc->KyreciPSI[nst].M[kx*Ny+ky] *= wfc->Ky.M[ky];
          }// for ky
        }// kx

      }// for nst

      // Kxrec(k) -> DxPSI(r)
      ft_2D(wfc->KxreciPSI,wfc->DxPSI,2,wfc->minx,wfc->miny,wfc->kxmin,wfc->kymin,wfc->dx,wfc->dy);
      // Kyrec(k) -> DyPSI(r)
      ft_2D(wfc->KyreciPSI,wfc->DyPSI,2,wfc->minx,wfc->miny,wfc->kxmin,wfc->kymin,wfc->dx,wfc->dy);

    
      for(nst=0;nst<wfc->nstates;nst++){
        wfc->DtreciPSI[nst] = 0.0;

        for(nst1=0;nst1<wfc->nstates;nst1++){
        
          psi1.dot(wfc->Dx[nst][nst1],wfc->DxPSI[nst1]);
          psi2.dot(wfc->Dy[nst][nst1],wfc->DyPSI[nst1]);
          psi = psi1 + psi2;

          cfft1_2D(psi,psi1,wfc->minx,wfc->miny,wfc->kxmin,wfc->kymin,wfc->dx,wfc->dy);

          // psi1 now contains the nst1-th contribution to time derivative in reciprocal space of the reciPSI[nst] (up to constant)
          
          wfc->DtreciPSI[nst] += psi1;

        }// for nst1
      }// for nst
      
      // Finally, update the wavefunction in reciprocal space - reciPSI:

      for(nst=0;nst<wfc->nstates;nst++){
        wfc->DtreciPSI[nst] *= -2.0*M_PI*(hbar/m0)*tau; // !!!!
        wfc->reciPSI[nst] += wfc->DtreciPSI[nst];
      }
    

    }// for nmts

    // At this point we have current state of the wavefunctions in the reciprocal space     
    //---------------------------------------------------
*/

    // Diagonal part    
    for(nst=0;nst<nstates;nst++){   psi.dot(reciPSI[nst],expK[nst]);  reciPSI[nst] = psi;   }// for nst


    // PSI(k)=reciPSI -> PSI(r)
    ft_2D(reciPSI,PSI,2,xmin,ymin,kxmin,kymin,dx,dy);

    //----------------- exp(-0.5*dt*i/hbar*H_loc) ------------------------
    for(nst=0;nst<nstates;nst++){   psi.dot(PSI[nst],expH[nst]);  PSI[nst] = psi;   }// for nst



}// void Wfcgrid::propagate_exact_2D(double dt, int Nmts)




}// namespace libwfcgrid
}// namespace libdyn


