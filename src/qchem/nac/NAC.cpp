/*********************************************************************************
* Copyright (C) 2014 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
 \file NAC.cpp
 \brief Implementation of the non-adiabatic couplings 

*/
/****************************************************************************
  This file contains following functions:

  void compute_overlap_nac(Electronic* ela1, Electronic* ela2, int indx_min1, int indx_max1, vector<int>& basis_fo1,
                           Electronic* elb1, Electronic* elb2, int indx_min2, int indx_max2, vector<int>& basis_fo2,
                           vector<AO>& basis_ao1, vector<AO>& basis_ao2,
                           double dt, vector<double*>& aux, int naux, MATRIX* ovlp, MATRIX* nac, MATRIX* ene, int nac_opt)

  void collect_matrices(int Nfrag, vector<MATRIX*>& x, vector<int>& min_orbs, vector<int>& max_orbs,
                        MATRIX* X, vector<int>& frags,vector<int>& orbs)


****************************************************************************/

#include "Electronic.h"


void compute_overlap_nac(Electronic* ela1, Electronic* ela2, int indx_min1, int indx_max1, vector<int>& basis_fo1,
                         Electronic* elb1, Electronic* elb2, int indx_min2, int indx_max2, vector<int>& basis_fo2,
                         vector<AO>& basis_ao1, vector<AO>& basis_ao2,
                         double dt, vector<double*>& aux, int naux, MATRIX* ovlp, MATRIX* nac, MATRIX* ene, int nac_opt){
/*!
// nac_opt = 0 - use Tully-Hamess-Schiffer formula - Hermitian NAC
// nac_opt = 1 - add non-othogonality correction - the resulting NAC will be non-Hermitian

// ela1  = psi_A(t)
// ela2  = psi_A(t+dt)
// elb1  = psi_B(t)
// elb2  = psi_B(t+dt)
// basis_fo1 = A  - for both t and t+dt
// basis_fo2 = B  - for both t and t+dt
// bais_ao1  = t    - global to both A and B
// bais_ao2  = t+dt - global to both A and B

// NAC between i-th orbitals on fragment A and j-th orbital on fragment B
// S_AB_ij =  hbar*(<A_i(t+dt)|B_j(t+dt)> + <A_i(t)|B_j(t)> + <A_i(t+dt)|B_j(t)> + <A_i(t)|B_j(t+dt)>)/(4*dt)
*/

  double hbar = 0.658218; ///> units eV * fs


  if(ela1->Norb!=ela2->Norb){ cout<<"Error in compute_nac: dimensions of the two input electronic structures are different\n"; exit(0); }
  if(elb1->Norb!=elb2->Norb){ cout<<"Error in compute_nac: dimensions of the two input electronic structures are different\n"; exit(0); }
  if(indx_min1<0){ cout<<"Error in compute_nac: indx_min1 must be not less than 0 current value is "<<indx_min1<<"\n"; exit(0); }
  if(indx_min2<0){ cout<<"Error in compute_nac: indx_min2 must be not less than 0 current value is "<<indx_min2<<"\n"; exit(0); }
  if(indx_max1>=ela1->Norb){ cout<<"Error in compute_nac: indx_max1 must be not less than "<<ela1->Norb<<" current value is "<<indx_max1<<"\n"; exit(0); }
  if(indx_max2>=elb2->Norb){ cout<<"Error in compute_nac: indx_max2 must be not less than "<<elb2->Norb<<" current value is "<<indx_max2<<"\n"; exit(0); }


  int i,j,a,b,ii1,jj1, ii2,jj2;

  int DF = 0; ///> Debug flag - available only to coder
 

  if(DF){
    cout<<"ela1->Norb= "<<ela1->Norb<<endl;
    cout<<"elb2->Norb= "<<elb2->Norb<<endl;
  }
 
  // Precompute 
  // All matrices Sab00, Sabtt, Sab0t, and Sabt0 are of dimensions: A->Norb x B->Norb
  // <xi(t)|xi(t)>
  MATRIX* Sab00; Sab00 = new MATRIX(ela1->Norb,elb1->Norb);
  VECTOR dIdA,dIdB;
  for(a=0;a<ela1->Norb;a++){
    for(b=0;b<elb1->Norb;b++){

      Sab00->M[a*elb1->Norb+b] = OVERLAP_INTEGRAL(basis_ao1[basis_fo1[a]],basis_ao1[basis_fo2[b]],0,dIdA,dIdB,aux,naux);
          
    }// for b
  }// for a
  if(DF){ cout<<"Inside of compute_overlap <xi(t)|xi(t)> = \n Sab00=\n"<<*Sab00<<endl; }


  // <xi(t+dt)|xi(t+dt)>
  MATRIX* Sabtt; Sabtt = new MATRIX(ela2->Norb,elb2->Norb);
  for(a=0;a<ela2->Norb;a++){
    for(b=0;b<elb2->Norb;b++){

      Sabtt->M[a*elb2->Norb+b] = OVERLAP_INTEGRAL(basis_ao2[basis_fo1[a]],basis_ao2[basis_fo2[b]],0,dIdA,dIdB,aux,naux);
          
    }// for b
  }// for a
  if(DF){ cout<<"Inside of compute_overlap <xi(t+dt)|xi(t+dt)> = \n Sabtt=\n"<<*Sabtt<<endl; }


  // <xi_a(t)|xi_b(t+dt)>
  MATRIX* Sab0t; Sab0t = new MATRIX(ela1->Norb,elb2->Norb);
  for(a=0;a<ela1->Norb;a++){
    for(b=0;b<elb2->Norb;b++){

      Sab0t->M[a*elb2->Norb+b] = OVERLAP_INTEGRAL(basis_ao1[basis_fo1[a]],basis_ao2[basis_fo2[b]],0,dIdA,dIdB,aux,naux);
          
    }// for b
  }// for a
  if(DF){ cout<<"Inside of compute_overlap <xi(t)|xi(t+dt)> = \n Sab0t=\n"<<*Sab0t<<endl; }


  MATRIX* Sabt0; Sabt0 = new MATRIX(ela2->Norb,elb1->Norb);
  for(a=0;a<ela1->Norb;a++){
    for(b=0;b<elb2->Norb;b++){

      Sabt0->M[a*elb1->Norb+b] = OVERLAP_INTEGRAL(basis_ao2[basis_fo1[a]],basis_ao1[basis_fo2[b]],0,dIdA,dIdB,aux,naux);
          
    }// for b
  }// for a
  if(DF){ cout<<"Inside of compute_overlap <xi(t+dt)|xi(t)> = \n Sabt0=\n"<<*Sabt0<<endl; }



  *ovlp = 0.0;
  *nac = 0.0;
  *ene = 0.0;

  for(i=indx_min1;i<=indx_max1;i++){
    for(int j=indx_min2;j<=indx_max2;j++){

      double res = 0.0;
      double res0 = 0.0;
      double res1 = 0.0;
      ii1 = i;  jj1 = j;  ii2 = i;  jj2 = j;

      //====================== Energy =====================
      if(i==j){
        res1 =  0.5 * (ela1->E_alp->M[i*ela1->Norb+i] + ela2->E_alp->M[i*ela2->Norb+i] );
      }

      
      for(a=0;a<ela1->Norb;a++){    // over basis of fragment A at time t
        for(b=0;b<elb1->Norb;b++){  // over basis of fragment B at time t



          //===================== Overlaps ===================
          /// Alpha channel contribution - interference can cause negative overlaps, so we need to exclude these terms
          /// but keep in mind that they exist
          /// (<A_i(t)|B_j(t+dt)> + <A_i(t+dt)|B_j(t)>)
///          res += (ela1->C_alp->M[a*ela1->Norb+ii1] * elb2->C_alp->M[b*elb2->Norb+jj2] * Sab0t->M[a*elb1->Norb+b]
///                + ela2->C_alp->M[a*ela2->Norb+ii2] * elb1->C_alp->M[b*elb1->Norb+jj1] * Sabt0->M[a*elb1->Norb+b]
///                 );                
 
          /// (<A_i(t+dt)|B_j(t+dt)> + <A_i(t)|B_j(t)>)
          res += (ela2->C_alp->M[a*ela2->Norb+ii2] * elb2->C_alp->M[b*elb2->Norb+jj2] * Sabtt->M[a*elb1->Norb+b]
                + ela1->C_alp->M[a*ela1->Norb+ii1] * elb1->C_alp->M[b*elb1->Norb+jj1] * Sab00->M[a*elb1->Norb+b]
                 );

          //======================== NAC =====================

          /// Alpha channel contribution
          /// Standard Tully-Hammess-Schiffer contribution:  
          /// d_AB_ij =  hbar*(<A_i(t)|B_j(t+dt)> - <A_i(t+dt)|B_j(t)>)/(2*dt)
          if(nac_opt==0 || nac_opt==1){

            res0 += (ela1->C_alp->M[a*ela1->Norb+ii1] * elb2->C_alp->M[b*elb2->Norb+jj2] * Sab0t->M[a*elb1->Norb+b]
                   - ela2->C_alp->M[a*ela2->Norb+ii2] * elb1->C_alp->M[b*elb1->Norb+jj1] * Sabt0->M[a*elb1->Norb+b]
                    );

          }
                 
          /// Less trivial contribution, due to non-orthogonality:
          /// hbar*(<A_i(t+dt)|B_j(t+dt)> - <A_i(t)|B_j(t)>)/(2*dt)
          /// This leads to non-anti-symmetric NAC matrix
          if(nac_opt==1){

            res0 += (ela2->C_alp->M[a*ela2->Norb+ii2] * elb2->C_alp->M[b*elb2->Norb+jj2] * Sabtt->M[a*elb1->Norb+b]
                   - ela1->C_alp->M[a*ela1->Norb+ii1] * elb1->C_alp->M[b*elb1->Norb+jj1] * Sab00->M[a*elb1->Norb+b]
                    );

          }


        }// for b
      }// for a

      ovlp->M[(i-indx_min1)*(indx_max2-indx_min2+1)+(j-indx_min2)] = (0.5/dt)*res;
       nac->M[(i-indx_min1)*(indx_max2-indx_min2+1)+(j-indx_min2)] = (0.5*hbar/dt)*res0;
       ene->M[(i-indx_min1)*(indx_max2-indx_min2+1)+(j-indx_min2)] = res1;




    }// for j
  }// for i

  /// Delete temporary matrices - this may be not very efficient, due to constant need for memory reallocation
  /// but it is simple
  delete Sab00;
  delete Sab0t;
  delete Sabt0;
  delete Sabtt;
  
}


void collect_matrices(int Nfrag, vector<MATRIX*>& x, vector<int>& min_orbs, vector<int>& max_orbs,MATRIX* X, vector<int>& frags,vector<int>& orbs){
// Matrix X should already be allocated!

//  cout<<"In collect_matrices\n";

  frags = vector<int>(Nfrag);
  orbs = vector<int>(Nfrag);
  *X = 0.0;

  int i = 0;
  for(int f1=0;f1<Nfrag;f1++){   
    for(int i1=min_orbs[f1];i1<=max_orbs[f1];i1++){
      frags[i] = f1;
      orbs[i] = i1;
      i++;
    }// for i1
  }// for f1

  int Xsize = i; // dimension of the X matrix

  int a, b;
  a = 0;
  for(int f1=0;f1<Nfrag;f1++){      // fragments
    for(int i1=min_orbs[f1];i1<=max_orbs[f1];i1++){ // all orbitals in the fragment f1

      b = 0;      
      for(int f2=0;f2<Nfrag;f2++){   // fragment B

        int x_i = f1*Nfrag + f2;               // index of the blok-matrix for fragments f1 and f2, assuming special order of matrices in x!!!
        int szf2 = max_orbs[f2]-min_orbs[f2]+1;// how many orbitals in the fragment f2

        for(int i2=min_orbs[f2];i2<=max_orbs[f2];i2++){

          X->M[a*Xsize+b] = x[x_i]->M[(i1-min_orbs[f1])*szf2+(i2-min_orbs[f2])];
          b++;

        }// for i2
      }// for f2

      a++;
    }// for i1
  }// for f1

//  cout<<"End of collect\n";


}

