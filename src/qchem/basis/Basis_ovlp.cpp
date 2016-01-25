/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Basis_ovlp.cpp
  \brief The file implements function for updating the AO overlap matrix
    
*/

#include "Basis.h"


/// libqchem namespace
namespace libqchem{

/// libbasis namespace
namespace libbasis{



void update_overlap_matrix(int x_period,int y_period,int z_period,const VECTOR& t1, const VECTOR& t2, const VECTOR& t3,
                           vector<AO>& basis_ao, MATRIX& Sao){
/**
  \brief Update the oberlap matrix (in AO basis): <AO(i)|AO(j)>
  \param[in] x_period Then number of periodic shells in X direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] y_period Then number of periodic shells in Y direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] z_period Then number of periodic shells in Z direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] t1 The periodicity vector along a crystal direction ("X")
  \param[in] t2 The periodicity vector along b crystal direction ("Y")
  \param[in] t3 The periodicity vector along c crystal direction ("Z")
  \param[in] basis_ao The list of all AOs (basis)
  \param[out] Sao The output overlap matrix

  This function can also take periodic images of the system into account
*/


  int i,j,n,I,J;
  VECTOR dIdA,dIdB,TV, Rij,Rij0;
  double dist, dist_min;  
  int opt_x,opt_y,opt_z;

  int Norb = basis_ao.size();

  Sao = 0.0;
 
  for(i=0;i<Norb;i++){
    for(j=i;j<Norb;j++){      

      for(int nx=-x_period;nx<=x_period;nx++){
        for(int ny=-y_period;ny<=y_period;ny++){
          for(int nz=-z_period;nz<=z_period;nz++){

            // This summation corresponds to k = 0 (Gamma-point)    
            TV = nx*t1 + ny*t2 + nz*t3;

            if(i==j){
              AO tmp_ao(basis_ao[i]); 
              tmp_ao.shift_position(TV);
              Sao.M[i*Norb+j] += gaussian_overlap(basis_ao[i],tmp_ao); //,0,dIdA,dIdB,aux,n_aux,TV);
            }
            else{

              basis_ao[j].shift_position(TV);
              Sao.M[i*Norb+j] += gaussian_overlap(basis_ao[i],basis_ao[j]); //,0,dIdA,dIdB,aux,n_aux,TV);
              basis_ao[j].shift_position(-TV);

            }// i != j

          }// for nz
        }// for ny
      }// for nx 

      Sao.M[j*Norb+i] = Sao.M[i*Norb+j];

    }// j
  }// i



/*   !!!!!!!!!!!!!! Don't need the below part untill we get to the gradients !!!!!!!!!!!

  else{

    for(n=0;n<Nnucl;n++){
      *dSao_dx[n] = 0.0;
      *dSao_dy[n] = 0.0;
      *dSao_dz[n] = 0.0;
    }

  
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){

        Sao->M[i*Norb+j] = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);

        dSao_dx[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.x;
        dSao_dx[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.x;

        dSao_dy[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.y;
        dSao_dy[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.y;

        dSao_dz[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.z;
        dSao_dz[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.z;
*/
/*
      // Tests
      double dx = 0.0001;
      R[basis_ao[i].at_indx] += VECTOR(dx,0.0,0.0);
      double s2 = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);
      R[basis_ao[i].at_indx] -= VECTOR(2.0*dx,0.0,0.0);
      double s1 = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);
      R[basis_ao[i].at_indx] += VECTOR(dx,0.0,0.0);
      OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);

      cout<<"Overlap:  der(numer)= "<<0.5*(s2-s1)/dx<<" der(anal)= "<<dIdA.x<<endl;
*/
     

//    }// for j
//  }// for i
//  }// else : method!=eht0

//  for(int n=0;n<Nnucl;n++){
//    cout<<"dSao_dx["<<n<<"]= \n"<<*dSao_dx[n]<<endl;
//  }

//  cout<<"Overlap matrix in AO basis:\n";
//  cout<<*Sao<<endl;
  // So far use my inverter:
//  Sao->Inverse(Sao_inv);

}





}//namespace libbasis
}//namespace libqchem
