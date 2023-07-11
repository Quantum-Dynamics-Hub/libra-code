/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Basis_nac.cpp
  \brief The file implements function for creating derivative coupling matrix
    
*/

#include "Basis.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libqobjects;

/// libbasis namespace
namespace libbasis{



void update_derivative_coupling_matrix
(int x_period,int y_period,int z_period,const VECTOR& t1, const VECTOR& t2, const VECTOR& t3,
 vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
 vector<AO>& basis_ao, int c, MATRIX& Dao_x, MATRIX& Dao_y, MATRIX& Dao_z
){
/**
  \brief Update the derivative coupling matrix (in AO basis): <AO(i)|d/dR_c|AO(j)>
  \param[in] x_period Then number of periodic shells in X direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] y_period Then number of periodic shells in Y direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] z_period Then number of periodic shells in Z direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] t1 The periodicity vector along a crystal direction ("X")
  \param[in] t2 The periodicity vector along b crystal direction ("Y")
  \param[in] t3 The periodicity vector along c crystal direction ("Z")
  \param[in] atom_to_ao_map The mapping from atomic indices to AO indices: atom_to_ao_map[n][i] - is the global index of i-th AO on the n-th atom
  \param[in] ao_to_atom_map The mapping from the global AO index to the atomic index: ao_to_atom_map[I] - is the index of the atom on which 
  AO with the global index I is located
  \param[in] basis_ao The list of all AOs (basis)
  \param[in] c The index of the nuclear DOF (atom index) to which the electronic states are coupled - only this derivative matrix is considered here
  \param[out] Dao_x The derivative coupling matrix in AO basis, coupled to x coordinate of the atom with index c
  \param[out] Dao_y The derivative coupling matrix in AO basis, coupled to y coordinate of the atom with index c
  \param[out] Dao_z The derivative coupling matrix in AO basis, coupled to z coordinate of the atom with index c

  This function also takes periodic images of the system into account
*/

  int i,j,n,I,J;
  VECTOR dIdA,dIdB,TV, Rij,Rij0;
  double dist, dist_min;  
  int opt_x,opt_y,opt_z;

  int Norb = basis_ao.size();

  Dao_x = 0.0;
  Dao_y = 0.0;
  Dao_z = 0.0;
 
  for(j=0;j<Norb;j++){      

    if(c==ao_to_atom_map[j]){

      for(i=0;i<Norb;i++){

        VECTOR Dao; Dao = 0.0; // sum from all periodic contributions

        for(int nx=-x_period;nx<=x_period;nx++){
          for(int ny=-y_period;ny<=y_period;ny++){
            for(int nz=-z_period;nz<=z_period;nz++){
        
              // This summation corresponds to k = 0 (Gamma-point)    
              TV = nx*t1 + ny*t2 + nz*t3;
        
              //if(i==j){      }
              if(ao_to_atom_map[i]==ao_to_atom_map[j]){  } // orbitals on the same atom - no NAC
              else{
        
                basis_ao[j].shift_position(TV);
        
                VECTOR dao;
                dao = derivative_coupling_integral(&basis_ao[i], &basis_ao[j]); //<i|d/dR_a|j> where a - is the atom of orbital j
                Dao += dao;
        
                basis_ao[j].shift_position(-TV);
        
              }// i != j
        
            }// for nz
          }// for ny
        }// for nx 

        Dao_x.M[i*Norb+j] = Dao.x;
        Dao_y.M[i*Norb+j] = Dao.y;
        Dao_z.M[i*Norb+j] = Dao.z;

      }// for i
    }// if c==ao_to_atom_map[j]
  }// for j




}




}//namespace libbasis
}//namespace liblibra
