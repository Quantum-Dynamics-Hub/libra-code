/*********************************************************************************
* Copyright (C) 2016-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file SD.cpp
  \brief The file implements the SD class (for representing Slater Determinants) and its members
    
*/

#include "SD.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libmeigen;
using namespace libspecialfunctions;


/// libqobjects namespace
namespace libqobjects{


SD::SD(){ 
/**
  \brif The default constructor.
*/

  mo = NULL; 
  N_bas = 0; 
  N = 0;     
  N_alp = 0; 
  N_bet = 0; 


}

SD::SD(int _N_bas, int _N_alp,int _N_bet){
/**
  \brif The constructor.that allocates memory 
  \param[in] _N_bas The number of basis functions for each MO
  \param[in] _N_alp The number of alpha spin-orbitals (first _N_alp coloumns)
  \param[in] _N_bet The number of beta spin-orbitals (last _N_bet coloumns)
*/

  N_bas = _N_bas;
  N_alp = _N_alp;
  N_bet = _N_bet;
  N = N_alp + N_bet;
  mo = new CMATRIX(N_bas, N);

  orb_indx_alp = vector<int>(N_alp);
  orb_indx_bet = vector<int>(N_bet);
  spin = vector<int>(N);

  for(int i=0;i<N_alp;i++){  orb_indx_alp[i] = i;   spin[i] = 1; }
  for(int i=0;i<N_bet;i++){  orb_indx_bet[i] = i;   spin[N_alp+i] = -1; }

}

SD::SD(CMATRIX& _mo_pool_alp, CMATRIX& _mo_pool_bet, vector<int>& _orb_indx_alp, vector<int>& _orb_indx_bet){

  set(_mo_pool_alp, _mo_pool_bet, _orb_indx_alp, _orb_indx_bet);

}


void SD::set(CMATRIX& _mo_pool_alp, CMATRIX& _mo_pool_bet, vector<int>& _orb_indx_alp, vector<int>& _orb_indx_bet){
/**
  \brif The constructor.that allocates memory and assignes values
  \param[in] _mo_pool_alp A set of all MO orbitals (spatial components) for alpha spin channel
  \param[in] _mo_pool_bet A set of all MO orbitals (spatial components) for beta spin channel
  \param[in] _orb_indx_alp The list of orbital global MO pool indices that for the MOs that should be included into present determinant
  as the spin-up orbitals. The length of this vector will define N_alp
  \param[in] _orb_indx_bet The list of orbital global MO pool indices that for the MOs that should be included into present determinant
  as the spin-down orbitals. The length of this vector will define N_bet

  Example:
  
  Consider two configurations (U - electron up, D - electron down)

     PHI_0          PHI_1

  3 --------     --------
  2 --------     --  D --
  1 --U D --     --U   --
  0 --U D --     --U D --

  The pools of orbitals (_mo_pool_alp) and (_mo_pool_bet) can be   Nbas x 4 matrices (including all 4 orbitals)
  The _orb_indx_alp and _orb_indx_bet variables defining the configurations should then be:
            
                     PHI_0      PHI_1

  _orb_indx_alp    [ 0, 1 ]    [ 0, 1 ] 
  _orb_indx_bet    [ 0, 1 ]    [ 0, 2 ]

*/

  if(_mo_pool_alp.n_rows!=_mo_pool_bet.n_rows){
    std::cout<<"Error in SD::SD : The number of rows in the _mo_pool_alp matrix ("<<_mo_pool_alp.n_rows
             <<") is not equal to the number of rows in the _mo_pool_bet matrix ("<<_mo_pool_bet.n_rows<<" )\n";
    exit(0);
  }
  N_bas = _mo_pool_alp.n_rows;

  orb_indx_alp = _orb_indx_alp;
  orb_indx_bet = _orb_indx_bet;

  N_alp = _orb_indx_alp.size();
  N_bet = _orb_indx_bet.size();
  N = N_alp + N_bet;

  mo = new CMATRIX(N_bas, N);
  spin = vector<int>(N);

  int i,j;
  for(i=0;i<N_alp;i++){  
    if(orb_indx_alp[i]>=_mo_pool_alp.n_cols){
      std::cout<<"Index of an alpha orbital is beyond the range of MO orbital indices available in the pool provided\n";
      std::cout<<"Index of an alpha orbital = "<<orb_indx_alp[i]<<endl;
      std::cout<<"Max index available in the pool provided = "<<_mo_pool_alp.n_cols<<endl;
      exit(0);
    }

    for(j=0;j<N_bas;j++){  mo->set(j,i, _mo_pool_alp.get(j,orb_indx_alp[i]));    }
    spin[i] = 1; 
  }
  for(i=0;i<N_bet;i++){  
    if(orb_indx_bet[i]>=_mo_pool_bet.n_cols){
      std::cout<<"Index of a beta orbital is beyond the range of MO orbital indices available in the pool provided\n";
      std::cout<<"Index of a beta orbital = "<<orb_indx_bet[i]<<endl;
      std::cout<<"Max index available in the pool provided = "<<_mo_pool_bet.n_cols<<endl;
      exit(0);
    }

    for(j=0;j<N_bas;j++){  mo->set(j,N_alp+i, _mo_pool_bet.get(j,orb_indx_bet[i]));    }
    spin[N_alp + i] = -1; 
  }

}

SD::SD(const SD& sd){ 
/** \brief Copy constructor
*/


  N_bas = sd.N_bas;
  N     = sd.N;
  N_alp = sd.N_alp;
  N_bet = sd.N_bet;
    
  orb_indx_alp = sd.orb_indx_alp; 
  orb_indx_bet = sd.orb_indx_bet; 

  spin = sd.spin;

  mo = new CMATRIX(N_bas, N);
  *mo = *sd.mo;

}

void SD::operator=(const SD& sd){ 
/** \brief Assignment operator
*/

  N_bas = sd.N_bas;
  N     = sd.N;
  N_alp = sd.N_alp;
  N_bet = sd.N_bet;
    
  orb_indx_alp = sd.orb_indx_alp; 
  orb_indx_bet = sd.orb_indx_bet; 

  spin = sd.spin;

  *mo = *sd.mo;
  
//  return *this; // return reference to allow chaining: A = B = C =... !!! No: so crap doesn't happen in PYthon

}

SD::~SD(){
/** \brief Destructor
*/
  if(orb_indx_alp.size()>0){  orb_indx_alp.clear(); }  
  if(orb_indx_bet.size()>0){  orb_indx_bet.clear(); }

  if(spin.size()>0){  spin.clear(); }  

  delete mo;
  mo = NULL;

  N_bas = 0; 
  N = 0;     
  N_alp = 0; 
  N_bet = 0; 


}




complex<double> overlap_sd(CMATRIX& Smo, vector<int>& SD1, vector<int>& SD2){
/**
##
# \brief Compute the overlap of the normalized determinants SD1 and SD2
#
# \param[in] Smo A matrix of spin-orbital overlaps: Smo_ij = <psi_i(t)|psi_j(t')> 
# \param[in] SD1 represents a determinant SD1 = 1/sqrt(N!) * |psi_i1(1) psi_i2(2) psi_i3(3) ... psi_iN(N)| - at time t
# \param[in] SD2 represents a determinant SD2 = 1/sqrt(N!) * |psi_i1(1) psi_i2(2) psi_i3(3) ... psi_iN(N)| - at time t'
# 
#  Note #1: the psi_ functions in the two sets are not orthogonal to each other, so
#  Smo_ij is not delta_ij  !!!
#
#  Note #2: the psi_ functions are the spin-orbitals, not just the spatial components
#  they are enumerated by integers [0, 1, 2, etc.] the meaning of these numbers are to be 
#  defined by the one computing the Smo matrix, but the indices should be consistent with
#  the one used in Smo
#
*/

    if(SD1.size() != SD2.size()){
      std::cout<<"Can not compute an overlap of Slater determinants of different size\n";
      exit(0);
    }

    int sz = SD1.size();

    CMATRIX smo(sz, sz);
    pop_submatrix(Smo, smo, SD1, SD2);
    
    complex<double> ovlp = det(smo); 

    return ovlp;  // this is a complex number, in general

}

CMATRIX overlap_sd(CMATRIX& Smo, vector< vector<int> >& SD_basis){

    int bas_sz = SD_basis.size();

    CMATRIX res(bas_sz, bas_sz);
    for(int i=0;i<bas_sz;i++){
      for(int j=0;j<bas_sz;j++){

        res.set(i,j, overlap_sd(Smo, SD_basis[i], SD_basis[j]));

      }// for j
    }// for i

    return res;
}



}/// libqobjects
}/// liblibra

