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
  \file SD.h
  \brief The file describes a Slater Determinant (SD) class that represents an antisimmetrized N-electron orbital
    
*/

#ifndef SD_H
#define SD_H

#include "../math_linalg/liblinalg.h"
#include "../math_meigen/libmeigen.h"
#include "../math_specialfunctions/libspecialfunctions.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libmeigen;
using namespace libspecialfunctions;

namespace libqobjects{


class SD{
/**
  Objects of this class are Slater Determinants
  Basically, it is a complex matrix of size Nbas x N,

  where N - is the number of active electrons (and spin-orbitals included in the determinant)
  and Nbas - is the number of basis states in which each spin-orbital is represented

  SD = det{ psi_1, psi_2, psi_3, ..., psi_N }  

  For convenience, the alpha and beta orbitals are sorted into two groups, so that
  first N_alp columns correspond to alpha-spin orbitals, and the following
  N_bet - to beta-spin orbitals. Of course N = N_alp + N_bet

  When considering excitations, keep in mind that the position of the orbitals involved in the
  excitation should remain unchanged - we only change them with different MO coefficients (corresponding to given excitation)

*/

public:

  CMATRIX* mo;      ///< The storage of the MOs included in the SD

  int N_bas;        ///< The number of basis functions used to represent each MO spin-orbital
  int N;            ///< The total number of electrons (spin-orbital) included in the wavefunction of SD type
  int N_alp;        ///< The number of alpa electrons (alpha spin-orbitals)
  int N_bet;        ///< The number of beta electrons (beta spin-orbitals)
  vector<int> orb_indx_alp;  ///< Indices of the spatial orbitals that enter the present determinant as 
                             ///< alpha spin-orbitals. The indexing is w.r.t. the external pool of alpha orbitals
                             ///< The length of this array corresponds to N_alp
  vector<int> orb_indx_bet;  ///< Indices of the spatial orbitals that enter the present determinant as 
                             ///< beta spin-orbitals. The indexing is w.r.t. the external pool of beta orbitals
                             ///< The length of this array corresponds to N_bet
  vector<int> spin; ///< The list of indices defining the spin of each orbital in this SD, in the same order as the MOs 
                    ///< are ordered in coloumns in the mo matrix. Values: 1 - alpha, -1 - beta, 0 - 2-component spinor
                    ///< The length of this array is N = N_alp + N_bet
                    ///< By default, we will be ordering the spin-orbitals such that we first list all alpha-spins, and then
                    ///< all beta-spins
  
  //----------------- Function members --------------------
  SD(); ///< default c-tor 
  SD(int _N_bas, int _N_alp,int _N_bet); ///< c-tor
  SD(CMATRIX& _mo_pool_alp, CMATRIX& _mo_pool_bet, vector<int>& _orb_indx_alp, vector<int>& _orb_indx_bet); ///< c-tor

  SD(const SD&); ///< cc-tor
  void operator=(const SD&); ///< assignment
  ~SD(); ///< destructor


  friend int operator == (const SD& sd1, const SD& sd2){
    int i;
    int res = ((sd1.N==sd2.N) && (sd1.N_bas==sd2.N_bas) && (sd1.N_alp==sd2.N_alp) && (sd1.N_bet==sd2.N_bet));
    if (res==1){
      for(i=0;i<sd1.N;i++){    res *= (sd1.spin[i]==sd2.spin[i]);   }
      for(i=0;i<sd1.N_alp;i++){    res *= (sd1.orb_indx_alp[i]==sd2.orb_indx_alp[i]);   }
      for(i=0;i<sd1.N_bet;i++){    res *= (sd1.orb_indx_bet[i]==sd2.orb_indx_bet[i]);   }      
    }
    return res;
  }

  CMATRIX get(){ return *mo; }
  void set(CMATRIX& _mo){ *mo = _mo; }  ///< Only update the actual MO data
  void set(CMATRIX& _mo_pool_alp, CMATRIX& _mo_pool_bet, vector<int>& _orb_indx_alp, vector<int>& _orb_indx_bet);
  double normalization_factor(){  return (1.0/sqrt(FACTORIAL(N))); }


};


complex<double> overlap_sd(CMATRIX& Smo, vector<int>& SD1, vector<int>& SD2);
CMATRIX overlap_sd(CMATRIX& Smo, vector< vector<int> >& SD_basis);

typedef std::vector<SD> SDList; ///< This is the data type for representing vector of SD objects


}/// namespace libqobjects
}/// namespace liblibra

#endif // SD_H
