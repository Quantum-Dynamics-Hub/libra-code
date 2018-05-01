/*********************************************************************************
* Copyright (C) 2015-2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file tsh_methods_ida.cpp
  \brief The file implements the Instantaneous Decoherence at Attempted hops (ID-A)
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{



int ida(CMATRIX& Coeff, int old_st, int new_st, double E_old, double E_new, double T, double ksi){
/**
  \brief Instantaneous decoherence at attempted hops (ID-A)

  \param[in,out] Coeff The matrix of coefficients (amplitudes of the basis excited states). This matrix is modified at every
  decoherence event
  \param[in] old_st The old state before an attempted hop
  \param[in] new_st The new state after an attempted hop (no matter what the algorithm - FSSH, GFSH, MSSH, or something else)
  \param[in] E_old The energy corresponding to the old state
  \param[in] E_new The energy corresponding to the new state
  \param[in] T is the temperature of the system
  \param[in] ksi a random number from a uniform distribution of the [0,1] interval - needed to decide the hop acceptance/decoherence

  The function returns the index of the new state after hop rejection/decoherence criteria
  The function also modifies the amplitudes of the coherent superposition

*/

  const double kb = 3.166811429e-6;  // Boltzmann constant: Hartree/K
  int istate = old_st;               // by default the resulting state is the old state
  
  if(new_st != old_st){  // attempted hop

    // Now apply energy considerations
    double dE = (E_new - E_old);
    double boltz_f = 1.0;

    if(dE>0.0){
      double argg = dE/(kb*T);
      if(argg > 50.0){ boltz_f = 0.0; }
      else{            boltz_f = exp(-argg); }

      if(ksi<boltz_f){
        istate = new_st;  // accepted hop

        // Collapse the wavefunction onto the new state
        Coeff *= 0.0;
        Coeff.set(new_st, 1.0, 0.0);
      }
      else{
        // Unsuccessful hop - collapse wfc back to the original state
        Coeff *= 0.0;
        Coeff.set(old_st, 1.0, 0.0);
      }
    }// dE > 0
  }// new_st != old_st

  return istate;

}


}// namespace libdyn
}// liblibra

