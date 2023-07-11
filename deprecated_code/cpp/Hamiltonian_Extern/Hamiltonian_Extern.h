/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Hamiltonian_Extern.h
  \brief The file describes the external Hamiltonian class - for interface with 3-rd party codes
    
*/

#ifndef HAMILTONIAN_EXTERN_H
#define HAMILTONIAN_EXTERN_H

#include "../Hamiltonian_Generic/Hamiltonian.h"
#include "../../math_linalg/liblinalg.h"
#include "../../math_meigen/libmeigen.h"


/// liblibra namespace
namespace liblibra{


/// libhamiltonian namespace
namespace libhamiltonian{

using namespace libhamiltonian_generic;
using namespace liblinalg;
using namespace libmeigen;

/// libhamiltonian_extern namespace
namespace libhamiltonian_extern{


class Hamiltonian_Extern : public Hamiltonian{
/**
  This is a class derived from the generic Hamiltonian class, so it inherits many of its properties
  //  We will use the Hamiltonian's members directly (see in Hamiltonian.h )
*/

  int adiabatic_opt; ///< defines how to perform adiabatic calculations: either use bound adiabatic matrices (0, default), or to
                     ///< perform diabatic -> adiabatic transformation (1)
  int vibronic_opt;  ///< defines how to perform vibronic Hamiltonian calculations: either use bound vibronic matrix (0, default), or to
                     ///< perform honest computations using electronic Hamiltonian and derivative couplings (1)


  // Bind status (bs_)
  int bs_ham_dia;   ///< bind status of diabatic Hamiltonian
  int bs_d1ham_dia; ///< bind status of first-order derivatives of the diabatic Hamiltonian
  int bs_d2ham_dia; ///< bind status of second-order derivatives of the diabatic Hamiltonian
  int bs_ham_adi;   ///< bind status of adiabatic Hamiltonian                                
  int bs_d1ham_adi; ///< bind status of first-order derivatives of the adiabatic Hamiltonian                     
  int bs_ham_vib;   ///< bind status of the vibronic Hamiltonian

  CMATRIX* ham_vib;  ///< Vibronic Hamiltonian 
  
   
public:

  /// Constructor: only allocates memory and sets up related variables
  Hamiltonian_Extern(int _nelec, int _nnucl);

  /// Destructor
  ~Hamiltonian_Extern();

  /// Set parameters
  void set_adiabatic_opt(int);
  void set_vibronic_opt(int vib_opt);

  void bind_ham_dia(MATRIX& _ham_dia);
  void bind_d1ham_dia(vector<MATRIX>& _d1ham_dia);
  void bind_d2ham_dia(vector<MATRIX>& _d2ham_dia);
  void bind_ham_adi(MATRIX& _ham_adi);
  void bind_d1ham_adi(vector<MATRIX>& _d1ham_adi);
  void bind_ham_vib(CMATRIX& _ham_vib);


  // Perform actual computations - this will construct the internals of the object of this type
  void compute_diabatic();
  void compute_adiabatic();

  std::complex<double> Hvib(int i,int j);

};

typedef std::vector<Hamiltonian_Extern> Hamiltonian_ExternList;  /// data type for keeping a list of external Hamiltonians of their derived classes


}// namespace libhamiltonian_extern
}// namespace libmodel
}// liblibra

#endif // HAMILTONIAN_EXTERN_H
