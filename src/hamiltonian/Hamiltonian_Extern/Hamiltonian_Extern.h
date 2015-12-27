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

#ifndef HAMILTONIAN_EXTERN_H
#define HAMILTONIAN_EXTERN_H

#include "../Hamiltonian_Generic/Hamiltonian.h"



namespace libhamiltonian{

using namespace libhamiltonian_generic;
using namespace libmmath;

namespace libhamiltonian_extern{


class Hamiltonian_Extern : public Hamiltonian{

  //  We will use the Hamiltonian's members directly (see in Hamiltonian.h )

  // Diabatic representation:
  // MATRIX* ham_dia;           // Hamiltonian in diabatic representation
  // vector<MATRIX*> d1ham_dia; // derivatives of the ham_dia w.r.t. all atomic DOFs: q0, q1, .. qN 
  // vector<MATRIX*> d2ham_dia; // derivatives of the ham_dia w.r.t. all atomic DOFs: q00, q01, ..., q0N, q10, q11, ... qNN

  // Adiabatic representation:
  // MATRIX* ham_adi;           // Hamiltonian in adiabatic representation
  // vector<MATRIX*> d1ham_adi; // first order derivative couplings: <i|d/dR|j> - is computed from the transformation coefficients


  int adiabatic_opt; // defines how to perform adiabatic calculations: either use bound adiabatic matrices (0, default), or to
                     // perform diabatic -> adiabatic transformation (1)

  // Bind status (bs_)
  int bs_ham_dia;
  int bs_d1ham_dia;
  int bs_d2ham_dia;
  int bs_ham_adi;
  int bs_d1ham_adi;

  
   
public:

  // Constructor: only allocates memory and sets up related variables
  Hamiltonian_Extern(int _nelec, int _nnucl);

  // Destructor
  ~Hamiltonian_Extern();

  // Set parameters
  void set_adiabatic_opt(int);

  void bind_ham_dia(MATRIX& _ham_dia);
  void bind_d1ham_dia(vector<MATRIX>& _d1ham_dia);
  void bind_d2ham_dia(vector<MATRIX>& _d2ham_dia);
  void bind_ham_adi(MATRIX& _ham_adi);
  void bind_d1ham_adi(vector<MATRIX>& _d1ham_adi);


  // Perform actual computations - this will construct the internals of the object of this type
  void compute_diabatic();
  void compute_adiabatic();


};

typedef std::vector<Hamiltonian_Extern> Hamiltonian_ExternList;


}// namespace libhamiltonian_extern
}// namespace libmodel

#endif // HAMILTONIAN_EXTERN_H
