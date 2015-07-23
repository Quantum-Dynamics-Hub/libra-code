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

#ifndef HAMILTONIAN_ATOMISTIC_H
#define HAMILTONIAN_ATOMISTIC_H

#include "../Hamiltonian_Generic/Hamiltonian.h"
#include "Hamiltonian_MM/Hamiltonian_MM.h"


namespace libhamiltonian{
namespace libhamiltonian_atomistic{

using namespace libmmath;
using namespace libhamiltonian_generic;
using namespace libhamiltonian_mm;


class Hamiltonian_Atomistic : public Hamiltonian{

 
  
public:

  // Data members:
  listHamiltonian_MM*  mm_ham;   // mm part


  // Constructor: only allocates memory and sets up related variables
  Hamiltonian_Atomistic(int, int);

  // Destructor
  ~Hamiltonian_Atomistic();
/*
  // Set properties
  void set_rep(int rep_);

  // Set parameters
  void set_params(vector<double>& params_);
  void set_params(boost::python::list params_);
  void set_q(vector<double>& q_);
  void set_q(boost::python::list q_);
  void set_v(vector<double>& v_);
  void set_v(boost::python::list v_);

  // Perform actual computations - this will construct the internals of the object of this type
  void compute();
*/
  void compute_diabatic();
  void compute_adiabatic();




};

typedef std::vector<Hamiltonian_Atomistic> Hamiltonian_AtomisticList;


}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian

#endif // HAMILTONIAN_ATOMISTIC_H
