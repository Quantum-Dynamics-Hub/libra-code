/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Verlet.cpp
  \brief The file implements the velocity Verlet algorithm for nuclear dynamics integration
    
*/

#include "Dynamics.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;
using namespace libhamiltonian;
namespace bp = boost::python;


/// libdyn namespace 
namespace libdyn{



void Verlet(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params){
/**
  \brief One step of Velocity verlet algorithm for nuclear DOF
  \param[in] Integration time step
  \param[in,out] mol Describes the nuclear DOF. Changes during the integration.
  \param[in] el Describes electronic DOF. Does not change during the integration.
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

*/

  p = p + ham.forces_adi().real() * 0.5*dt;
  q = q + invM*p*dt;

  ham.compute_diabatic(py_funct, bp::object(q), params);
  ham.compute_adiabatic(1);

  p = p + ham.forces_adi().real() * 0.5*dt;


}// Verlet



}// namespace libdyn
}// liblibra


