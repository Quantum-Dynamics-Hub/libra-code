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
  \file Excitation.h
  \brief The file describes the class that stores the parameters controlling the calculations. Also the 
  auxiliary classes and functions are defined
    
*/

#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H


//#include <complex>
#include <vector>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>


//#include <boost/python.hpp>
//#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


/// liblibra namespace
namespace liblibra{

using namespace std;

/// libcontrol_parameters namespace
namespace libcommon_types{



class excitation{
/**
  The class representing the Slater determinant excitation
  Example:
  for convenience, in the input spin components is given by a letter: A = 1, B = -1
  { {0, 0}, {1, -1}, {1, 1}, {1, -1} } = double excitations from HOMO to LUMO with conservation of spin multiplicity
*/

  public:
  int size;                ///< 1 = single, 2 = double, etc
  vector<int> from_orbit;  ///< indices of the orbitals from which excitation is performed: -1 = HOMO-1, 0 = HOMO, 1 = LUMO, etc.
  vector<int> from_spin;   ///< same as above, but for spin: 1 = spin up, -1 = spin down
  vector<int> to_orbit;    ///< indices of the orbitals to which excitation is performed: -1 = HOMO-1, 0 = HOMO, 1 = LUMO, etc.
  vector<int> to_spin;     ///< same as above, but for spin: 1 = spin up, -1 = spin down

  double Energy;           ///< electronic energy of this configuration
  double f;                ///< oscillator strength


  /// default constructor
  excitation() { size = 0; Energy = 0.0; f = 0.0; }
  
  excitation(int _f_o, int _f_s, int _t_o, int _t_s){ 
  /** Constructor to generate single excitations  
  
  \param[in] _f_o "from orbital" - index of the source orbital
  \param[in] _f_s "from spin" - spin index (1 - alpha, -1 - beta) of the source orbital
  \param[in] _t_o "to orbital" - index of the target orbital
  \param[in] _t_s "to spin" - spin index (1 - alpha, -1 - beta) of the target orbital

  */
    size = 1;
    Energy = 0.0;
    f = 0.0;
    from_orbit = vector<int>(1,_f_o);
    from_spin = vector<int>(1,_f_s);
    to_orbit = vector<int>(1,_t_o);
    to_spin = vector<int>(1,_t_s);
    cout<<"Constructing single excitation: "<<_f_o<<" "<<_f_s<<" -> "<<_t_o<<" "<<_t_s<<endl;
  }



};

}// libcommon_types
}// liblibra

#endif //  COMMON_TYPES
