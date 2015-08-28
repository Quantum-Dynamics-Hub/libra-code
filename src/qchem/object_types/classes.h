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

#ifndef classes_h
#define classes_h


#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <complex>
#include <string.h>
#include "Mathematics.h"
#include "units.h"
using namespace std;



class Atom{
  public:

  std::string atom_type;    // type name of this atom
  VECTOR R;                 // coordinates of this atom

  // Constructor 
  Atom(){ ;; }
  Atom(std::string type_name_,VECTOR& r){
    atom_type = type_name_;
    R = r;
  }
  Atom(std::string type_name_,double x,double y,double z){ 
    atom_type = type_name_;
    R = VECTOR(x,y,z);
  }

};

class Cell{
  public:

  double a, b, c;          // lattice parameters
  double alp,bet,gam;      // lattice angles

  VECTOR tv1,tv2,tv3;      // lattice vectors

  // Constructor
  Cell(){ ;; }
  Cell(VECTOR& t1,VECTOR& t2,VECTOR& t3){
    tv1 = t1; tv2 = t2; tv3 = t3;
    a = tv1.length(); 
    b = tv2.length();
    c = tv3.length();
  }

};


class excitation{

  public:
  int size;                // 1 = single, 2 = double, etc
  vector<int> from_orbit;  // indices of the orbitals from which excitation is performed: -1 = HOMO-1, 0 = HOMO, 1 = LUMO, etc.
  vector<int> from_spin;   // same as above, but for spin: 1 = spin up, -1 = spin down
  vector<int> to_orbit;    // indices of the orbitals to which excitation is performed: -1 = HOMO-1, 0 = HOMO, 1 = LUMO, etc.
  vector<int> to_spin;     // same as above, but for spin: 1 = spin up, -1 = spin down

  double Energy;           // electronic energy of this configuration
  double f;                // oscillator strength

  // Example of exitation:
  // for convenience, in the input spin components is given by a letter: A = 1, B = -1
  // { {0, 0}, {1, -1}, {1, 1}, {1, -1} } = double excitations from HOMO to LUMO with conservation of spin multiplicity

  // dafault constructor
  excitation() { size = 0; Energy = 0.0; f = 0.0; }
  
  // single excitations
  excitation(int _f_o, int _f_s, int _t_o, int _t_s){ 
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




#endif // classes_h