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

#include "converters.h"

/// liblibra namespace
namespace liblibra{


namespace libconverters{


void system_to_nuclear(System& syst, Nuclear& nucl){

  // Simple data converter - no object creation (nor reduction of dimension)
  if(nucl.nnucl != 3.0 * syst.Number_of_atoms){
    cout<<"Error: in system_to_nuclear: dimensions of System and Nuclear objects are not compatible\n";  
    exit(0);
  }
  

  for(int i=0;i<syst.Number_of_atoms;i++){

    nucl.mass[3*i+0] = syst.Atoms[i].Atom_RB.rb_mass; 
    nucl.mass[3*i+1] = syst.Atoms[i].Atom_RB.rb_mass; 
    nucl.mass[3*i+2] = syst.Atoms[i].Atom_RB.rb_mass; 

    nucl.q[3*i+0] = syst.Atoms[i].Atom_RB.rb_cm.x; 
    nucl.q[3*i+1] = syst.Atoms[i].Atom_RB.rb_cm.y; 
    nucl.q[3*i+2] = syst.Atoms[i].Atom_RB.rb_cm.z; 

    nucl.p[3*i+0] = syst.Atoms[i].Atom_RB.rb_p.x; 
    nucl.p[3*i+1] = syst.Atoms[i].Atom_RB.rb_p.y; 
    nucl.p[3*i+2] = syst.Atoms[i].Atom_RB.rb_p.z; 

    nucl.f[3*i+0] = syst.Atoms[i].Atom_RB.rb_force.x; 
    nucl.f[3*i+1] = syst.Atoms[i].Atom_RB.rb_force.y; 
    nucl.f[3*i+2] = syst.Atoms[i].Atom_RB.rb_force.z; 

/*
    nucl.f[3*i+0] = 0;
    nucl.f[3*i+1] = 1;
    nucl.f[3*i+2] = 2;
*/

  }// for i


}


void nuclear_to_system(Nuclear& nucl,System& syst){
  // Simple data converter - no object creation (nor reduction of dimension)

  if(nucl.nnucl != 3.0 * syst.Number_of_atoms){
    cout<<"Error: in system_to_nuclear: dimensions of System and Nuclear objects are not compatible\n";  
    exit(0);
  }
  

  for(int i=0;i<syst.Number_of_atoms;i++){

    syst.Atoms[i].Atom_RB.rb_mass = nucl.mass[3*i+0];

    syst.Atoms[i].Atom_RB.rb_cm.x = nucl.q[3*i+0]; 
    syst.Atoms[i].Atom_RB.rb_cm.y = nucl.q[3*i+1]; 
    syst.Atoms[i].Atom_RB.rb_cm.z = nucl.q[3*i+2];

    syst.Atoms[i].Atom_RB.rb_p.x = nucl.p[3*i+0]; 
    syst.Atoms[i].Atom_RB.rb_p.y = nucl.p[3*i+1]; 
    syst.Atoms[i].Atom_RB.rb_p.z = nucl.p[3*i+2]; 

    syst.Atoms[i].Atom_RB.rb_force.x = nucl.f[3*i+0];
    syst.Atoms[i].Atom_RB.rb_force.y = nucl.f[3*i+1];
    syst.Atoms[i].Atom_RB.rb_force.z = nucl.f[3*i+2];

  }// for i
}


void IndexError() { PyErr_SetString(PyExc_IndexError, "Index out of range"); }
void KeyError() { PyErr_SetString(PyExc_KeyError, "Key not found"); }

//boost::python::dict map_to_dict(std::map<std::string, double> map_){  return <std::string, double>map_to_dict_templ(map_); }




}// namespace libconverters
}// liblibra


