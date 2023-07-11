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

#include "ForceField.h"

/// liblibra namespace
namespace liblibra{


namespace libforcefield{



int ForceField::set_ff_charges(int sz, vector<string> types, VECTOR** r, double** q){

  vector<int> indxs; 
  int res = 1;

  for(int i=0;i<sz;i++){
    //-------------- Start with looking the whole record ---------------------
    // Find index of Atom_Record corresponding to the force field type types[i]
    int ff_atom1_indx = Atom_Record_Index(types[i]);
    if(ff_atom1_indx==-1){ 
      *q[i] = 0.0;
      cout<<"Warning: In ForceField::set_ff_charges : Can not fing the atom type ff_atom type in force field "<<ForceField_Name<<"\n";
      cout<<"Setting charge to zero\n";
    }
    else{ 
      if(Atom_Records[ff_atom1_indx].is_Atom_partial_charge){ *q[i] = Atom_Records[ff_atom1_indx].Atom_partial_charge;   }
    }
  }// for i

  return res;
}


}// namespace libforcefield
}// namespace liblibra

