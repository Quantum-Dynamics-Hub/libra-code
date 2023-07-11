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




int ForceField::set_ff_epsilon_and_sigma(int sz, vector<string> types,double** epsilon, double** sigma){

  int res = 1;
  for(int i=0;i<sz;i++){
    //-------------- Start with looking the whole record ---------------------
    // Find index of Atom_Record corresponding to the force field type types[i]
    int ff_atom1_indx = Atom_Record_Index(types[i]);
    if(ff_atom1_indx==-1){
      *epsilon[i] = 0.0; *sigma[i] = 3.0 * Angst;
      cout<<"Warning: In ForceField::set_ff_epsilon_and_sigma : Can not fing the atom type ff_atom type in force field "<<ForceField_Name<<"\n";
      cout<<"Setting epsilon to zero and sigma to 3.0\n";
    }
    else{
      map<string,double> d;
      get_vdw_parameters(types[i],types[i],"none",d);
      *epsilon[i] = d["epsilon"];
      *sigma[i] = d["sigma"];
    }
  }// for i

  return res;

}


}// namespace libforcefield
}// namespace liblibra

