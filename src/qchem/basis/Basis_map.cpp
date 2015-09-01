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

#include "Basis.h"

/****************************************************************************

  This file contains following functions:

  void map_atoms_and_orbitals(int Nnucl, const vector<AO>&  basis_ao, vector<vector<int> >& at_orbitals)
  void show_mapping(const vector<vector<int> >& at_orbitals)

****************************************************************************/

namespace libqchem{
namespace libbasis{


void show_mapping(const vector<vector<int> >& at_orbitals){

  for(int n=0;n<at_orbitals.size();n++){
  
    cout<<"List n= "<<n<<" has following entries: ";
    for(int i=0;i<at_orbitals[n].size();i++){  cout<<at_orbitals[n][i]<<"  ";   }
    cout<<endl;

  }// for n

}


}//namespace libbasis
}//namespace libqchem


