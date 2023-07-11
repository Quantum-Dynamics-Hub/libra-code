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
  \file Basis_map.cpp
  \brief The file implements functions for mapping atom-orbital indices
    
*/

#include "Basis.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libqobjects;


/// libbasis namespace
namespace libbasis{


void show_mapping(const vector<vector<int> >& at_orbitals){
/** 
  \brief Printing the mapping information
  \param[in] at_orbital The mapping table: at_orbitals[n] is the list of the global indices of the AOs attached to the atom n
*/

  for(int n=0;n<at_orbitals.size();n++){
  
    cout<<"List n= "<<n<<" has following entries: ";
    for(int i=0;i<at_orbitals[n].size();i++){  cout<<at_orbitals[n][i]<<"  ";   }
    cout<<endl;

  }// for n

}


}//namespace libbasis
}//namespace liblibra


