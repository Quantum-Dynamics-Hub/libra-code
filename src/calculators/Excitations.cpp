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
/**
  \file Excitations.cpp
  \brief The file implements functions for creation of the excitation objects for excited state calculations
    
*/

#include "Excitations.h"
#include "Bands.h"

/// libcalculators namespace
namespace libcalculators{

void excite(int I, int J, vector< pair<int,double> >& occ_ini, vector< pair<int,double> >& occ_fin){
/**
  \brief Create a list representing a I --> J ecxitation

  basically, vector< pair<int,double> > datastructure represents Slater determinant (by selecting ordering of orbitals in
  the active space)

  \param[in]  occ_ini  initial occupation list 
  \param[out] occ_fin  final occupation list
*/


  int Norb = occ_ini.size();

  if(I<0){ std::cout<<"Error: source orbital index("<<I<<") can not be negative\n";  exit(0); }
  if(J<0){ std::cout<<"Error: target orbital index("<<J<<") can not be negative\n";  exit(0); }
  if(I>=Norb){ std::cout<<"Error: source orbital index("<<I<<") exceeds the number of orbitals in the active space("<<Norb<<")\n";  exit(0); }
  if(J>=Norb){ std::cout<<"Error: target orbital index("<<J<<") exceeds the number of orbitals in the active space("<<Norb<<")\n";  exit(0); }


  if(occ_fin.size()>0){  occ_fin.clear(); }                    // clean result

  for(int i=0;i<Norb;i++){  occ_fin.push_back(occ_ini[i]);  }  // copy initial to final

  // Swap populations I-th and J-th orbitals
  double pop = occ_fin[I].second; 
  occ_fin[I].second = occ_fin[J].second;
  occ_fin[J].second = pop;
   
}

boost::python::list excite(int I, int J, boost::python::list occ_ini){
/**
  \brief Create a list representing a I --> J ecxitation (Python-friendly)

  basically, vector< pair<int,double> > datastructure represents Slater determinant (by selecting ordering of orbitals in
  the active space)
  Returns the final occupation list

  \param[in]  occ_ini  initial occupation list 
*/

  vector< pair<int,double> > occ_i;
  vector< pair<int,double> > occ_f;

  convert_1(occ_ini, occ_i);

  excite(I,J,occ_i, occ_f);

  return convert_2(occ_f);

}


}// namespace libcalculators
