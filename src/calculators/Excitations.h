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
  \file Excitations.h
  \brief The file defines functions for creation of the objects for excited state calculations
    
*/

#ifndef EXCITATIONS_H
#define EXCITATIONS_H

#include "../math_linalg/liblinalg.h"
#include "../common_types/libcommon_types.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libcommon_types;

/// libcalculators namespace
namespace libcalculators{

void excite(int I, int J, vector< pair<int,double> >& occ_ini, vector< pair<int,double> >& occ_fin);
boost::python::list excite(int I, int J, boost::python::list occ_ini);

void excite(int Norb, excitation& ex, 
            int Nocc_alp, vector< pair<int,double> >& occ_alp,
            int Nocc_bet, vector< pair<int,double> >& occ_bet);


}// libcalculators
}// liblibra

#endif // EXCITATIONS_H
