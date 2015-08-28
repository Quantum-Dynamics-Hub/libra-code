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

#ifndef CONVERTERS_H
#define CONVERTERS_H

#include "../chemobjects/libchemobjects.h"
#include "../dyn/libdyn.h"

using namespace libchemobjects;
using namespace libdyn;


namespace libconverters{


void system_to_nuclear(System& syst, Nuclear& nucl);
void nuclear_to_system(Nuclear& nucl, System& syst);

void system_to_vector_q(System& syst, vector<double>& q);
void system_to_vector_p(System& syst, vector<double>& p);




}// namespace libconverters

#endif // CONVERTERS_H

