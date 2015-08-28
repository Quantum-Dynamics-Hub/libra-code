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

#ifndef DIPOLE_H
#define DIPOLE_H

#include "MOAO.h"
#include "Nuclear.h"
#include "Electronic.h"
#include "Basis.h"
#include "Memory.h"

#include <stdlib.h>
#include <time.h>
#include <sstream>

using namespace std;

void update_dipole_matrix(int,int,int,const VECTOR&,const VECTOR&,const VECTOR&,
                          vector<int>&,vector<AO>&,
                          MATRIX*,MATRIX*,MATRIX*,vector<double*>&,int,Nuclear&);

void compute_dipole_moments(Control_Parameters&, Model_Parameters&, Nuclear&,
                            vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&,
                            Electronic*, Memory*,  VECTOR&,MATRIX*,MATRIX*,MATRIX*);



#endif // DIPOLE_H
