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

#include "Engine.h"
#include "MOAO.h"
#include "Control_Parameters.h"
#include "Model_Parameters.h"
#include "Nuclear.h"
#include "Basis.h"
#include "Electronic.h"
#include "Hamiltonian.h"
#include "Memory.h"

#include <stdlib.h>
#include <time.h>
#include <sstream>

using namespace std;

void guess(Control_Parameters&, Model_Parameters&, Nuclear&,
           vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&, 
           Electronic*,Memory*);

