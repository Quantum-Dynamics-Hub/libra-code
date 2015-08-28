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

#ifndef SCF_H
#define SCF_H

#include "MOAO.h"
#include "Control_Parameters.h"
#include "Model_Parameters.h"
#include "Electronic.h"
#include "Nuclear.h"
#include "Memory.h"


double scf(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&, vector<int>&,vector<AO>&,vector<vector<int> >&, Electronic*,Electronic*, Memory*);

double scf_oda(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&, vector<int>&,vector<AO>&,vector<vector<int> >&, Electronic*,Electronic*, Memory*);
double scf_oda_disk(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&, vector<int>&,vector<AO>&,vector<vector<int> >&, Electronic*,Electronic*, Memory*);

double scf_diis_fock(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&, vector<int>&,vector<AO>&,vector<vector<int> >&, Electronic*,Electronic*, Memory*);
double scf_diis_dm(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&, vector<int>&,vector<AO>&,vector<vector<int> >&, Electronic*,Electronic*, Memory*);


#endif // SCF_H
