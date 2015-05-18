#ifndef EXCITATIONS_H
#define EXCITATIONS_H

#include "MOAO.h"
#include "Control_Parameters.h"
#include "Model_Parameters.h"
#include "Electronic.h"
#include "Nuclear.h"
#include "Memory.h"

#include "Engine.h"
#include "Hamiltonian.h"
#include "units.h"
#include "classes.h"


//class excitation; - see "classes.h" file

void excite(int, excitation&, int, vector< pair<int,double> >&, int, vector< pair<int,double> >&);
void compute_excitations(Control_Parameters&, Model_Parameters&, Nuclear&,
                         vector<int>&, vector<int>&, vector<AO>&,vector<vector<int> >&, Electronic*,Electronic*, Memory*);






#endif // EXCITATIONS_H
