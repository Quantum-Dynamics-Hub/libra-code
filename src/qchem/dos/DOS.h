#ifndef DENSITY_OF_STATES_H
#define DENSITY_OF_STATES_H

#include "MOAO.h"
#include "Control_Parameters.h"
#include "Model_Parameters.h"
#include "Electronic.h"
#include "Nuclear.h"
#include "Memory.h"


void compute_dos(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                 vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                 Electronic* el, Memory* mem);


#endif // DENSITY_OF_STATES_H
