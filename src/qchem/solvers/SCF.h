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
