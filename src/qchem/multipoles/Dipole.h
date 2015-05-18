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
