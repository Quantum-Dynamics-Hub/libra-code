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

