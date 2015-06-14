#ifndef MULLIKEN_H
#define MULLIKEN_H

#include "../mmath/libmmath.h"
using namespace libmmath;


namespace libcalculators{


void update_Mull_orb_pop(MATRIX*, MATRIX*, vector<double>&, vector<double>&);

void update_Mull_charges(vector<int>&, vector<int>&, vector<vector<int> >&,vector<double>&,
                         vector<double>&, vector<double>&, vector<double>&, vector<double>&);


}// namespace libcalculators

#endif // MULLIKEN_H

