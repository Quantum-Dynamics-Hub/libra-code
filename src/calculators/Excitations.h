#ifndef EXCITATIONS_H
#define EXCITATIONS_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libcalculators{

void excite(int I, int J, vector< pair<int,double> >& occ_ini, vector< pair<int,double> >& occ_fin);


}// libcalculators

#endif // EXCITATIONS_H
