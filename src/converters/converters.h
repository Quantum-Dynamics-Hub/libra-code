#ifndef CONVERTERS_H
#define CONVERTERS_H

#include "../chemobjects/libchemobjects.h"
#include "../dyn/libdyn.h"

using namespace libchemobjects;
using namespace libdyn;


namespace libconverters{


void system_to_nuclear(System& syst, Nuclear& nucl);
void nuclear_to_system(Nuclear& nucl, System& syst);

void system_to_vector_q(System& syst, vector<double>& q);
void system_to_vector_p(System& syst, vector<double>& p);




}// namespace libconverters

#endif // CONVERTERS_H

