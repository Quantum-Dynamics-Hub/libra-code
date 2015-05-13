#include "Nuclear.h"

namespace libdyn{
namespace libnuclear{


Nuclear::Nuclear(){ n_nucl = 0; }

Nuclear::Nuclear(int _n_nucl){  
    n_nucl = _n_nucl;
    mass = vector<double>(n_nucl,2000.0); 
    q = vector<double>(n_nucl,0.0);
    p = vector<double>(n_nucl,0.0);
    f = vector<double>(n_nucl,0.0);
    
}


}// namespace libnuclear
}// namespace libdyn


