#ifndef ELECTRONIC_H
#define ELECTRONIC_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libdyn{
namespace libelectronic{

class Electronic{
  // This class describes a time-dependent state (multiconfigurational in sense that 
  // it can have time-dependent contributions from different stationary states )
  

  void rnd_phase(double&, double&, double);

  public:

  int nstates;             // number of stationary (basis) states
  int istate;              // index of current basis state (for stochastic)
  std::vector<double> q;   // MMTS variables for all basis states Re(c)
  std::vector<double> p;   // MMTS variables        --- // ---    Im(c)


  //--------------- Class functions ---------------

  void init(int,int);
  void init(int);

  Electronic();
  Electronic(int);
  Electronic(int,int);

};

}// namespace libelectronic
}// namespace libdyn

#endif // ELECTRONIC_H
