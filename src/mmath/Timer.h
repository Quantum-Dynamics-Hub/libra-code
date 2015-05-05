#ifndef TIMER_H
#define TIMER_H

#include <time.h>

namespace libmmath{

class Timer{

 time_t t1,t2; // start and end points
 time_t acc;   // accumulator 

public:

  Timer(){ acc = 0.0; }

  inline void start(){ t1 = clock(); }
  inline void stop(){ t2 = clock(); acc += (t2-t1); }
  inline double show(){  return (acc/((double)CLOCKS_PER_SEC)); }


};

}// namespace libmmath

#endif // TIMER_H
