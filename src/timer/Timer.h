/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Timer.h
  \brief The file describes and implements the Timer class for benchmarking performance
    
*/

#ifndef TIMER_H
#define TIMER_H

#include <time.h>

/// liblibra namespace
namespace liblibra{


class Timer{
/**
  The Timer class which can be used for benchmarking purposes
*/

 time_t t1,t2; // start and end points
 time_t acc;   // accumulator 

public:

  Timer(){ acc = 0.0; } ///< Constructor: resets accumulated time to zero

  inline void start(){ t1 = clock(); }  ///< Start: saves the time of the start
  inline double stop(){
  /** Stop: gets the time of call and computes the time difference w.r.t to the start time. The difference is returned
  but is also added to the internal accumulator
  */
    t2 = clock(); acc += (t2-t1); return ((t2-t1)/((double)CLOCKS_PER_SEC));
  }
  inline double show(){  return (acc/((double)CLOCKS_PER_SEC)); }  ///< Returns the time accumulated so far (in between start/stop) calls


};

}// namespace liblibra

#endif // TIMER_H
