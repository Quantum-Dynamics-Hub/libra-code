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
  inline double stop(){ t2 = clock(); acc += (t2-t1); return ((t2-t1)/((double)CLOCKS_PER_SEC)); }
  inline double show(){  return (acc/((double)CLOCKS_PER_SEC)); }


};

}// namespace libmmath

#endif // TIMER_H
