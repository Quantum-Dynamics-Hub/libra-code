/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
* Some code is adapted from the ErgoSCF code (see the copyright info below)
*
*********************************************************************************/
/**
  \file Timer.h
  \brief The file describes and implements the Timer class for benchmarking performance
    
*/

/* Ergo, version 3.3, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2013 Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn?Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */



#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <string.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>


#include <stdexcept>
//#include "output.h"
#include "../realtype.h"

/// liblibra namespace
namespace liblibra{

using namespace libergoscf;
using namespace std;

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




/** Time-measuring class. Measures the time between the
    construction of the object and the call of the print method. */
class TimeMeter {

public:

  double startTimeCPU_sys;
  double startTimeCPU_usr;
  double startTimeWall;

  double endTimeCPU_sys;
  double endTimeCPU_usr;
  double endTimeWall;

  double secondsTakenWall;
  double secondsTakenCPU_usr;
  double secondsTakenCPU_sys;


  double get_start_time_wall_seconds() const {
    return startTimeWall;
  }

  static double get_wall_seconds() {
    struct timeval tv;
    if(gettimeofday(&tv, NULL) != 0)
      throw std::runtime_error("Error in get_wall_seconds(), in gettimeofday().");
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
  }

  static void get_current_cpu_times(double & seconds_usr, double & seconds_sys) {
    struct rusage usage;
    if(getrusage (RUSAGE_SELF, &usage) != 0)
	throw std::runtime_error("Error in get_current_cpu_times(), in getrusage().");
    seconds_usr = usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec / 1000000;
    seconds_sys = usage.ru_stime.tv_sec + (double)usage.ru_stime.tv_usec / 1000000;
  }

  TimeMeter() {
    startTimeWall = get_wall_seconds();
    get_current_cpu_times(startTimeCPU_usr, startTimeCPU_sys);

    endTimeCPU_sys = startTimeCPU_sys;
    endTimeCPU_usr = startTimeCPU_usr;
    endTimeWall = startTimeWall;

    secondsTakenCPU_sys = 0.0;
    secondsTakenCPU_usr = 0.0;
    secondsTakenWall = 0.0; 
  }

  void print(int area, const char *routine) {
    endTimeWall = get_wall_seconds();
    secondsTakenWall = endTimeWall - startTimeWall;

    get_current_cpu_times(endTimeCPU_usr, endTimeCPU_sys);  
    secondsTakenCPU_usr = endTimeCPU_usr - startTimeCPU_usr;
    secondsTakenCPU_sys = endTimeCPU_sys - startTimeCPU_sys;

//    do_output(LOG_CAT_TIMINGS, area, "%s took %9.2f usr cpu s  %9.2f sys cpu s  %9.2f wall s", 
//		routine, secondsTakenCPU_usr, secondsTakenCPU_sys, secondsTakenWall);

    /// AVA: For now
    cout<<*routine<<" took "<<secondsTakenCPU_usr<<" usr cpu s "<<secondsTakenCPU_sys<<" sys cpu s "
        <<secondsTakenWall<<" wall s\n";
  }


  
}; // TimeMeter class



}// namespace liblibra

#endif // TIMER_H
