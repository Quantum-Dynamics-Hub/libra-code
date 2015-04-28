/*********************************************************************************
* Copyright (C) 2012 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef RANDOM_H
#define RANDOM_H

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;


// Uniform distribution
double uniform(double a,double b);

// Exponential distribution
double exponential(double lambda);

// Normal (Gaussian) distribution
double normal();

// Gamma distribution
double gamma(double a);

// Beta distribution
double beta(double a,double b);

// Poisson distribution
int poiss(double lambda,double t);
void poiss(double lambda,double maxT,double dt,vector< pair<double,int> >& out);
int poiss1(double lambda);
int poiss2(double lambda);


#endif // RANDOM_H
