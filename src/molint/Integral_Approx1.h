/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef INTEGRAL_APPROX1_H
#define INTEGRAL_APPROX1_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#endif 

/// liblibra namespace
namespace liblibra{

using namespace std;
 
namespace libmolint{



//============================================================
//=== auxiliary functions according to Rosen =================
double A_plus_1(int n,double alpha);
double A_minus_1(int n,double alpha);
double B(int n,double alpha);
double D(int m,int n,int p);
double K2ab(int m,int n,double alpha,double betha,double R);
double K2aa(int m,int n,double alpha,double betha,double R);

//============================================================
//=== According to O. Kitao and T. Ogawa =====================
double Jab(int m,int n,double alpha,double betha,double R);
double Jab(int m,int n,double alpha,double betha);
double dJab_dqa(int m,int n,double alpha,double betha,double R);

//=============================================================
//==== Parameterized way to calculate Coulomb integral ========
double Coulomb_Integral(double R,int n_i, double Jii, double ksi_i, std::string type_i, double q_i,
                                 int n_j, double Jjj, double ksi_j, std::string type_j, double q_j,
                                 double epsilon, int mode);


}// namespace libmolint
}// namespace liblibra


#endif // INTEGRAL_APPROX1_H

