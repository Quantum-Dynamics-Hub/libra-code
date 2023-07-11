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

#ifndef QEQ_H
#define QEQ_H

#include "Utility.h"

// -------------- Some auxiliary structures -----------
struct _GMP{

  std::string type;// atom element

  int    n;   // quantum number
  double r;   // covalent radius in Angstoms
  double ksi; // exponent (a.u.)^-1
  double xi;  // electronegativities
  double J;   // idempotential of a given atom

};
struct _MOL{

        vector<string> type;
        vector<VECTOR> R; // coordinates
        vector<double> Q; // initial carges
        vector<_GMP> GMP; // GMP electronegativities

};

class Qeq{

  //***************** Variables **************
  std::string idempotential;        int is_idempotential;   // Way to get idempotential values (parameters)
  // allowed values = rappe_goddard, parr_pearson, bakowies_thiel
  int solution_method;              int is_solution_method; // Solution method (which algorithm to use)
  int integral_method;              int is_integral_method; // How we calculate Jab integrals
  // allowed values = 0, 1, 2, 3
  double epsilon;                   int is_epsilon;         // Dielectric constant
  double threshold;                 int is_threshold;       // Convergence precision
  int MaxCount;                     int is_MaxCount;        // Maximal number of iterations

  //--------- Auxiliary internal functions ------------
  void init_variables();     // Initializes dynamical variables
  void copy_content(const Qeq&); // Copies the content which is defined

public:
  //----------- Basic class operations ---------------------------
  // Defined in Qeq.cpp
  Qeq();            // constructor
  Qeq(const Qeq&);  // copy-constructor
 ~Qeq();            // destructor

  Qeq& operator=(const Qeq&); // assignment operator
  void show_info();
  void set(object);

  //-------------- Getters, setters -----------------------------
  // Defined in Qeq_aux.cpp

  //------------------- Methods ---------------------------------
  // Defined in Qeq_methods.cpp


};

/**************************************************************
   Auxiliary functions: Defined in Qeq_aux.cpp
***************************************************************/

//============================================================
//=== According to Rosen =====================================

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
//==== Parameterized way to calculate Coulomb integral ==========
double Coulomb_Integral(double& R,int i,int j,double& qi,double& qj,vector<_GMP>& GMP,double& epsilon,int& mode);


/****************************************************************
     Method functions : Defined in Qeq_methods.cpp
****************************************************************/
void Qeq_Init(std::string,map<std::string,_GMP>& GMP);
void Solve(_MOL& mol,double epsilon,int mode);



#endif // QEQ_H
