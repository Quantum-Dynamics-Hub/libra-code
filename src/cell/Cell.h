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

#ifndef CELL_H
#define CELL_H

#include "../mmath/libmmath.h"

using namespace libmmath;

namespace libcell{


struct triple{
  int is_central;
  int n1,n2,n3;
};

struct quartet{
  int is_central;
  int j;        // index of the atom with which another one is interacting
  int n1,n2,n3; // translation vector of atom j
};

struct excl_scale{
  int at_indx1, at_indx2;
  double scale;
};


class Cell{

// Auxiliary variables
  double a2,b2,c2,T12,T13,T23;
  double Xa,Xb,Xc;
  double Ya,Yb,Yc;
  double Za,Zb,Zc;
  double M,N,P,Q,R;
  double A,B,C,D,C2;
  VECTOR g1,g2,g3; // Reciprocal vectors

// Data variables
  VECTOR t1;       int is_t1;
  VECTOR t2;       int is_t2;
  VECTOR t3;       int is_t3;
  double Roff;     int is_Roff;

public:

  
  Cell(){
    is_Roff = 0;
    is_t1 = 0;
    is_t2 = 0;
    is_t3 = 0;
  }
  Cell(VECTOR& t1_,VECTOR& t2_,VECTOR& t3_, double Roff_){
    // Assign basic data
    t1 = t1_;      is_t1 = 1;
    t2 = t2_;      is_t2 = 1;
    t3 = t3_;      is_t3 = 1;
    Roff = Roff_;  is_Roff = 1;

    // Calculate auxiliary variables
    a2 = t1.length2();
    b2 = t2.length2();
    c2 = t3.length2();
    T12 = t1*t2;
    T13 = t1*t3;
    T23 = t2*t3;

    // Another round of auxiliary variables 
    // only those which depend on the cell shape
    P = (a2 - T13*T13/c2);
    M = (T13*T23/c2 - T12);
    Za = c2;
    Xa = (b2 - T23*T23/c2);
    Ya = (P*Xa - M*M);

    VECTOR g;
    g.cross(t2,t3);   g1 = g/(t1*g);
    g.cross(t3,t1);   g2 = g/(t2*g);
    g.cross(t1,t2);   g3 = g/(t3*g);
  }
  void init(VECTOR& t1_,VECTOR& t2_,VECTOR& t3_, double Roff_){
    // Assign basic data
    t1 = t1_;      is_t1 = 1;
    t2 = t2_;      is_t2 = 1;
    t3 = t3_;      is_t3 = 1;
    Roff = Roff_;  is_Roff = 1;

    // Calculate auxiliary variables
    a2 = t1.length2();
    b2 = t2.length2();
    c2 = t3.length2();
    T12 = t1*t2;
    T13 = t1*t3;
    T23 = t2*t3;

    // Another round of auxiliary variables 
    // only those which depend on the cell shape
    P = (a2 - T13*T13/c2);
    M = (T13*T23/c2 - T12);
    Za = c2;
    Xa = (b2 - T23*T23/c2);
    Ya = (P*Xa - M*M);

    VECTOR g;
    g.cross(t2,t3);   g1 = g/(t1*g);
    g.cross(t3,t1);   g2 = g/(t2*g);
    g.cross(t1,t2);   g3 = g/(t3*g);
  }

  void brute_force(VECTOR& rij, int degree, vector<triple>& res,triple& central_translation);
  void calculate(VECTOR&,vector<triple>&,triple& central_translation);
  void update_vlist(int sz,VECTOR* r, vector< vector<quartet> >& at_neib, vector<triple>& central_translation);


};


}// namespace libcell

#endif // CELL_H
