/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Cell.h
  \brief The file describes classes and functions for periodic calculations
    
*/

#ifndef CELL_H
#define CELL_H


#include "../math_linalg/liblinalg.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libcell namespace
namespace libcell{


struct triple{
/**
  \brief The structure of 3 integers and a flag
    
*/

  int is_central;  ///< the flag telling if this replica of the unit cell is central
  int n1,n2,n3;    ///< integers describing the position of the given cell in the 3D grid of replicas of the 
                   ///< original unit cell
};

struct quartet{
/**
  \brief The structure of 3 integers and a flag, and atom index
    
*/

  int is_central;  ///< the flag telling if this replica of the unit cell is central
  int j;           ///< index of the atom with which another one is interacting
  int n1,n2,n3;    ///< translation vector (in units of the cell dimensions) of atom j
};

struct excl_scale{
/**
  \brief The structure for exclusion factors for 2-body interactions
    
*/

  int at_indx1, at_indx2; ///< indices of the interacting atoms (this defines the excluded(or not) pair)
  double scale;           ///< the corresponding scaling factor for the interaction
};


class Cell{
/**
  \brief This class describes the properties of the unit cell
    
*/


// Auxiliary variables
  double a2,b2,c2,T12,T13,T23;
  double Xa,Xb,Xc;
  double Ya,Yb,Yc;
  double Za,Zb,Zc;
  double M,N,P,Q,R;
  double A,B,C,D,C2;
  VECTOR g1,g2,g3; ///< Reciprocal vectors

// Data variables
  VECTOR t1;       int is_t1;   ///< basis vector 1 (a direction) and the status flag
  VECTOR t2;       int is_t2;   ///< basis vector 2 (a direction) and the status flag
  VECTOR t3;       int is_t3;   ///< basis vector 3 (a direction) and the status flag
  double Roff;     int is_Roff; ///< neighbour distance cut off and the status flag

public:

  
  Cell(){  /// Default constructor
    is_Roff = 0;
    is_t1 = 0;
    is_t2 = 0;
    is_t3 = 0;
  }
  Cell(VECTOR& t1_,VECTOR& t2_,VECTOR& t3_, double Roff_){
    /// Constructor with initialization

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
    /// Initialization

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
}// liblibra

#endif // CELL_H
