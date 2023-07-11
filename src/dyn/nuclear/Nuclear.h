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
  \file Nuclear.h
  \brief The file describes Nuclear class and functions for propagation of nuclear DOF
    
*/

#ifndef NUCLEAR_H
#define NUCLEAR_H

#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libdyn namespace
namespace libdyn{

/// libnuclear namespace
namespace libnuclear{


class Nuclear{
/** 
  \brief This clas represents a point in a N-dimensional phase space, N = 6*Nnucl
  plus some properties are associated with each projection (coordinate)

*/
  public:

  //------------- Data members ---------------------
  int nnucl;                      ///< the number of nuclear DOFs

  // Atomic properties
  vector<double> mass;            ///< (generalized) masses of particles

  // Dynamic variables
  vector<double> q;               ///< nuclear DOFs - these are generalized scalar coordinates 
  vector<double> p;               ///< conjugate momenta
  vector<double> f;               ///< forces
  vector<int> ctyp;               ///< coordinate type, for instance one can set 0,1,2 to be x,y,z projections of all atoms, so
                                  ///< in this case ctyp will be: [0, 1, 2, 0, 1, 2, 0, 1, 2. ..]
                                  ///< or one can set 0 to be the x coordinates of atoms in given fragment, so the ordering may look like:
                                  ///< [0, 0, 0, 0, 0, ... 1, 1, 1, 1, 1, ... ]



  //------------- Constructors ----------------------
  Nuclear();
  Nuclear(int _nnucl);
  Nuclear(const Nuclear&);
 ~Nuclear();

  void propagate_p(int i,double dt);
  void propagate_p(double dt);
  void propagate_p(double dt,vector<int>& active);

  void scale_p(int i,double scl);
  void scale_p(double scl);
  void scale_p(double scl,vector<int>& active);

  void propagate_q(int i,double dt);
  void propagate_q(double dt);
  void propagate_q(double dt,vector<int>& active);

  void scale_q(int i,double scl);
  void scale_q(double scl);
  void scale_q(double scl,vector<int>& active);

  friend bool operator == (const Nuclear& n1, const Nuclear& n2){
    return &n1 == &n2;
  }
  friend bool operator != (const Nuclear& n1, const Nuclear& n2){
    return !(n1 == n2);  // only compare addresses
  }



  
};

typedef std::vector< Nuclear > NuclearList; ///< Type containing the vector of Nuclear objects



} // libnuclear
} // libdyn
}// liblibra

#endif // NUCLEAR_H
