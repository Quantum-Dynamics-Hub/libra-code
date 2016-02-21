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
 \file DIIS.h
 \brief The file describes the DIIS class and related functions
        
*/


#ifndef DIIS_H
#define DIIS_H

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

#include "../mmath/libmmath.h"
using namespace libmmath;


/// libsolvers namespace
namespace libsolvers{


class DIIS{
/**
  This is the class that handles DIIS (direct inversion of the iterative space) method
*/

public:

  DIIS(int _N_diis_max,int Norb);  ///< Constructor

  void update_diis_coefficients();
  void add_diis_matrices(MATRIX* X, MATRIX* err);
  void add_diis_matrices(MATRIX& X, MATRIX& err);

  void extrapolate_matrix(MATRIX* X_ext);
  void extrapolate_matrix(MATRIX& X_ext);



  int N_diis_max;                ///< Length of DIIS history (size of the lists)  

  int N_diis;                    ///< current # of matrices stored
  int N_diis_eff;                ///< effective size of the DIIS matrix (such that it is full-rank)
  vector<MATRIX*> diis_X;        ///< diis iteration of objective matrices (typically Fock matrices)
  vector<MATRIX*> diis_err;      ///< diis error matrices
  vector<double>  diis_c;        ///< diis extrapolation coefficients


  boost::python::list get_diis_X();
  boost::python::list get_diis_err();
  boost::python::list get_diis_c();

};

}// libsolvers namespace

#endif // DIIS_H 
