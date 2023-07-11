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
/**
  \file liblinalg.h
  \brief The file describes Python export function and data types
    
*/

#ifndef LIB_LINALG_H
#define LIB_LINALG_H

#include "permutations.h"
#include "base_matrix.h"  
#include "IMATRIX.h"                               
#include "MATRIX.h"                               
#include "CMATRIX.h"
#include "MATRIX3x3.h" 
#include "QUATERNION.h"  
#include "VECTOR.h"
#include "FT.h"
#include "Mathematics.h"
#include "PyCopy.h"


/// liblibra namespace
namespace liblibra{


/// liblinalg namespace
namespace liblinalg{


typedef std::vector<int> intList;  ///< data type for holding the list of integers
typedef std::vector<float> floatList;  ///< data type for holding the list of floats
typedef std::vector<double> doubleList;  ///< data type for holding the list of doubles
typedef std::vector<std::complex<double> > complexList;  ///< data type for holding the list of complex values

typedef std::vector< std::vector<int> > intList2;       ///< data type for holding the list of lists of integers
typedef std::vector< std::vector<float> > floatList2;   ///< data type for holding the list of lists of floats
typedef std::vector< std::vector<double> > doubleList2; ///< data type for holding the list of lists of doubles
typedef std::vector< std::vector<std::complex<double> > > complexList2;  ///< data type for holding the list of lists of complex values

typedef std::vector< std::vector<  std::vector<int> > > intList3;    ///< data type for holding the list of lists of lists of integers
typedef std::vector< std::vector<  std::vector<float> > > floatList3;    ///< data type for holding the list of lists of lists of floats
typedef std::vector< std::vector<  std::vector<double> > > doubleList3;    ///< data type for holding the list of lists of lists of double
typedef std::vector< std::vector<  std::vector< complex<double> > > > complexList3;    ///< data type for holding the list of lists of lists of complex


typedef std::vector<vector<int> > intMap;  ///< data type for holding the table of integers
typedef std::vector<vector<float> > floatMap;  ///< data type for holding the table of floats
typedef std::vector<vector<double> > doubleMap;  ///< data type for holding the table of doubles
typedef std::vector<vector<std::complex<double> > > complexMap;  ///< data type for holding the table of complex values





void export_linalg_objects();


}// namespace liblinalg
}// namespace liblibra


#endif// LIB_LINALG_H
