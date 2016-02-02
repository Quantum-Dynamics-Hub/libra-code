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
  \file liblinalg.h
  \brief The file describes Python export function and data types
    
*/

#ifndef LIB_LINALG_H
#define LIB_LINALG_H


#include "Units.h"

#include "CMATRIX.h"
#include "MATRIX.h"                                
#include "MATRIX3x3.h"                             
#include "QUATERNION.h"                            
#include "SMATRIX.h"                               
#include "VECTOR.h" 
                               
#include "PyCopy.h"
//#include "Utility.h"



/// libmmath namespace
namespace libmmath{

/// liblinalg namespace
namespace liblinalg{


typedef std::vector<int> intList;  ///< data type for holding the list of integers
typedef std::vector<float> floatList;  ///< data type for holding the list of floats
typedef std::vector<double> doubleList;  ///< data type for holding the list of doubles
typedef std::vector<std::complex<double> > complexList;  ///< data type for holding the list of complex values

typedef std::vector<vector<int> > intMap;  ///< data type for holding the table of integers
typedef std::vector<vector<float> > floatMap;  ///< data type for holding the table of floats
typedef std::vector<vector<double> > doubleMap;  ///< data type for holding the table of doubles
typedef std::vector<vector<std::complex<double> > > complexMap;  ///< data type for holding the table of complex values



void export_linalg_objects();


}// namespace liblinalg
}// namespace libmmath


#endif// LIB_LINALG_H
