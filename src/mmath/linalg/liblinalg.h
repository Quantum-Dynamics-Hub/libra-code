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




namespace libmmath{
namespace liblinalg{


typedef std::vector<int> intList;
typedef std::vector<float> floatList;
typedef std::vector<double> doubleList;
typedef std::vector<std::complex<double> > complexList;

typedef std::vector<vector<int> > intMap;
typedef std::vector<vector<float> > floatMap;
typedef std::vector<vector<double> > doubleMap;
typedef std::vector<vector<std::complex<double> > > complexMap;



void export_linalg_objects();


}// namespace liblinalg
}// namespace libmmath


#endif// LIB_LINALG_H
