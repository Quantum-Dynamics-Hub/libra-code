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
#include "Utility.h"




namespace libmmath{
namespace liblinalg{


typedef std::vector<int> intList;
typedef std::vector<float> floatList;
typedef std::vector<double> doubleList;
typedef std::vector<std::complex<double> > complexList;


void export_linalg_objects();


}// namespace liblinalg
}// namespace libmmath


#endif// LIB_LINALG_H
