#ifndef LIBMMATH_H
#define LIBMMATH_H

#include "Units.h"

#include "CMATRIX.h"
#include "DATA.h"                                  
#include "GRAPH.h"                                 
#include "MATRIX.h"                                
#include "MATRIX3x3.h"                             
#include "QUATERNION.h"                            
#include "random.h"                                
#include "SMATRIX.h"                               
#include "SpecialFunctions.h"                      
#include "VECTOR.h" 
                               
#include "Timer.h"
#include "PyCopy.h"
#include "Utility.h"


typedef std::vector<int> intList;
typedef std::vector<float> floatList;
typedef std::vector<double> doubleList;
typedef std::vector<std::complex<double> > complexList;


void export_Mathematics_objects();



#endif // LIBMMATH_H

