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

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libdata.h"
#include "../linalg/PyCopy.h"

using namespace boost::python;


namespace libmmath{
namespace libdata{


void export_Data_objects(){


int (DATA::*ScaleData1)(double)                = &DATA::ScaleData;
int (DATA::*ScaleData2)(double,double)         = &DATA::ScaleData;

boost::python::list (DATA::*expt_Calculate_Distribution)(boost::python::list Interval) = &DATA::Calculate_Distribution;


  class_<DATA>("DATA",init<>())
      .def(init<boost::python::list>())
      .def("__copy__", &generic__copy__<DATA>)
      .def("__deepcopy__", &generic__deepcopy__<DATA>)

//      .def("Calculate_Estimators",&DATA::Calculate_Estimators)
//      .def("Calculate_MiniMax",&DATA::Calculate_MiniMax)

      .def("LinearTransformData", &DATA::LinearTransformData)
      .def("ScaleData", ScaleData1)
      .def("ScaleData", ScaleData2)
      .def("ShiftData", &DATA::ShiftData)
      .def("NormalizeData",  &DATA::NormalizeData)

      .def_readwrite("Data",&DATA::Data)

      .def_readwrite("ave",&DATA::ave)
      .def_readwrite("var",&DATA::var)
      .def_readwrite("sd",&DATA::sd)
      .def_readwrite("se",&DATA::se)
      .def_readwrite("mse",&DATA::mse)
      .def_readwrite("mae",&DATA::mae)
      .def_readwrite("rmse",&DATA::rmse)

      .def_readwrite("min",&DATA::min)
      .def_readwrite("min_indx",&DATA::min_indx)
      .def_readwrite("max",&DATA::max)
      .def_readwrite("max_indx",&DATA::max_indx)

      .def_readwrite("scale_factor",&DATA::scale_factor)
      .def_readwrite("shift_amount",&DATA::shift_amount)


      .def("Calculate_Distribution",expt_Calculate_Distribution)
       
  ;

  class_< DATAList >("DATAList")
      .def(vector_indexing_suite< DATAList >())
  ;



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygdata){
#else
BOOST_PYTHON_MODULE(libdata){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Data_objects();

}

}// libdata
}// libmmath


