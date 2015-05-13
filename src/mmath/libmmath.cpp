#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libmmath.h"

using namespace boost::python;
//using namespace libmmath;
//using libmmath::VECTOR;
//using libmmath::MATRIX;
//using libmmath::CMATRIX;
//using libmmath::MATRIX3x3;
//using libmmath::QUATERNION;
//using libmmath::DATA;

namespace libmmath{

using namespace libspecialfunctions;
using namespace liblinalg;
using namespace libgraph;
using namespace liboperators;

void export_Mathematics_objects(){


int (DATA::*ScaleData1)(double)                = &DATA::ScaleData;
int (DATA::*ScaleData2)(double,double)         = &DATA::ScaleData;


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


  ;

  class_< DATAList >("DATAList")
      .def(vector_indexing_suite< DATAList >())
  ;


  class_<Timer>("Timer",init<>())
      .def("__copy__", &generic__copy__<Timer>)
      .def("__deepcopy__", &generic__deepcopy__<Timer>)

      .def("start", &Timer::start)
      .def("stop", &Timer::stop)
      .def("show", &Timer::show)

  ;




}// export_Mathematics_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmmath){
#else
BOOST_PYTHON_MODULE(libmmath){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
  export_linalg_objects();
  export_SpecialFunctions_objects();
  export_GRAPH_objects();
  export_Operators_objects();
  export_Mathematics_objects();

}


}// libmmath

