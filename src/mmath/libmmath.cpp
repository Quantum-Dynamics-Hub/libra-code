#include <boost/python.hpp>
#include "libmmath.h"
//#include "Mathematics_objects.h"
#include "PyCopy.h"
#include "Utility.h"

using namespace boost::python;

using namespace libmmath;
using libmmath::VECTOR;
using libmmath::MATRIX;
using libmmath::MATRIX3x3;
using libmmath::QUATERNION;
using libmmath::DATA;



void    (VECTOR::*cross1)(VECTOR&,VECTOR&) = &VECTOR::cross;
//VECTOR  (*cross2)(const double, const VECTOR& ,const VECTOR&) = &cross;



double (MATRIX::*get1)(int)            = &MATRIX::get;
double (MATRIX::*get2)(int,int)        = &MATRIX::get;
void   (MATRIX::*set1)(int,double)     = &MATRIX::set;
void   (MATRIX::*set2)(int,int,double) = &MATRIX::set; 
void   (MATRIX::*Inverse1)(MATRIX*)    = &MATRIX::Inverse;
void   (MATRIX::*Inverse2)(MATRIX&)    = &MATRIX::Inverse;



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmmath){
#else
BOOST_PYTHON_MODULE(libmmath){
#endif

//----------- MATHEMATICS.h -------------------------

  class_<VECTOR>("VECTOR",init<>())
      .def(init<double,double,double>())
      .def(init<const VECTOR&>())
      .def("__copy__", &generic__copy__<VECTOR>) 
      .def("__deepcopy__", &generic__deepcopy__<VECTOR>)
      .def_readwrite("x",&VECTOR::x)
      .def_readwrite("y",&VECTOR::y)
      .def_readwrite("z",&VECTOR::z)
      .def_readwrite("is_transposed",&VECTOR::is_transposed)
      .def(self+self)
      .def(self+double())
      .def(self-self)
      .def(self*double())
      .def(double()*self)
      .def(self/double())
      .def(self*self)
//      .def(self*MATRIX())
      .def(MATRIX()*self)
      .def(self==other<VECTOR>())
      .def(self!=other<VECTOR>())

      .def("length",&VECTOR::length)
      .def("length2",&VECTOR::length2)
      .def("unit",&VECTOR::unit)
      .def("normalize",&VECTOR::normalize)
      .def("cross",cross1)
//      .def("cross",cross2)
      
      
  ;

  class_<MATRIX>("MATRIX",init<>())      
      .def("__copy__", &generic__copy__<MATRIX>)
      .def("__deepcopy__", &generic__deepcopy__<MATRIX>)

      .def("get",get1)
      .def("get",get2)
      .def("set",set1)
      .def("set",set2)

      .def_readwrite("num_of_cols",&MATRIX::num_of_cols)
      .def_readwrite("num_of_rows",&MATRIX::num_of_rows)
      .def_readwrite("num_of_elems",&MATRIX::num_of_elems)
      .def(init<int,int>())
      .def(init<const MATRIX&>())
      .def(self+self)
      .def(self-self)
      .def(self*self)
      .def(self/double())
      .def(self*double())
      .def(double()*self)

      .def("Init",&MATRIX::Init)
      .def("InitSquareMatrix",&MATRIX::InitSquareMatrix)
      .def("Init_Unit_Matrix",&MATRIX::Init_Unit_Matrix)
      .def("Load_Matrix_From_File",&MATRIX::Load_Matrix_From_File)
      .def("show_num_of_rows",&MATRIX::show_num_of_rows)
      .def("show_num_of_cols",&MATRIX::show_num_of_cols)
      .def("show_num_of_elems",&MATRIX::show_num_of_elems)
      .def("RightRotation",&MATRIX::RightRotation)
      .def("LeftRotation",&MATRIX::LeftRotation)
      .def("Ortogonalization",&MATRIX::Ortogonalization)
      .def("Inverse",Inverse1)  
      .def("Inverse",Inverse2)    

  ;

  class_<MATRIX3x3>("MATRIX3x3",init<>())      
      .def(init<const VECTOR&, const VECTOR&, const VECTOR&>())
      .def(init<const MATRIX3x3&>())
      .def("__copy__", &generic__copy__<MATRIX3x3>)
      .def("__deepcopy__", &generic__deepcopy__<MATRIX3x3>)


//      .def("get",get1)
//      .def("get",get2)
//      .def("set",set1)
//      .def("set",set2)

      .def_readwrite("xx",&MATRIX3x3::xx)
      .def_readwrite("xy",&MATRIX3x3::xy)
      .def_readwrite("xz",&MATRIX3x3::xz)
      .def_readwrite("yx",&MATRIX3x3::yx)
      .def_readwrite("yy",&MATRIX3x3::yy)
      .def_readwrite("yz",&MATRIX3x3::yz)
      .def_readwrite("zx",&MATRIX3x3::zx)
      .def_readwrite("zy",&MATRIX3x3::zy)
      .def_readwrite("zz",&MATRIX3x3::zz)

      .def(self+self)
      .def(self-self)
      .def(self*self)
      .def(self/double())
      .def(self*double())
      .def(double()*self)

  ;


  class_<QUATERNION>("QUATERNION",init<>())
      .def(init<double,double,double,double>())
      .def("__copy__", &generic__copy__<QUATERNION>)
      .def("__deepcopy__", &generic__deepcopy__<QUATERNION>)

      .def_readwrite("Lt",&QUATERNION::Lt)
      .def_readwrite("Lx",&QUATERNION::Lx)
      .def_readwrite("Ly",&QUATERNION::Ly)
      .def_readwrite("Lz",&QUATERNION::Lz)
      .def(self+self)
      .def(self-self)
      .def(self*self)
      .def(self*double())
      .def(double()*self)
      .def(self+=self)
      .def(self-=self)

  ;


  class_<DATA>("DATA",init<>())
      .def(init<boost::python::list>())
      .def("Calculate_Estimators",&DATA::Calculate_Estimators)
      .def("Calculate_MiniMax",&DATA::Calculate_MiniMax)

  ;


}
