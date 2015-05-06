#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libmmath.h"
//#include "Mathematics_objects.h"
#include "PyCopy.h"
#include "Utility.h"

using namespace boost::python;

using namespace libmmath;
using libmmath::VECTOR;
using libmmath::MATRIX;
using libmmath::CMATRIX;
using libmmath::MATRIX3x3;
using libmmath::QUATERNION;
using libmmath::DATA;


void export_Mathematics_objects(){


  class_< intList >("intList")
      .def(vector_indexing_suite< intList >())
  ;

  class_< floatList >("floatList")
      .def(vector_indexing_suite< floatList >())
  ;

  class_< doubleList >("doubleList")
      .def(vector_indexing_suite< doubleList >())
  ;

  class_< complexList >("complexList")
      .def(vector_indexing_suite< complexList >())
  ;



  


void    (VECTOR::*cross1)(VECTOR&,VECTOR&) = &VECTOR::cross;
//VECTOR  (*cross2)(const double, const VECTOR& ,const VECTOR&) = &cross;


double (MATRIX::*get1)(int)            = &MATRIX::get;
double (MATRIX::*get2)(int,int)        = &MATRIX::get;
void   (MATRIX::*set1)(int,double)     = &MATRIX::set;
void   (MATRIX::*set2)(int,int,double) = &MATRIX::set; 
void   (MATRIX::*Inverse1)(MATRIX*)    = &MATRIX::Inverse;
void   (MATRIX::*Inverse2)(MATRIX&)    = &MATRIX::Inverse;

complex<double> (CMATRIX::*get3)(int)            = &CMATRIX::get;
complex<double> (CMATRIX::*get4)(int,int)        = &CMATRIX::get;
void   (CMATRIX::*set3)(int,double,double)       = &CMATRIX::set;
void   (CMATRIX::*set4)(int,int,double,double)   = &CMATRIX::set; 

void (CMATRIX::*tridiagonalize1)(CMATRIX& T)              = &CMATRIX::tridiagonalize;
void (CMATRIX::*tridiagonalize2)(CMATRIX& T,CMATRIX& H)   = &CMATRIX::tridiagonalize;


int (DATA::*ScaleData1)(double)                = &DATA::ScaleData;
int (DATA::*ScaleData2)(double,double)         = &DATA::ScaleData;


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

  class_< VECTORList >("VECTORList")
      .def(vector_indexing_suite< VECTORList >())
  ;



  class_<MATRIX>("MATRIX",init<>())      
      .def(init<int,int>())
      .def(init<const MATRIX&>())
      .def("__copy__", &generic__copy__<MATRIX>)
      .def("__deepcopy__", &generic__deepcopy__<MATRIX>)

      .def("get",get1)
      .def("get",get2)
      .def("set",set1)
      .def("set",set2)

      .def_readwrite("num_of_cols",&MATRIX::num_of_cols)
      .def_readwrite("num_of_rows",&MATRIX::num_of_rows)
      .def_readwrite("num_of_elems",&MATRIX::num_of_elems)
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

  class_< MATRIXList >("MATRIXList")
      .def(vector_indexing_suite< MATRIXList >())
  ;


  class_<MATRIX3x3>("MATRIX3x3",init<>())      
      .def(init<const VECTOR&, const VECTOR&, const VECTOR&>())
      .def(init<const MATRIX3x3&>())
      .def("__copy__", &generic__copy__<MATRIX3x3>)
      .def("__deepcopy__", &generic__deepcopy__<MATRIX3x3>)


      .def("get",get1)
      .def("get",get2)
      .def("set",set1)
      .def("set",set2)

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

  class_< MATRIX3x3List >("MATRIX3x3List")
      .def(vector_indexing_suite< MATRIX3x3List >())
  ;


  class_<CMATRIX>("CMATRIX",init<>())      
      .def(init<int,int>())
      .def(init<const CMATRIX&>())
      .def("__copy__", &generic__copy__<CMATRIX>)
      .def("__deepcopy__", &generic__deepcopy__<CMATRIX>)

      .def("get",get3)
      .def("get",get4)
      .def("set",set3)
      .def("set",set4)

      .def_readwrite("num_of_cols",&CMATRIX::n_cols)
      .def_readwrite("num_of_rows",&CMATRIX::n_rows)
      .def_readwrite("num_of_elems",&CMATRIX::n_elts)
      .def(self+self)
      .def(self-self)
      .def(self*self)
      .def(self/double())
      .def(self*double())
      .def(double()*self)

      .def("show", &CMATRIX::show)
      .def("conj", &CMATRIX::conj)   // return complex conjugate matrix
      .def("T", &CMATRIX::T)         // return transposed matrix
      .def("H", &CMATRIX::H)         // return Hermitian-conjugate matrix
      .def("load_identity", &CMATRIX::load_identity)
      .def("dot", &CMATRIX::dot)     // dot product of two matrices
      .def("col", &CMATRIX::col)     // return given column of the matrix 
      .def("row", &CMATRIX::row)     // return given row of the matrix 

      .def("QR", &CMATRIX::QR)       // QR for general (Hermitian or symmetric) matrices
      .def("QR1",&CMATRIX::QR1)      // QR for tridiagonal matrices

      .def("tridiagonalize", tridiagonalize1)  // for Hermitian or symmetric CMATRIX - only resulting tridiagonal CMATRIX
      .def("tridiagonalize", tridiagonalize2)  // ---//---  also keep track of Householder transformation matrices

      .def("eigen" , &CMATRIX::eigen )   // interface
      .def("eigen0", &CMATRIX::eigen0)   // Schur decomposition or Jacobi rotation
      .def("eigen1", &CMATRIX::eigen1)   // only eigenvalues - fast
      .def("eigen2", &CMATRIX::eigen2)   // also eigenvectors, slower - keep track of transformation matrixes
      .def("eigen3", &CMATRIX::eigen3)   // also eigenvectors - solve for each eigenvector independently

      .def("inverse",        &CMATRIX::inverse)        // interface
      .def("direct_inverse", &CMATRIX::direct_inverse)



/*
      .def("Init",&CMATRIX::Init)
      .def("InitSquareMatrix",&CMATRIX::InitSquareMatrix)
      .def("Init_Unit_Matrix",&CMATRIX::Init_Unit_Matrix)
      .def("Load_Matrix_From_File",&CMATRIX::Load_Matrix_From_File)
      .def("show_num_of_rows",&CMATRIX::show_num_of_rows)
      .def("show_num_of_cols",&CMATRIX::show_num_of_cols)
      .def("show_num_of_elems",&CMATRIX::show_num_of_elems)
      .def("RightRotation",&CMATRIX::RightRotation)
      .def("LeftRotation",&CMATRIX::LeftRotation)
      .def("Ortogonalization",&CMATRIX::Ortogonalization)
//      .def("Inverse",Inverse1)  
//      .def("Inverse",Inverse2)    
*/
  ;

  class_< CMATRIXList >("CMATRIXList")
      .def(vector_indexing_suite< CMATRIXList >())
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

  class_< QUATERNIONList >("QUATERNIONList")
      .def(vector_indexing_suite< QUATERNIONList >())
  ;



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



  // Now introduce normal functions:
  def("FAST_POW", FAST_POW);

  def("sinh_",sinh_);  // sinh(x)/x
  def("sin_", sin_);   // sin(x)/x
  def("ERF",ERF);      // error function
  def("ERFC",ERFC);    // complementary error function

  def("FACTORIAL",FACTORIAL); // n!
  def("DRACTORIAL",DFACTORIAL);  // n!!
  def("BINOM",BINOM);  // C_n^i

  def("LEGENDRE",LEGENDRE);
  def("CHEBYSHEV",CHEBYSHEV);
  def("LAGUERRE",LAGUERRE);
  def("HERMITE",HERMITE); 

  def("Ellipe",Ellipe);
  def("Ellipe2",Ellipe2);
  def("Jacobi_Elliptic", Jacobi_Elliptic);
  def("Km", Km);
  def("Ellint",Ellint);


}// export_Mathematics_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmmath){
#else
BOOST_PYTHON_MODULE(libmmath){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Mathematics_objects();

}
