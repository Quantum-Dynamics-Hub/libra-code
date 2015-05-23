#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libspecialfunctions.h"

using namespace boost::python;

//using namespace libmmath;
//using namespace libmmath::libspecialfunctions;
//using namespace libmmath::liblinalg;
//using libmmath::VECTOR;
//using libmmath::MATRIX;
//using libmmath::CMATRIX;
//using libmmath::MATRIX3x3;
//using libmmath::QUATERNION;
//using libmmath::DATA;



namespace libmmath{
namespace libspecialfunctions{
//namespace liblinalg{

void export_SpecialFunctions_objects(){

  boost::python::list (*expt_binomial_expansion)(int, int, double, double, int) = &binomial_expansion;

  // Now introduce normal functions:
  def("FAST_POW", FAST_POW);

  def("sinh_",sinh_);  // sinh(x)/x
  def("sin_", sin_);   // sin(x)/x
  def("ERF",ERF);      // error function
  def("ERFC",ERFC);    // complementary error function
  def("gamma_lower", gamma_lower);  // lower gamma function divided by the power
  def("Fn", Fn);
  def("gaussian_int", gaussian_int);  
  def("gaussian_norm", gaussian_norm);

  def("FACTORIAL",FACTORIAL); // n!
  def("DRACTORIAL",DFACTORIAL);  // n!!
  def("BINOM",BINOM);  // C_n^i
  def("binomial_expansion", expt_binomial_expansion);

  def("LEGENDRE",LEGENDRE);
  def("CHEBYSHEV",CHEBYSHEV);
  def("LAGUERRE",LAGUERRE);
  def("HERMITE",HERMITE); 

  def("Ellipe",Ellipe);
  def("Ellipe2",Ellipe2);
  def("Jacobi_Elliptic", Jacobi_Elliptic);
  def("Km", Km);
  def("Ellint",Ellint);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygspecialfunctions){
#else
BOOST_PYTHON_MODULE(libspecialfunctions){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
  export_SpecialFunctions_objects();

}

}// libspecialfunctions
}// libmmath
