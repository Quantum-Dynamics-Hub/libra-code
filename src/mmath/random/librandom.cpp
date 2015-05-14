#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "librandom.h"

using namespace boost::python;


namespace libmmath{
namespace librandom{


void export_Random_objects(){


//  double (*expt_scale1)(double, double) = &expt_scale;
//  def("scale", expt_scale1);

  class_<Random>("Random",init<>())
//      .def("__copy__", &generic__copy__<Random>)
//      .def("__deepcopy__", &generic__deepcopy__<Random>)

      .def("uniform",&Random::uniform)
      .def("p_uniform",&Random::p_uniform)

      .def("exponential",&Random::exponential)
      .def("p_exponential",&Random::p_exponential)

      .def("normal",&Random::normal)
      .def("p_normal",&Random::p_normal)

      .def("gamma",&Random::gamma)
      .def("p_gamma",&Random::p_gamma)

      .def("beta",&Random::beta)
      .def("p_beta",&Random::p_beta)

      .def("poiss1",&Random::poiss1)
      .def("poiss2",&Random::poiss2)
      .def("p_poiss",&Random::p_poiss)


  // Poisson distribution
//  int poiss(double lambda,double t);
//  void poiss(double lambda,double maxT,double dt,vector< pair<double,int> >& out);


  ;


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygrandom){
#else
BOOST_PYTHON_MODULE(librandom){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Random_objects();

}

}// librandom
}// libmmath


