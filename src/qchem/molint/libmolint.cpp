#define BOOST_PYTHON_MAX_ARITY 20
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libmolint.h"

using namespace boost::python;
using namespace libmmath;


namespace libqchem{
namespace libmolint{




void export_molint_objects(){


  // Basic overlaps
  // 1D Gaussians
  double (*expt_gaussian_overlap_ref)(int,double,double, int,double,double) = &gaussian_overlap_ref;
  double (*expt_gaussian_overlap_1D_v1)(int,double,double, int,double,double) = &gaussian_overlap;
  double (*expt_gaussian_overlap_1D_v2)(int,double,double, int,double,double, int) = &gaussian_overlap;
  boost::python::list (*expt_gaussian_overlap_1D_v3)(int,double,double, int,double,double, int, int) = &gaussian_overlap;


  // 3D Gaussians
  double (*expt_gaussian_overlap_3D_v1)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&  ) = &gaussian_overlap;
  double (*expt_gaussian_overlap_3D_v2)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&, int ) = &gaussian_overlap;
  boost::python::list (*expt_gaussian_overlap_3D_v3)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&, int, int ) = &gaussian_overlap;


  // Moments 
  // 1D Moments
  double (*expt_gaussian_moment_ref)(int,double,double, int,double, double, int,double, double ) = &gaussian_moment;
  double (*expt_gaussian_moment_1D_v1)(int,double,double, int,double,double, int,double,double ) = &gaussian_moment;
  double (*expt_gaussian_moment_1D_v2)(int,double,double, int,double,double, int,double,double, int ) = &gaussian_moment;
  boost::python::list (*expt_gaussian_moment_1D_v3)(int,double,double, int,double,double, int,double,double, int, int ) = &gaussian_moment;
 

  // 3D Moments
  double (*expt_gaussian_moment_3D_v1)(int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&          ) = &gaussian_moment;
  double (*expt_gaussian_moment_3D_v2)(int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&, int     ) = &gaussian_moment;
  boost::python::list (*expt_gaussian_moment_3D_v3)(int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&,
                                       int,int,int,double, const VECTOR&, int, int) = &gaussian_moment;



  // Pseudopotentials
  // 3D pseudopotentials
  double (*expt_pseudopot02_v1)(double,double,double,const VECTOR&,
                                int,int,int,double,const VECTOR&,
                                int,int,int,double,const VECTOR&                 ) = &pseudopot02;
  boost::python::list (*expt_pseudopot02_v2)(double, double, double, const VECTOR&,
                                int,int,int,double, const VECTOR&,
                                int,int,int,double, const VECTOR&,  int, int      ) = &pseudopot02;


  // Kinetic energy integrals (KEI)
  // 1D KEI
  double (*expt_kinetic_integral_1D_v1)(int,double, double, int, double, double)      = &kinetic_integral;
  double (*expt_kinetic_integral_1D_v2)(int,double, double, int, double, double, int) = &kinetic_integral;
  boost::python::list (*expt_kinetic_integral_1D_v3)(int,double, double, int, double, double, int, int) = &kinetic_integral;


  // 3D Gaussians
  double (*expt_kinetic_integral_3D_v1)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&          ) = &kinetic_integral;
  double (*expt_kinetic_integral_3D_v2)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&, int     ) = &kinetic_integral;
  boost::python::list (*expt_kinetic_integral_3D_v3)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&, int, int) = &kinetic_integral;


  // Nuclear attraction integrals (3D)
  double (*expt_nuclear_attraction_integral_v1)(int,int, int, double, VECTOR&,
                        int,int,int,double,VECTOR&,VECTOR&                  ) = nuclear_attraction_integral;
  double (*expt_nuclear_attraction_integral_v2)(int,int, int, double, VECTOR&,
                        int,int,int,double,VECTOR&,VECTOR&, int             ) = nuclear_attraction_integral;
  boost::python::list (*expt_nuclear_attraction_integral_v3)(int,int, int, double, VECTOR&,
                        int,int,int,double,VECTOR&,VECTOR&, int, int        ) = nuclear_attraction_integral;





  // ============ Now export functions =============
  // ==== Overlaps ====
  def("gaussian_overlap_ref", expt_gaussian_overlap_ref);
  def("gaussian_overlap", expt_gaussian_overlap_1D_v1);
  def("gaussian_overlap", expt_gaussian_overlap_1D_v2);
  def("gaussian_overlap", expt_gaussian_overlap_1D_v3);

  def("gaussian_overlap", expt_gaussian_overlap_3D_v1);
  def("gaussian_overlap", expt_gaussian_overlap_3D_v2);
  def("gaussian_overlap", expt_gaussian_overlap_3D_v3);

  def("A_coefficient_general", A_coefficient_general); 
  def("generate_coefficients", generate_coefficients);

  def("sto_norm", sto_norm);
  def("sto_overlap", sto_overlap);
  def("sto_overlap_fast", sto_overlap_fast);


  // ==== Moments ====                                     
  def("gaussian_moment_ref", expt_gaussian_moment_ref);
  def("gaussian_moment", expt_gaussian_moment_1D_v1);
  def("gaussian_moment", expt_gaussian_moment_1D_v2);
  def("gaussian_moment", expt_gaussian_moment_1D_v3);

  def("gaussian_moment", expt_gaussian_moment_3D_v1);
  def("gaussian_moment", expt_gaussian_moment_3D_v2);
  def("gaussian_moment", expt_gaussian_moment_3D_v3);


  // ==== Pseudopotentials ====
  def("pseudopot02", expt_pseudopot02_v1);
  def("pseudopot02", expt_pseudopot02_v2);


  // ==== Kinetic energy integrals ====
  def("kinetic_integral", expt_kinetic_integral_1D_v1);
  def("kinetic_integral", expt_kinetic_integral_1D_v2);
  def("kinetic_integral", expt_kinetic_integral_1D_v3);

  def("kinetic_integral", expt_kinetic_integral_3D_v1);
  def("kinetic_integral", expt_kinetic_integral_3D_v2);
  def("kinetic_integral", expt_kinetic_integral_3D_v3);


  // ==== Nuclear attraction integral ====
  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_v1);
  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_v2);
  def("nuclear_attraction_integral", expt_nuclear_attraction_integral_v3);



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmolint){
#else
BOOST_PYTHON_MODULE(libmolint){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_molint_objects();

}


}// namespace libmolint
}// namespace libqchem

