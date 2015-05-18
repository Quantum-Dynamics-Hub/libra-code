#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libmolint.h"

using namespace boost::python;
using namespace libmmath;


namespace libqchem{
namespace libmolint{


void export_molint_objects(){


  // Basic functions
  double (*expt_gaussian_int)(int, double) = &gaussian_int;
  double (*expt_gaussian_norm)(int,double) = &gaussian_norm;

  // Basic overlaps
  // 1D Gaussians
  double (*expt_gaussian_overlap_1D)(int,double, double, int, double, double) = &gaussian_overlap;

  // 3D Gaussians
  double (*expt_gaussian_overlap_3D_v1)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&          ) = &gaussian_overlap;

  double (*expt_gaussian_overlap_3D_v2)(int,int,int,double,VECTOR&,
                                        int,int,int,double,VECTOR&,
                                        int, VECTOR&, VECTOR&               ) = &gaussian_overlap;



  // Moments - Essentially the generalized overlaps
  // 1D version
  double (*expt_gaussian_moment_1D)(int,double,double, int,double, double, int,double, double ) = &gaussian_moment;

  // 3D version
  double (*expt_gaussian_moment_3D_v1)(int,int,int,double,VECTOR&,
                                       int,int,int,double,VECTOR&,
                                       int,int,int,double,VECTOR&         ) = &gaussian_moment;

/*
// Don't know why I can't export the function below, but it is not really needed anyways, so
// just comment for now
  double (*expt_gaussian_moment_3D_v2)(int nx, int ny,  int nz,  double alp, VECTOR& R,
                       int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                       int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                       int is_normalize, 
                       VECTOR& dIdA, VECTOR& dIdB
                      ) = &gaussian_moment;
*/


  // Pseudopotentials - essentially a combination of moments and overlaps

  double (*expt_pseudopot02_v1)(double,double,double,VECTOR&,
                                int,int,int,double,VECTOR&,
                                int,int,int,double,VECTOR&                 ) = &pseudopot02;

/*
  double (*expt_pseudopot02_v2)(double C0, double C2, double alp, VECTOR& R,
                   int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                   int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                   int is_normalize, 
                   VECTOR& dIdA, VECTOR& dIdB
                  ) = &pseudopot02;
*/


  // ============ Now export functions =============

  def("gaussian_int", expt_gaussian_int);
  def("gaussian_norm", expt_gaussian_norm);
  def("gaussian_overlap", expt_gaussian_overlap_1D);
  def("gaussian_overlap", expt_gaussian_overlap_3D_v1);
  def("gaussian_overlap", expt_gaussian_overlap_3D_v2);
  def("gaussian_moment", expt_gaussian_moment_1D);
  def("gaussian_moment", expt_gaussian_moment_3D_v1);
//  def("gaussian_moment", expt_gaussian_moment_3D_v2);
  def("pseudopot02", expt_pseudopot02_v1);
//  def("pseudopot02", expt_pseudopot02_v2);



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

