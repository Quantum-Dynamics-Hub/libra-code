#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libcalculators.h"

namespace libcalculators{


void export_calculators_objects(){


  double (*expt_fermi_population_v1)
  (double e,double ef,double degen, double kT) = &fermi_population;

  double (*expt_fermi_integral_v1)
  (vector<double>& bnds, double ef, double degen, double kT) = &fermi_integral;

  double (*expt_fermi_energy_v1)
  (vector<double>& bnds,double Nel,double degen, double kT, double etol) = &fermi_energy;

 
  def("fermi_population", expt_fermi_population_v1);
  def("fermi_integral", expt_fermi_integral_v1);
  def("fermi_energy", expt_fermi_energy_v1);


}// export_calculators_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygcalculators){
#else
BOOST_PYTHON_MODULE(libcalculators){
#endif

  export_calculators_objects();

}


}// namespace libcalculators




