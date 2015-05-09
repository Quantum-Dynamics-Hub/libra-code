#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libpot.h"

using namespace libpot;
using namespace boost::python;


void export_Pot_objects(){


double (*Vdw_LJ_1)(VECTOR&, VECTOR&, VECTOR&, VECTOR&, double, double) =  &Vdw_LJ;

double (*Vdw_LJ_2)(VECTOR* r,VECTOR* g,VECTOR* m,VECTOR* f,MATRIX3x3& at_stress, 
                  MATRIX3x3& fr_stress, MATRIX3x3& ml_stress,
                  int sz,double* epsilon, double* sigma,
                  int nexcl, int* excl1, int* excl2, double* scale,
                  MATRIX3x3* box,int rec_deg,int pbc_deg,
                  double etha,int is_cutoff, double R_on, double R_off,
                  int& time, vector< vector<triple> >& images, vector<triple>& central_translation,
                  double* dr2, double dT, int& is_update
                 ) = &Vdw_LJ;



  def("SWITCH", SWITCH);
  def("DOUBLE_SWITCH", DOUBLE_SWITCH);
  def("Bond_Harmonic", Bond_Harmonic);
  def("Bond_Quartic", Bond_Quartic);
  def("Bond_Morse", Bond_Morse);

  def("Angle_Harmonic", Angle_Harmonic);
  def("Angle_Fourier", Angle_Fourier);
  def("Angle_Fourier_General", Angle_Fourier_General);
  def("Angle_Fourier_Special", Angle_Fourier_Special);
  def("Angle_Harmonic_Cos", Angle_Harmonic_Cos);
  def("Angle_Harmonic_Cos_General", Angle_Harmonic_Cos_General);
  def("Angle_Cubic", Angle_Cubic);

  def("Stretch_Bend_Harmonic", Stretch_Bend_Harmonic);
  def("Dihedral_General", Dihedral_General);
  def("Dihedral_Fourier", Dihedral_Fourier);

  def("OOP_Fourier", OOP_Fourier);
  def("OOP_Wilson", OOP_Wilson);
  def("OOP_Harmonic", OOP_Harmonic);


  def("Vdw_LJ", Vdw_LJ_1);

  def("Vdw_Buffered14_7", Vdw_Buffered14_7);
  def("Vdw_Morse", Vdw_Morse);
  def("Elec_Coulomb", Elec_Coulomb);
//  def("Gay_Berne", Gay_Berne);
  def("Girifalco12_6", Girifalco12_6);

//  def("Elec_Ewald3D", Elec_Ewald3D);
//  def("Vdw_LJ", Vdw_LJ_2);
//  def("Vdw_LJ1", Vdw_LJ1);
//  def("Vdw_LJ2_no_excl", Vdw_LJ2_no_excl);
//  def("Vdw_LJ2_excl", Vdw_LJ2_excl);
//  def("LJ_Coulomb", LJ_Coulomb);
//  def("", );


} // export_Pot_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygpot){
#else
BOOST_PYTHON_MODULE(libpot){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Pot_objects();

}


