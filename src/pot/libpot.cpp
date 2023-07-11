/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libpot.h"

/// liblibra namespace
namespace liblibra{


using namespace boost::python;
using namespace libcell;
using namespace liblinalg;

namespace libpot{


void export_Pot_objects(){

  double (*Vdw_LJ_1)(VECTOR&, VECTOR&, VECTOR&, VECTOR&, double, double) =  &Vdw_LJ;
/*

double (*Vdw_LJ_2)(VECTOR* r,VECTOR* g,VECTOR* m,VECTOR* f,MATRIX3x3& at_stress, 
                  MATRIX3x3& fr_stress, MATRIX3x3& ml_stress,
                  int sz,double* epsilon, double* sigma,
                  int nexcl, int* excl1, int* excl2, double* scale,
                  MATRIX3x3* box,int rec_deg,int pbc_deg,
                  double etha,int is_cutoff, double R_on, double R_off,
                  int& time, vector< vector<triple> >& images, vector<triple>& central_translation,
                  double* dr2, double dT, int& is_update
                 ) = &Vdw_LJ;
*/

boost::python::list (*expt_SWITCH)(VECTOR, VECTOR, double, double) = &SWITCH;
boost::python::list (*expt_DOUBLE_SWITCH)(double,double,double) = &DOUBLE_SWITCH;

boost::python::list (*expt_Bond_Harmonic_v1)(VECTOR,VECTOR,double,double) = &Bond_Harmonic; 
boost::python::list (*expt_Bond_Harmonic_v2)(VECTOR,VECTOR,double,double,int) = &Bond_Harmonic; 

boost::python::list (*expt_Bond_Quartic)(VECTOR,VECTOR,double,double) = &Bond_Quartic; 
boost::python::list (*expt_Bond_Morse)(VECTOR,VECTOR,double,double,double) = &Bond_Morse; 

boost::python::list (*expt_Angle_Harmonic)(VECTOR,VECTOR,VECTOR,double,double) = &Angle_Harmonic;
boost::python::list (*expt_Angle_Fourier)(VECTOR,VECTOR,VECTOR,double,double,double,double,int) = &Angle_Fourier;
boost::python::list (*expt_Angle_Fourier_General)(VECTOR,VECTOR,VECTOR,double,double,double,double) = &Angle_Fourier_General;
boost::python::list (*expt_Angle_Fourier_Special)(VECTOR,VECTOR,VECTOR,double,int) = &Angle_Fourier_Special;
boost::python::list (*expt_Angle_Harmonic_Cos)(VECTOR,VECTOR,VECTOR,double,double,int) = &Angle_Harmonic_Cos;
boost::python::list (*expt_Angle_Harmonic_Cos_General)(VECTOR,VECTOR,VECTOR,double,double) = &Angle_Harmonic_Cos_General;
boost::python::list (*expt_Angle_Cubic)(VECTOR,VECTOR,VECTOR,double,double) = &Angle_Cubic;


double (*expt_Elec_Ewald3D_v1)(vector<VECTOR>& r, vector<double>& q, MATRIX3x3& box, double epsilon,
                    vector<VECTOR>& f, MATRIX3x3& at_stress,
                    int rec_deg,int pbc_deg, double etha, double R_on, double R_off   
                   ) = &Elec_Ewald3D;


double (*expt_VdW_Ewald3D_v1)(vector<VECTOR>& r, vector<int>& types, int max_type, vector<double>& Bij, MATRIX3x3& box, /* Inputs */ 
                   vector<VECTOR>& f, MATRIX3x3& at_stress,  /* Outputs*/
                   int rec_deg,int pbc_deg, double etha, double R_on, double R_off    /* Parameters */                   
                   ) = &VdW_Ewald3D;

double (*expt_VdW_Ewald3D_v2)(vector<VECTOR>& r, vector<double>& q, MATRIX3x3& box, 
                   vector<VECTOR>& f, MATRIX3x3& at_stress,  
                   int rec_deg,int pbc_deg, double etha, double R_on, double R_off   
                   ) = &VdW_Ewald3D;




  def("SWITCH", expt_SWITCH);
  def("DOUBLE_SWITCH", expt_DOUBLE_SWITCH);

  def("Bond_Harmonic", expt_Bond_Harmonic_v1);
  def("Bond_Harmonic", expt_Bond_Harmonic_v2);
  def("Bond_Quartic", expt_Bond_Quartic);
  def("Bond_Morse", expt_Bond_Morse);

  def("Angle_Harmonic", expt_Angle_Harmonic);
  def("Angle_Fourier", expt_Angle_Fourier);
  def("Angle_Fourier_General", expt_Angle_Fourier_General);
  def("Angle_Fourier_Special", expt_Angle_Fourier_Special);
  def("Angle_Harmonic_Cos", expt_Angle_Harmonic_Cos);
  def("Angle_Harmonic_Cos_General", expt_Angle_Harmonic_Cos_General);
  def("Angle_Cubic", expt_Angle_Cubic);

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

  def("Elec_Ewald3D", expt_Elec_Ewald3D_v1);
  def("VdW_Ewald3D", expt_VdW_Ewald3D_v1);
  def("VdW_Ewald3D", expt_VdW_Ewald3D_v2);


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

//  export_Mathematics_objects();
//  export_Cell_objects();
  export_Pot_objects();

}

}// namespace libpot
}// liblibra

