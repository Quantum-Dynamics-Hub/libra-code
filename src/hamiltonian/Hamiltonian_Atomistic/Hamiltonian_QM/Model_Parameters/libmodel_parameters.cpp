#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libmodel_parameters.h"


using namespace boost::python;

namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_qm{
namespace libmodel_parameters{


void export_Model_Parameters_objects(){

  class_<HF_integrals>("HF_integrals",init<>())
      .def("set_JK_values", &HF_integrals::set_JK_values)
      .def("get_JK_values", &HF_integrals::get_JK_values)

  ;

  class_<EHT_K>("EHT_K",init<>())
      .def("set_PSPS_value", &EHT_K::set_PSPS_value)

      .def("set_PPa_value", &EHT_K::set_PPa_value)
      .def("set_PP0_value", &EHT_K::set_PP0_value)
      .def("set_PP1_value", &EHT_K::set_PP1_value)
      .def("set_PP2_value", &EHT_K::set_PP2_value)
      .def("get_PPa_value", &EHT_K::get_PPa_value)
      .def("get_PP0_value", &EHT_K::get_PP0_value)
      .def("get_PP1_value", &EHT_K::get_PP1_value)
      .def("get_PP2_value", &EHT_K::get_PP2_value)

      .def("set_K_value",  &EHT_K::set_K_value)
      .def("set_K1_value", &EHT_K::set_K1_value)
      .def("set_K2_value", &EHT_K::set_K2_value)
      .def("set_K4_value", &EHT_K::set_K4_value)
      .def("get_K_value",  &EHT_K::get_K_value)
      .def("get_K1_value", &EHT_K::get_K1_value)
      .def("get_K2_value", &EHT_K::get_K2_value)
      .def("get_K4_value", &EHT_K::get_K4_value)

      .def("set_C0_value", &EHT_K::set_C0_value)
      .def("set_C1_value", &EHT_K::set_C1_value)
      .def("set_C2_value", &EHT_K::set_C2_value)
      .def("set_C3_value", &EHT_K::set_C3_value)
      .def("set_C4_value", &EHT_K::set_C4_value)
      .def("get_C0_value", &EHT_K::get_C0_value)
      .def("get_C1_value", &EHT_K::get_C1_value)
      .def("get_C2_value", &EHT_K::get_C2_value)
      .def("get_C3_value", &EHT_K::get_C3_value)
      .def("get_C4_value", &EHT_K::get_C4_value)

  ;

  class_<mEHT_K>("mEHT_K",init<>())

      .def_readwrite("size",   &mEHT_K::size)
      .def_readwrite("eht_K",  &mEHT_K::eht_K)
      .def_readwrite("eht_K1", &mEHT_K::eht_K1)
      .def_readwrite("eht_K2", &mEHT_K::eht_K2)
      .def_readwrite("eht_K3", &mEHT_K::eht_K3)
      .def_readwrite("eht_K4", &mEHT_K::eht_K4)

      .def_readwrite("eht_C0", &mEHT_K::eht_C0)
      .def_readwrite("eht_C1", &mEHT_K::eht_C1)
      .def_readwrite("eht_C2", &mEHT_K::eht_C2)
      .def_readwrite("eht_C3", &mEHT_K::eht_C3)
      .def_readwrite("eht_C4", &mEHT_K::eht_C4)

      .def_readwrite("eht_PPa", &mEHT_K::eht_PPa)
      .def_readwrite("eht_PP0", &mEHT_K::eht_PP0)
      .def_readwrite("eht_PP1", &mEHT_K::eht_PP1)
      .def_readwrite("eht_PP2", &mEHT_K::eht_PP2)


      .def("set_mapping", &mEHT_K::set_mapping)
      .def("set_mapping1", &mEHT_K::set_mapping1)

      .def("get_K_value",  &mEHT_K::get_K_value)
      .def("get_K1_value", &mEHT_K::get_K1_value)
      .def("get_K2_value", &mEHT_K::get_K2_value)
      .def("get_K3_value", &mEHT_K::get_K3_value)
      .def("get_K4_value", &mEHT_K::get_K4_value)

      .def("get_C0_value", &mEHT_K::get_C0_value)
      .def("get_C1_value", &mEHT_K::get_C1_value)
      .def("get_C2_value", &mEHT_K::get_C2_value)
      .def("get_C3_value", &mEHT_K::get_C3_value)
      .def("get_C4_value", &mEHT_K::get_C4_value)

  ;


  class_<Element>("Element",init<>())
//      .def(_s)
//  void _set(std::string en,int z){ elt_name = en; Z = z;  }
//  void _set(std::string en,int z,int nv){ elt_name = en; Z = z; Nval = nv; }
//  void _set(std::string en,int z,int nv,double zeff){ elt_name = en; Z = z; Nval = nv; Zeff = zeff; }
//  void _set(std::string en,int z,int nv,double zeff,map<std::string, double> ip){ elt_name = en; Z = z; Nval = nv; IP = ip; Zeff = zeff;}
//  void set_mass(double m_){ mass = m_; }


      .def_readwrite("elt_name",   &Element::elt_name)
      .def_readwrite("Z",   &Element::Z)
      .def_readwrite("PQN",   &Element::PQN)
      .def_readwrite("Nval",   &Element::Nval)
      .def_readwrite("Zeff",   &Element::Zeff)
      .def_readwrite("mass",   &Element::mass)

      .def_readwrite("IP",   &Element::IP)
      .def_readwrite("EA",   &Element::EA)

      .def_readwrite("Nquant",   &Element::Nquant)
      .def_readwrite("Nzeta",   &Element::Nzeta)
      .def_readwrite("zetas",   &Element::zetas)
      .def_readwrite("coeffs",   &Element::coeffs)

      .def_readwrite("J_param1",   &Element::J_param1)
      .def_readwrite("J_param2",   &Element::J_param2)
      .def_readwrite("J_param3",   &Element::J_param3)
      .def_readwrite("J_param4",   &Element::J_param4)

      .def_readwrite("G1",   &Element::G1)
      .def_readwrite("F2",   &Element::F2)
      .def_readwrite("beta0",  &Element::beta0)
      .def_readwrite("Zeta",   &Element::Zeta)

  ;





  class_<Control_Parameters>("Control_Parameters",init<>())
      .def("__copy__", &generic__copy__<Control_Parameters>)
      .def("__deepcopy__", &generic__deepcopy__<Control_Parameters>)

      .def_readwrite("runtype", &Control_Parameters::runtype)
      .def_readwrite("hamiltonian", &Control_Parameters::hamiltonian)
      .def_readwrite("spin_method", &Control_Parameters::spin_method)
      .def_readwrite("DF", &Control_Parameters::DF)

      .def_readwrite("guess_type", &Control_Parameters::guess_type)

      .def_readwrite("scf_algo", &Control_Parameters::scf_algo)
      .def_readwrite("use_disk", &Control_Parameters::use_disk)
      .def_readwrite("use_rosh", &Control_Parameters::use_rosh)
      .def_readwrite("do_annihilate", &Control_Parameters::do_annihilate)
      .def_readwrite("pop_opt", &Control_Parameters::pop_opt)
      .def_readwrite("use_diis", &Control_Parameters::use_diis)
      .def_readwrite("diis_max", &Control_Parameters::diis_max)
      .def_readwrite("diis_start_iter", &Control_Parameters::diis_start_iter)
      .def_readwrite("use_level_shift", &Control_Parameters::use_level_shift)
      .def_readwrite("shift_magnitude", &Control_Parameters::shift_magnitude)
      .def_readwrite("use_damping", &Control_Parameters::use_damping)
      .def_readwrite("damping_start", &Control_Parameters::damping_start)
      .def_readwrite("damping_const", &Control_Parameters::damping_const)
      .def_readwrite("etol", &Control_Parameters::etol)
      .def_readwrite("den_tol", &Control_Parameters::den_tol)
      .def_readwrite("Niter", &Control_Parameters::Niter)
      .def_readwrite("degen_tol", &Control_Parameters::degen_tol)

      .def_readwrite("parameters", &Control_Parameters::parameters)
      .def_readwrite("eht_params_format", &Control_Parameters::eht_params_format)
      .def_readwrite("eht_formula", &Control_Parameters::eht_formula)
      .def_readwrite("eht_sce_formula", &Control_Parameters::eht_sce_formula)
      .def_readwrite("eht_fock_opt", &Control_Parameters::eht_fock_opt)
      .def_readwrite("eht_electrostatics", &Control_Parameters::eht_electrostatics)


      .def_readwrite("compute_vertical_ip", &Control_Parameters::compute_vertical_ip)
      .def_readwrite("compute_vertical_ea", &Control_Parameters::compute_vertical_ea)

      .def_readwrite("md_dt", &Control_Parameters::md_dt)
      .def_readwrite("md_nsteps", &Control_Parameters::md_nsteps)

      .def_readwrite("opt_dt", &Control_Parameters::opt_dt)
      .def_readwrite("opt_nsteps", &Control_Parameters::opt_nsteps)

      .def_readwrite("compute_dipole", &Control_Parameters::compute_dipole)

      .def_readwrite("compute_dos", &Control_Parameters::compute_dos)
      .def_readwrite("dos_opt", &Control_Parameters::dos_opt)
      .def_readwrite("dos_prefix", &Control_Parameters::dos_prefix)

      .def_readwrite("compute_charge_density", &Control_Parameters::compute_charge_density)
      .def_readwrite("nx_grid", &Control_Parameters::nx_grid)
      .def_readwrite("ny_grid", &Control_Parameters::ny_grid)
      .def_readwrite("nz_grid", &Control_Parameters::nz_grid)
      .def_readwrite("charge_density_prefix", &Control_Parameters::charge_density_prefix)
      .def_readwrite("orbs", &Control_Parameters::orbs)

      .def_readwrite("nac_md_trajectory_filename", &Control_Parameters::nac_md_trajectory_filename)
      .def_readwrite("nac_prefix", &Control_Parameters::nac_prefix)
      .def_readwrite("nac_min_frame", &Control_Parameters::nac_min_frame)
      .def_readwrite("nac_max_frame", &Control_Parameters::nac_max_frame)
      .def_readwrite("nac_min_orbs", &Control_Parameters::nac_min_orbs)
      .def_readwrite("nac_max_orbs", &Control_Parameters::nac_max_orbs)
      .def_readwrite("nac_dt", &Control_Parameters::nac_dt)
      .def_readwrite("nac_opt", &Control_Parameters::nac_opt)

      .def_readwrite("scan_mov_at", &Control_Parameters::scan_mov_at)
      .def_readwrite("scan_ref_at", &Control_Parameters::scan_ref_at)
      .def_readwrite("scan_dir", &Control_Parameters::scan_dir)
      .def_readwrite("scan_dxmin", &Control_Parameters::scan_dxmin)
      .def_readwrite("scan_dxmax", &Control_Parameters::scan_dxmax)
      .def_readwrite("scan_dx", &Control_Parameters::scan_dx)


      .def_readwrite("compute_excitations", &Control_Parameters::compute_excitations)
      .def_readwrite("num_excitations", &Control_Parameters::num_excitations)
      .def_readwrite("excitations_opt", &Control_Parameters::excitations_opt)
      .def_readwrite("spectral_width", &Control_Parameters::spectral_width)
      .def_readwrite("excitations", &Control_Parameters::excitations)  

      .def_readwrite("t1", &Control_Parameters::t1)
      .def_readwrite("t2", &Control_Parameters::t2)
      .def_readwrite("t3", &Control_Parameters::t3)
      .def_readwrite("x_period", &Control_Parameters::x_period)
      .def_readwrite("y_period", &Control_Parameters::y_period)
      .def_readwrite("z_period", &Control_Parameters::z_period)

      .def_readwrite("Natoms", &Control_Parameters::Natoms)
      .def_readwrite("charge", &Control_Parameters::charge)
      .def_readwrite("spin", &Control_Parameters::spin)
      .def_readwrite("coordinates", &Control_Parameters::coordinates)

      .def_readwrite("fragments", &Control_Parameters::fragments)
      .def_readwrite("frag_size", &Control_Parameters::frag_size)
      .def_readwrite("frag_name", &Control_Parameters::frag_name)
      .def_readwrite("frag_charge", &Control_Parameters::frag_charge)

  ;


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmodel_parameters){
#else
BOOST_PYTHON_MODULE(libmodel_parameters){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Model_Parameters_objects();

}



}// namespace libmodel_parameters
}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


