#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libensemble.h"

using namespace boost::python;


namespace libdyn{
namespace libensemble{


void export_Ensemble_objects(){

//  def("SAC_Ham", expt_SAC_Ham1);
//  void (Electronic::*expt_propagate_electronic)(double,Hamiltonian&) = &Electronic::propagate_electronic;

  void (Ensemble::*expt_ham_set_v_v1)(int i, vector<double>& v) = &Ensemble::ham_set_v;
  void (Ensemble::*expt_ham_set_v_v2)(int i, boost::python::list v) = &Ensemble::ham_set_v;
  void (Ensemble::*expt_ham_set_v_v3)() = &Ensemble::ham_set_v;

  void (Ensemble::*expt_ham_set_ham_v1)(int i, std::string opt, int mopt) = &Ensemble::ham_set_ham;
  void (Ensemble::*expt_ham_set_ham_v2)(std::string opt, int mopt) = &Ensemble::ham_set_ham;
  void (Ensemble::*expt_ham_set_ham_v3)(int i, Hamiltonian& _ham) = &Ensemble::ham_set_ham;

  void (Ensemble::*expt_ham_set_rep_v1)(int i, int _rep) = &Ensemble::ham_set_rep;
  void (Ensemble::*expt_ham_set_rep_v2)(int _rep) = &Ensemble::ham_set_rep;

  void (Ensemble::*expt_el_propagate_electronic_v1)(int i, double dt) = &Ensemble::el_propagate_electronic;
  void (Ensemble::*expt_el_propagate_electronic_v2)(double dt) = &Ensemble::el_propagate_electronic;
  void (Ensemble::*expt_mol_propagate_q_v1)(int i, double dt) = &Ensemble::mol_propagate_q;
  void (Ensemble::*expt_mol_propagate_q_v2)(double dt) = &Ensemble::mol_propagate_q;
  void (Ensemble::*expt_mol_propagate_p_v1)(int i, double dt) = &Ensemble::mol_propagate_p;
  void (Ensemble::*expt_mol_propagate_p_v2)(double dt) = &Ensemble::mol_propagate_p;

  void (Ensemble::*expt_se_pop_v1)(vector<double>& v) = &Ensemble::se_pop;
  boost::python::list (Ensemble::*expt_se_pop_v2)() = &Ensemble::se_pop;

  void (Ensemble::*expt_sh_pop_v1)(vector<double>& v) = &Ensemble::sh_pop;
  boost::python::list (Ensemble::*expt_sh_pop_v2)() = &Ensemble::sh_pop;





  class_<Ensemble>("Ensemble",init<>())
      .def(init<int,int,int>())
      .def("__copy__", &generic__copy__<Ensemble>)
      .def("__deepcopy__", &generic__deepcopy__<Ensemble>)

      .def_readwrite("ntraj", &Ensemble::ntraj)
      .def_readwrite("nnucl", &Ensemble::nnucl)
      .def_readwrite("nelec", &Ensemble::nelec)
      .def_readwrite("is_active", &Ensemble::is_active)

      .def_readwrite("mol", &Ensemble::mol)
      .def_readwrite("el", &Ensemble::el)
//      .def_readwrite("ham", &Ensemble::ham)

//      .def("get_mol", &Ensemble::get_mol)
//      .def("get_el", &Ensemble::get_el)
//      .def("get_ham", &Ensemble::get_ham, return_value_policy<manage_new_object>())
//      .def("get_ham", &Ensemble::get_ham, return_value_policy<reference_existing_object>())
//      .def("set_mol", &Ensemble::set_mol)
//      .def("set_el", &Ensemble::set_el)
      .def("ham_set_ham", expt_ham_set_ham_v1)
      .def("ham_set_ham", expt_ham_set_ham_v2)
      .def("ham_set_ham", expt_ham_set_ham_v3)

      .def("ham_set_rep", expt_ham_set_rep_v1)
      .def("ham_set_rep", expt_ham_set_rep_v2)

      .def("ham_set_v", expt_ham_set_v_v1)
      .def("ham_set_v", expt_ham_set_v_v2)
      .def("ham_set_v", expt_ham_set_v_v3)

      .def("el_propagate_electronic", expt_el_propagate_electronic_v1)
      .def("el_propagate_electronic", expt_el_propagate_electronic_v2)

      .def("mol_propagate_q", expt_mol_propagate_q_v1)
      .def("mol_propagate_q", expt_mol_propagate_q_v2)
      .def("mol_propagate_p", expt_mol_propagate_p_v1)
      .def("mol_propagate_p", expt_mol_propagate_p_v2)


      .def("se_pop", expt_se_pop_v1)
      .def("se_pop", expt_se_pop_v2)
      .def("sh_pop", expt_sh_pop_v1)
      .def("sh_pop", expt_sh_pop_v2)

 
  ;




}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygensemble){
#else
BOOST_PYTHON_MODULE(libensemble){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Ensemble_objects();

}


}// namespace libensemble
}// namespace libdyn

