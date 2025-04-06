/*********************************************************************************
* Copyright (C) 2025 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libshamiltonian.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include <ATen/ATen.h>

#include "libshamiltonian.h"

/// liblibra namespace
namespace liblibra{


using namespace boost::python;
namespace py = boost::python;

/// libshamiltonian namespace
namespace libshamiltonian{


// Convert Python torch.Tensor -> C++ torch::Tensor
struct TorchTensorConverter {
    static void* convertible(PyObject* obj) {
        if (!THPVariable_Check(obj)) {
            return nullptr;
        }
        return obj;
    }

    static void construct(PyObject* obj, py::converter::rvalue_from_python_stage1_data* data) {
        // Get the underlying at::Tensor
        at::Tensor tensor = THPVariable_Unpack(obj);
        
        // Allocate storage and construct the tensor
        void* storage = (
            (py::converter::rvalue_from_python_storage<torch::Tensor>*)
            data
        )->storage.bytes;
        
        new (storage) torch::Tensor(std::move(tensor));
        data->convertible = storage;
    }
};

// Convert C++ torch::Tensor -> Python torch.Tensor
struct TorchTensorToPython {
    static PyObject* convert(const torch::Tensor& tensor) {
        return THPVariable_Wrap(tensor);
    }
};


torch::Tensor process_tensor(torch::Tensor input) {
    return input.mul(2).add(1);  // Example: y = 2x + 1
}

void export_shamiltonian_objects(){
/** 
  \brief Exporter of the libshamiltonian_generic classes and functions

*/

  def("process_tensor", &process_tensor, "Doubles input tensor and adds 1");

  class_<torch::Tensor>("Tensor",init<>())
      .def("__copy__", &generic__copy__<torch::Tensor>)
     .def("__deepcopy__", &generic__deepcopy__<torch::Tensor>)
  ;



  bp::object (*expt_tensor_to_python_v1)(const torch::Tensor& tensor) = &tensor_to_python;
  def("cpp2py", expt_tensor_to_python_v1);

  torch::Tensor (*expt_python_to_tensor_v1)(const bp::object& obj) = &python_to_tensor;
  def("py2cpp", expt_python_to_tensor_v1);


 void (sHamiltonian::*expt_compute_v1)(std::string property, bp::object py_funct, torch::Tensor q, bp::object params)
  = &sHamiltonian::compute;

 void (sHamiltonian::*expt_compute_v2)(bp::object py_funct, torch::Tensor q, bp::object params)
  = &sHamiltonian::compute;


  class_<sHamiltonian>("sHamiltonian",init<int,int,int>())
      .def(init<const sHamiltonian&>())
//      .def("__copy__", &generic__copy__<Hamiltonian>)
//      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian>)
      .def_readwrite("nbeads", &sHamiltonian::nbeads)
      .def_readwrite("nstates", &sHamiltonian::nstates)
      .def_readwrite("nnucl", &sHamiltonian::nnucl)
      .def_readwrite("ovlp_dia", &sHamiltonian::ovlp_dia)
      .def_readwrite("ham_dia", &sHamiltonian::ham_dia)
      .def_readwrite("nac_dia", &sHamiltonian::nac_dia)
      .def_readwrite("hvib_dia", &sHamiltonian::hvib_dia)
      .def_readwrite("dc1_dia", &sHamiltonian::dc1_dia)
      .def_readwrite("d1ham_dia", &sHamiltonian::d1ham_dia)
      .def_readwrite("d2ham_dia", &sHamiltonian::d2ham_dia)

      .def_readwrite("ham_adi", &sHamiltonian::ham_adi)
      .def_readwrite("nac_adi", &sHamiltonian::nac_adi)
      .def_readwrite("hvib_adi", &sHamiltonian::hvib_adi)
      .def_readwrite("dc1_adi", &sHamiltonian::dc1_adi)
      .def_readwrite("d1ham_adi", &sHamiltonian::d1ham_adi)
      .def_readwrite("d2ham_adi", &sHamiltonian::d2ham_adi)

      .def_readwrite("basis_transform", &sHamiltonian::basis_transform)
      .def_readwrite("time_overlap_dia", &sHamiltonian::time_overlap_dia)
      .def_readwrite("time_overlap_adi", &sHamiltonian::time_overlap_adi)

      .def("bind", &sHamiltonian::bind)
      //.def("get_tensor", &sHamiltonian::get_tensor)

      .def("compute", expt_compute_v1)
      .def("compute", expt_compute_v2)
      .def("dia2adi", &sHamiltonian::dia2adi)
  
  ;

  // Python -> C++
     py::converter::registry::push_back(
        &TorchTensorConverter::convertible,
        &TorchTensorConverter::construct,
        py::type_id<torch::Tensor>()
    );
    
  // C++ -> Python
    py::to_python_converter<torch::Tensor, TorchTensorToPython>();


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygshamiltonian){
#else
BOOST_PYTHON_MODULE(libshamiltonian){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_shamiltonian_objects();

}



}// namespace libshamiltonian
}// liblibra

