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
  def("t2p", expt_tensor_to_python_v1);

  torch::Tensor (*expt_python_to_tensor_v1)(const bp::object& obj) = &python_to_tensor;
  def("p2t", expt_python_to_tensor_v1);



  class_<sHamiltonian>("sHamiltonian",init<int,int,int>())
      .def(init<const sHamiltonian&>())
//      .def("__copy__", &generic__copy__<Hamiltonian>)
//      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian>)
      .def_readwrite("nbeads", &sHamiltonian::nbeads)
      .def_readwrite("nel_dof", &sHamiltonian::nel_dof)
      .def_readwrite("nn_dof", &sHamiltonian::nn_dof)
      .def_readwrite("eigen_algo", &sHamiltonian::eigen_algo)
      .def_readwrite("phase_corr_ovlp_tol", &sHamiltonian::phase_corr_ovlp_tol)
      .def_readwrite("ovlp_dia", &sHamiltonian::ovlp_dia)
      .def("get_ovlp_dia", &sHamiltonian::get_ovlp_dia)
      .def("set_tensor", &sHamiltonian::set_tensor)
      .def("get_tensor", &sHamiltonian::get_tensor)
  
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

