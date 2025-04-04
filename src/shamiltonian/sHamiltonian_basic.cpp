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
  \file sHamiltonian_basic.cpp
  \brief The file implements the generic methods of the sHamiltonian class: getters and setters
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <stdlib.h>
#endif 

//#include <torch/torch.h>
#include "sHamiltonian.h"
#include <torch/csrc/api/include/torch/python.h>

#include "../math_meigen/mEigen.h"

namespace bp = boost::python;

/// liblibra namespace
namespace liblibra{

using namespace libmeigen;

/// libnhamiltonian namespace 
namespace libshamiltonian{



bp::object tensor_to_python(const torch::Tensor& tensor) {
    // Convert torch::Tensor to a NumPy array, then to a Python object
    PyObject* array = THPVariable_Wrap(tensor);
    return bp::object(bp::handle<>(array));
}

torch::Tensor python_to_tensor(const bp::object& obj) {
    PyObject* py_obj = obj.ptr();
    if (!THPVariable_Check(py_obj)) {
        throw std::runtime_error("Expected a PyTorch tensor.");
    }
    return THPVariable_Unpack(py_obj);
}




sHamiltonian::sHamiltonian(const sHamiltonian& src){
  nbeads = src.nbeads;
  nel_dof = src.nel_dof;
  nn_dof = src.nn_dof;

  ovlp_dia = src.ovlp_dia;
}


sHamiltonian::sHamiltonian(int nbeads_, int nelec_, int nnucl_){

  nbeads = nbeads_;
  nel_dof = nelec_; 
  nn_dof = nnucl_;
 
  ovlp_dia = torch::zeros({ nbeads, nel_dof, nn_dof });

}



sHamiltonian::~sHamiltonian(){ 
  ;;
}


boost::python::object sHamiltonian::get_ovlp_dia(){

  return tensor_to_python(ovlp_dia);

}



}// namespace libshamiltonian
}// liblibra

