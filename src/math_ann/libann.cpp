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

//#include <string>
#if defined(USING_PCH)
#include "../pch.h"
#else

#include <boost/python.hpp>

#endif 

#include "libann.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libann namespace
namespace libann{


void export_NeuralNetwork_objects(){


  void (NeuralNetwork::*expt_init_weights_biases_uniform_v1)
  (Random& rnd, double left_w, double right_w, double left_b, double right_b) = &NeuralNetwork::init_weights_biases_uniform;
  void (NeuralNetwork::*expt_init_weights_biases_normal_v1)
  (Random& rnd, double scaling_w, double shift_w, double scaling_b, double shift_b) = &NeuralNetwork::init_weights_biases_normal;

  vector<MATRIX> (NeuralNetwork::*expt_propagate_v1)(MATRIX& input) = &NeuralNetwork::propagate;
  vector<MATRIX> (NeuralNetwork::*expt_derivatives_v1)(MATRIX& input) = &NeuralNetwork::derivatives;
  double (NeuralNetwork::*expt_back_propagate_v1)(vector<MATRIX>& Y, MATRIX& target) = &NeuralNetwork::back_propagate;
  double (NeuralNetwork::*expt_error_v1)(MATRIX& input, MATRIX& target) = &NeuralNetwork::error;

  vector<double> (NeuralNetwork::*expt_train_v1)
  (Random& rnd, bp::dict params, MATRIX& inputs, MATRIX& targets) = &NeuralNetwork::train;


  void (NeuralNetwork::*expt_save_v1)(std::string filename) = &NeuralNetwork::save;
  int (NeuralNetwork::*expt_load_v1)(std::string filename) = &NeuralNetwork::load;



  class_<NeuralNetwork>("NeuralNetwork",init<>())
      .def(init< vector<int>& >())
      .def(init< std::string >())
      .def("allocate",&NeuralNetwork::allocate)
      .def("show",&NeuralNetwork::show)
      .def("save",expt_save_v1)
      .def("load",expt_load_v1)

      .def("init_weights_biases_uniform",expt_init_weights_biases_uniform_v1)
      .def("init_weights_biases_normal",expt_init_weights_biases_normal_v1)

      .def("propagate",expt_propagate_v1)
      .def("derivatives",expt_derivatives_v1)
      .def("back_propagate",expt_back_propagate_v1)
      .def("error",expt_error_v1)
      .def("train",expt_train_v1)
         
      .def_readwrite("B",&NeuralNetwork::B)
      .def_readwrite("grad_b",&NeuralNetwork::grad_b)
      .def_readwrite("grad_b_old",&NeuralNetwork::grad_b_old)
      .def_readwrite("dB",&NeuralNetwork::dB)
      .def_readwrite("dBold",&NeuralNetwork::dBold)
      .def_readwrite("W",&NeuralNetwork::W)
      .def_readwrite("grad_w",&NeuralNetwork::grad_w)
      .def_readwrite("grad_w_old",&NeuralNetwork::grad_w_old)
      .def_readwrite("dW",&NeuralNetwork::dW)
      .def_readwrite("dWold",&NeuralNetwork::dWold)
      .def_readwrite("Nlayers",&NeuralNetwork::Nlayers)
      .def_readwrite("Npe",&NeuralNetwork::Npe)
      .def_readwrite("sz_x",&NeuralNetwork::sz_x)
      .def_readwrite("sz_y",&NeuralNetwork::sz_y)
       
      .enable_pickling()
  ;

}




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygann){
#else
BOOST_PYTHON_MODULE(libann){
#endif

  //to_python_converter<std::vector<libmmath::DATA>, VecToList<libmmath::DATA> >();

//  export_Mathematics_objects();  // also register mmath python functions!
  export_NeuralNetwork_objects();

}


}// namespace libann
}// namespace liblibra




