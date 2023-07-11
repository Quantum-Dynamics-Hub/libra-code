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

#ifndef NEURAL_NETWORK_H
#define NEURAL_NETWORK_H

#if defined(USING_PCH)
#include "../pch.h"
#else

#include <fstream>
#include <string>
#include <vector>
#include <boost/python.hpp>

#endif 

#include "../math_linalg/liblinalg.h"
#include "../math_data/libdata.h"
#include "../math_random/librandom.h"
#include "../math_specialfunctions/libspecialfunctions.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;
using namespace liblinalg;
using namespace libdata;
using namespace librandom;
using namespace libspecialfunctions;
using namespace std;


/// libann namespace
namespace libann{


class NeuralNetwork{
 
public:

  //-------------------------------------------------
  int sz_x, sz_y; // dimensionality of the input and output
  int Nlayers;
  vector<int> Npe;// Numper of PEs in a particular layers ~= LayerSize;

  // Bias-related
  vector<MATRIX> B;   
  vector<MATRIX> grad_b;      // current gradient:  dError/dB
  vector<MATRIX> grad_b_old;  // previous gradient:  dError/dB
  vector<MATRIX> dB;          // current delta
  vector<MATRIX> dBold;       // previous delta

  // Weight-related
  vector<MATRIX> W;
  vector<MATRIX> grad_w;      // gradient:  dError/dW
  vector<MATRIX> grad_w_old;  // previous gradient:  dError/dB
  vector<MATRIX> dW;          // current delta
  vector<MATRIX> dWold;       // previous delta

//  vector<MATRIX> NET;  // net values    NET[L] = W[L]*Y[L-1] + B[L]; where X is the input, the output of the layer L-1
//  vector<MATRIX> TF;   // transfer function   e.g.  Y[L] = tanh(NET[L]), element by element 
//  vector<MATRIX> dF;   // derivative of the transfer functions, e.g. dF[L] = 1- Y[L]*Y[L], element by element
  vector<MATRIX> D;
  vector<MATRIX> Delta;
 
  //--------------- Basic methods: NeuralNetwork.cpp ---------------------
  // Auxiliary
  void allocate(vector<int>& arch);

  // Constructors
  NeuralNetwork(){ Nlayers = 0;  }
  NeuralNetwork(vector<int>& arch);
  NeuralNetwork(std::string xml_filename);

  // Copy constructor
  NeuralNetwork(const NeuralNetwork&); 

  // Assignment
  NeuralNetwork& operator=(const NeuralNetwork&);

  // Destructor
  ~NeuralNetwork(){ ;; }

  // Printing out 
  int show();  

  // Saving/Loading to a property tree
  void save(boost::property_tree::ptree& pt,std::string path);
  void save(std::string filename);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);
  int load(std::string filename);


  //--------------- Simulation methods: NeuralNetwork_Algorithms.cpp ---------------------

  // Init weights and biases
  void init_weights_biases_uniform(Random& rnd, double left_w, double right_w, double left_b, double right_b);
  void init_weights_biases_normal(Random& rnd, double scaling_w, double shift_w, double scaling_b, double shift_b);

  // Basic calculations
  vector<MATRIX> propagate(MATRIX& input);
  vector<MATRIX> derivatives(MATRIX& input);

  // Back propagate - compute gradients w.r.t. weights and biases
  double back_propagate(vector<MATRIX>& Y, MATRIX& target);
  double error(MATRIX& input, MATRIX& target);

  // Training
  vector<double> train(Random& rnd, bp::dict params, MATRIX& inputs, MATRIX& targets);


};

void save(boost::property_tree::ptree& pt,std::string path,vector<NeuralNetwork>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<NeuralNetwork>& vt,int& status);


}// namespace libann
}// namespace liblibra

#endif
