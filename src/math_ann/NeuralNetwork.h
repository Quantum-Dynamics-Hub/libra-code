/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef NEURAL_NETWORK_H
#define NEURAL_NETWORK_H

#include <fstream>
#include <string>
#include <vector>
#include <boost/python.hpp>

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

//----------------------------------------------
class ANNData{

public:
  // This is how the input data (training patterns) should
  // look like:
  // This is single input-output pair of vectors

  vector<double> Input;    // may be of size n
  vector<double> Output;   // may be of size m
  vector<double> Derivs;   // must be of size n*m 

  ANNData& operator=(const ANNData&);

};
//----------------------------------------------



class NeuralNetwork{
     
  int ScaleDerivatives();
 
public:

  // --------------- Parameters of NN simulations ----------------

  std::string learning_method;     int is_learning_method;

  double learning_rate;            int is_learning_rate;
  double momentum_term;            int is_momentum_term;
  int epoch_size;                  int is_epoch_size;
  int iterations_in_cycle;         int is_iterations_in_cycle;
  double grad_weight;              int is_grad_weight;
  vector<double> weight_decay;     int is_weight_decay;
  double norm_exp;                 int is_norm_exp;
  double a_plus;                   int is_a_plus;
  double a_minus;                  int is_a_minus;

  std::string scale_method;
  int Iteration;
  int Cycle;
  int derivs_flag;

  //---------- Data and its properties ------------
  vector<ANNData> TrainData; // Thus TrainData[i] - is i-th pattern
  ANNData         Recall;
  vector<DATA>    Inputs;
  vector<DATA>    Outputs;
  vector<DATA>    Derivs;
    
  int num_of_patterns; //This is the total size of all data patterns
  int sz_x;
  int sz_y;
  // And here are the corresponding scaling parameters
  // for all inputs (x_) and outputs (y_)
  int sz_d; // size of derivatives

  //-------------------------------------------------
  int Nlayers;
  vector<int> Npe;// Numper of PEs in a particular layers ~= LayerSize;

  vector<MATRIX> B;
  vector<MATRIX> dB; // for momentum term
  vector<MATRIX> dBcurr; // for batch
  vector<MATRIX> dBold;
  vector<MATRIX> W;
  vector<MATRIX> dW; // for momentum term
  vector<MATRIX> dWcurr; // for batch
  vector<MATRIX> dWold; 
//  vector<MATRIX> NET;  // net values    NET[L] = W[L]*Y[L-1] + B[L]; where X is the input, the output of the layer L-1
//  vector<MATRIX> TF;   // transfer function   e.g.  Y[L] = tanh(NET[L]), element by element 
//  vector<MATRIX> dF;   // derivative of the transfer functions, e.g. dF[L] = 1- Y[L]*Y[L], element by element
  vector<MATRIX> D;
  vector<MATRIX> Delta;
 
  //--------------- Methods ---------------------
  // Default constructor
  NeuralNetwork(){           

   Iteration = 0;
   Cycle     = 0;

   is_learning_method = 0;
   is_learning_rate   = 0;
   is_momentum_term   = 0;
   is_epoch_size      = 0;
   is_iterations_in_cycle  = 0;
   is_grad_weight     = 0;
   is_weight_decay    = 0;
   is_norm_exp        = 0;
   is_a_plus          = 0;
   is_a_minus         = 0;
   derivs_flag = 0;

   scale_method = "none";
  }

  NeuralNetwork(vector<int>& arch);


  // Copy constructor
  NeuralNetwork(const NeuralNetwork&); 

  // Destructor
  ~NeuralNetwork(); 

  // Basic calculations
  vector<MATRIX> propagate(MATRIX& input);
  vector<MATRIX> derivatives(MATRIX& input);
  double back_propagate(vector<MATRIX>& Y, MATRIX& target);

  // Init weights and biases
  void init_weights_biases_uniform(Random& rnd, double left_w, double right_w, double left_b, double right_b);
  void init_weights_biases_normal(Random& rnd, double scaling_w, double shift_w, double scaling_b, double shift_b);


  // Training
  void train(Random& rnd, bp::dict params, MATRIX& inputs, MATRIX& targets);


  // Methods
  void CreateANN(boost::python::list);
  int ShowANN();
  
  void set(boost::python::object obj);
  void set(boost::python::dict d);

  int ExportANN(std::string filename);
  int ImportANN(std::string filename);

  void save(boost::property_tree::ptree& pt,std::string path);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);

  NeuralNetwork& operator=(const NeuralNetwork&);

  int ClearTrainingData();
  int AddTrainingData(boost::python::list patterns);


  int SetTrainingData(object,int);
  int ScaleTrainingData(int,int);
  int ScaleTrainingData(int,int,boost::python::list,boost::python::list);
  int NormalizeTrainingData(int,int);
  int NormalizeAndScaleTrainingData(int,int,boost::python::list minlim,boost::python::list maxlim);
  int NormalizeAndTransformTrainingData(int,int);
  int CropTrainingData(boost::python::list,boost::python::list);

  //void set_parameters(boost::python::dict params);

  void ANNTrain();
//  int Propagate(const MATRIX& input, MATRIX& result);
//  int Propagate(boost::python::list input, boost::python::list& result);
  int Propagate(const MATRIX& input, MATRIX& result, MATRIX& derivs);  
  int Propagate(boost::python::list input, boost::python::list& result, boost::python::list& derivs);
  void LearningHistory(std::string filename,std::string data_flag);


};

void save(boost::property_tree::ptree& pt,std::string path,vector<NeuralNetwork>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<NeuralNetwork>& vt,int& status);


}// namespace libann
}// namespace liblibra

#endif
