/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "NeuralNetwork.h"

/// liblibra namespace
namespace liblibra{

using namespace boost;

/// libann namespace
namespace libann{


vector<MATRIX> NeuralNetwork::propagate(MATRIX& input){
/**
  This is a new implementation of the ANN propagation:

  This function computes the output in each layer of the ANN,
  the one in the last layer is regarded as the predicted overall output = 
  the approximate of the target properties

  Args:
    input  - sz_input x num_patterns, each column is a sz_x - dimensional input

  Returns: 
    output - Nlayers matrices of the  Npe[L] x num_patterns size each 

    Predicted output in every layer


  Notation is from here:
    https://cnl.salk.edu/~schraudo/teach/NNcourse/backprop.html

*/

  int i, j, L;  

  if(input.n_rows!=sz_x){
    std::cout<<"Error: Size of the input "<<input.n_rows<<" does not match the ANN architecture "<<sz_x<<std::endl;
    exit(0);
  }
  int sz = input.n_cols; // number of patterns to handle at the same time


  /**

  L       0                   1                       ....             NL = Nlayers - 1

  W      [junk]              W[1]                                        W[NL]

  B      [junk]              B[1]                                        B[NL]

  Y      [Y[0]=input]      [ f(W[1]*Y[0] + B[1]) ]           [ output = f(W[NL]*Y[NL-1] + B[NL]) ]
  
  */

  vector<MATRIX> Y;

  /// L = 0
  Y.push_back(input);

  /// L = 1, 2... Nlayers -1
  for(L = 1; L < Nlayers; L++){

    Y.push_back(MATRIX(Npe[L], sz)); 
    Y[L] = W[L] * Y[L-1];

    /// Transfer function and bias
    for(i=0; i<Npe[L]; i++){ 
      for(j=0; j<sz; j++){ 
        Y[L].set(i, j, tanh( Y[L].get(i, j) +  B[L].get(i, 0)) );
      }
    }
    
  }// for L

  return Y;

}



vector<MATRIX> NeuralNetwork::derivatives(MATRIX& input){
/**
  This is a new implementation of the derivatives of the ANN outputs w.r.t. its inputs:

  Args:
    input  - sz_input x num_patterns, each column is a sz_x - dimensional input

  Returns: 
    output - num_patters matrices of the  Npe[L] Npe[0] size each , output.get(i,j) = doutput[j]/dinput[i]

  Notation is from here:
    https://cnl.salk.edu/~schraudo/teach/NNcourse/backprop.html

*/

  int i, j, p, L;  

  if(input.n_rows!=sz_x){
    std::cout<<"Error: Size of the input "<<input.n_rows<<" does not match the ANN architecture "<<sz_x<<std::endl;
    exit(0);
  }
  int sz = input.n_cols; // number of patterns to handle at the same time


  /**

  L       0                   1                       ....             NL = Nlayers - 1

  W      [junk]              W[1]                                        W[NL]

  B      [junk]              B[1]                                        B[NL]

  Y      [Y[0]=input]      [ f(W[1]*Y[0] + B[1]) ]           [ output = f(W[NL]*Y[NL-1] + B[NL]) ]
  
  */

  vector<MATRIX> res;


  vector<MATRIX> Y = propagate(input);

  /// ksi[L].get(i,j) - is basically dx_L[i]/dx_0[j] - temporary variables for a given pattern
  vector<MATRIX> ksi;

  for(L = 0; L < Nlayers; L++){  
    ksi.push_back(MATRIX(Npe[L], sz_x));
  }// for L

  /// Repeat for all patterns 
  for(p = 0; p < sz; p++){

    // Initialize
    ksi[0].Init_Unit_Matrix(1.0);  
    
    for(L = 1; L < Nlayers; L++){            

      MATRIX D_L(Npe[L],Npe[L]); 
      for(i = 0; i<Npe[L]; i++){  D_L.set(i,i, 1.0 - Y[L].get(i, p) * Y[L].get(i, p));   }

      ksi[L] = D_L * W[L] * ksi[L-1]; 

    }// for L

    res.push_back(ksi[Nlayers-1]);

  }// for p

  return res;

}




double NeuralNetwork::back_propagate(vector<MATRIX>& Y, MATRIX& target){
/**
  This is a new implementation of the ANN back propagation

  Args:
    Y - Nlayers matrices of the  Npe[L] x num_patterns size each. Predicted output in every layer
    target  Npe[Nlayers-1] x num_patterns - the expected output of the ANN

  Returns: 
    update dW: Nlayers matrices of the  Npe[L] x Npe[L-1] - the gradients of the weights  - average over all patterns
    update dB: Nlayers matrices of the  Npe[L] x 1 - the gradients of the biases  - average over all patterns

  Notation is from here:
    https://cnl.salk.edu/~schraudo/teach/NNcourse/backprop.html

*/

  int i, j, L;  

  if(target.n_rows!=sz_y){
    std::cout<<"Error: Size of the target output "<<target.n_rows<<" does not match the ANN architecture "<<sz_y<<std::endl;
    exit(0);
  }
  int sz = target.n_cols; // number of patterns to handle at the same time


  /**

  L       0                   1                       ....             NL = Nlayers - 1

  W      [junk]              W[1]                                        W[NL]

  B      [junk]              B[1]                                        B[NL]

  Y      [Y[0]=input]      [ f(W[1]*Y[0] + B[1]) ]           [ output = f(W[NL]*Y[NL-1] + B[NL]) ]

 deltas  [junk]           W^T[2]*delta[2] *f'(Y[1])                   target - output[NL]
  
  */


  vector<MATRIX> delta;

  /// Allocate deltas:
  for(L = 0; L < Nlayers; L++){  delta.push_back(MATRIX(Npe[L], sz));  }


  //====== Back-propagate deltas =============
  /// L = Nlayers-1  
  delta[Nlayers-1] = target - Y[Nlayers-1];

  /// L = Nlayers - 2, Nlayers-3, ..., 1
  for(L = Nlayers-2; L > 0; L--){

    delta[L] = W[L+1].T() * delta[L+1];

    /// This is the effect of the transfer function derivative
    for(i=0; i<Npe[L]; i++){ 
      for(j=0; j<sz; j++){ 
        delta[L].scale(i, j,  (1.0 - Y[L].get(i, j) * Y[L].get(i, j))  );
      }
    }
    
  }// for L      

 
  //========= Compute weight and bias gradients =======
  for(L = 1; L < Nlayers; L++){
    dW[L] = 0.0;
    dB[L] = 0.0;

    for(j = 0; j<sz; j++){
      dW[L] +=  delta[L].col(j) * Y[L-1].col(j).T();
      dB[L] +=  delta[L].col(j);
    }// for j

    dW[L] /= double(sz);
    dB[L] /= double(sz);

  }// for L


  // Compute error
  double err = 0.0;
  for(i=0; i<Npe[Nlayers-1]; i++){ 
    for(j=0; j<sz; j++){ 
      err += delta[Nlayers-1].get(i,j) * delta[Nlayers-1].get(i,j);
    }// for j
  }// for i
  err *= (0.5/double(sz));

  return err;

}




double NeuralNetwork::error(MATRIX& inputs, MATRIX& target){
/**
  This is a function to compute the error of the prediction

  Args:
    inputs - Npe[0] x num_patterns matrix of the inputs
    target - Npe[Nlayers-1] x num_patterns - the expected output of the ANN

  Returns: 
    the average error 

*/

  int i, j, L;  

  if(target.n_cols!=inputs.n_cols){    
    std::cout<<"Error: The number of patterns is different for inputs "<<inputs.n_cols<<" and targets "<<target.n_cols<<std::endl;
    exit(0);
  }

  int sz = target.n_cols; // number of patterns to handle at the same time

  /**

  L       0                   1                       ....             NL = Nlayers - 1

  W      [junk]              W[1]                                        W[NL]

  B      [junk]              B[1]                                        B[NL]

  Y      [Y[0]=input]      [ f(W[1]*Y[0] + B[1]) ]           [ output = f(W[NL]*Y[NL-1] + B[NL]) ]

 deltas  [junk]           W^T[2]*delta[2] *f'(Y[1])                   target - output[NL]
  
  */

  vector<MATRIX> Y;
  Y = propagate(inputs);


  MATRIX delta(Npe[Nlayers-1], sz);

  /// L = Nlayers-1  
  delta = target - Y[Nlayers-1];

  // Compute error
  double err = 0.0;
  for(i=0; i<Npe[Nlayers-1]; i++){ 
    for(j=0; j<sz; j++){ 
      err += delta.get(i,j) * delta.get(i,j);
    }// for j
  }// for i
  err *= (0.5/double(sz));

  return err;

}




void NeuralNetwork::init_weights_biases_normal(Random& rnd, double scaling_w, double shift_w, double scaling_b, double shift_b){

  for(int L =1; L < Nlayers; L++){
    for(int i=0; i<Npe[L]; i++){
      ///========== Weights ============
      for(int j=0; j<Npe[L-1]; j++){
        W[L].set(i, j, scaling_w * rnd.normal() + shift_w);
      }// for j
      ///========== Biases ============
      B[L].set(i, 0, scaling_b * rnd.normal() + shift_b );
    }// for i
  }// for L

}

void NeuralNetwork::init_weights_biases_uniform(Random& rnd, double left_w, double right_w, double left_b, double right_b){

  for(int L =1; L < Nlayers; L++){
    for(int i=0; i<Npe[L]; i++){
      ///========== Weights ============
      for(int j=0; j<Npe[L-1]; j++){
        W[L].set(i, j, rnd.uniform(left_w, right_w));
      }// for j
      ///========== Biases ============
      B[L].set(i, 0, rnd.uniform(left_b, right_b) );
    }// for i
  }// for L

}


vector<double> NeuralNetwork::train(Random& rnd, bp::dict params, MATRIX& inputs, MATRIX& targets){
/**
  See more details here:
  http://page.mi.fu-berlin.de/rojas/neural/chapter/K8.pdf
*/

  int i, j, epoch, L;

  ///============ Get the parameters ==================
  learning_rate = 0.0;
  momentum_term = 0.0;
  int num_epochs = 1;
  int steps_per_epoch = 1;
  int epoch_size = 1;
  int n_patterns = inputs.n_cols;
  int verbosity = 0;
  int is_error_collect_frequency = 0; 
  int error_collect_frequency = 1;
  vector<double> weight_decay(1, 0.0);   int is_weight_decay = 0;
  vector<double> bias_decay(1, 0.0);     int is_bias_decay = 0;

  std::string key;
  for(int i=0;i<len(params.values());i++){
    key = bp::extract<std::string>(params.keys()[i]);


    if(key=="learning_rate") { learning_rate = bp::extract<double>(params.values()[i]); }
    else if(key=="momentum_term") { momentum_term = bp::extract<double>(params.values()[i]); }
    else if(key=="num_epochs") { num_epochs = bp::extract<int>(params.values()[i]);   }
    else if(key=="steps_per_epoch") { steps_per_epoch = bp::extract<int>(params.values()[i]);   }
    else if(key=="epoch_size") { epoch_size = bp::extract<int>(params.values()[i]);   }
    else if(key=="verbosity") { verbosity = bp::extract<int>(params.values()[i]);   }
    else if(key=="error_collect_frequency") { 
      is_error_collect_frequency = 1;
      error_collect_frequency = bp::extract<int>(params.values()[i]);  
    }
    else if(key=="weight_decay"){
      is_weight_decay = 1;
      boost::python::list tmp = extract<boost::python::list>(params.values()[i]);      
      for(j=0; j<len(tmp); j++){  weight_decay.push_back( extract<double>(tmp[j]) );  }

      if(weight_decay.size()!=Nlayers){
        std::cout<<"Error: The number of weight decay constants should be equal to the number of weight matrices \n";
        std::cout<<"Now exiting...\n";
        exit(0);
      }// if      
    }
    else if(key=="bias_decay"){
      is_bias_decay = 1;
      boost::python::list tmp = extract<boost::python::list>(params.values()[i]);      
      for(j=0; j<len(tmp); j++){  bias_decay.push_back( extract<double>(tmp[j]) );  }

      if(bias_decay.size()!=Nlayers){
        std::cout<<"Error: The number of bias decay constants should be equal to the number of bias matrices \n";
        std::cout<<"Now exiting...\n";
        exit(0);
      }// if      
    }

  }


  if(!is_error_collect_frequency){ error_collect_frequency = steps_per_epoch;  }
  if(!is_weight_decay){ for(j=0; j<Nlayers-1; j++){  weight_decay.push_back(0.0); } }
  if(!is_bias_decay){ for(j=0; j<Nlayers-1; j++){  bias_decay.push_back(0.0); } }

  MATRIX input_subset(sz_x, epoch_size);
  MATRIX target_subset(sz_y, epoch_size);
  vector<int> subset(epoch_size);
  vector<int> inp_dim(sz_x); for(i=0; i<sz_x; i++){ inp_dim[i] = i; }
  vector<int> tar_dim(sz_y); for(i=0; i<sz_y; i++){ tar_dim[i] = i; }
  vector<MATRIX> Y;

  int counter = 0;
  double err_loc = 0.0;
  vector<double> err;
   
  for(epoch = 0; epoch < num_epochs; epoch++){    

    for(i = 0; i < steps_per_epoch; i++){
            
        // Make a random selection of the training patterns
        randperm(epoch_size, n_patterns, subset);
        
        // Extract the corresponding matrices from the inputs and outputs
        pop_submatrix(inputs,  input_subset,  inp_dim, subset);
        pop_submatrix(targets, target_subset, tar_dim, subset);

        // Update gradients and outputs
        Y = propagate(input_subset);
        err_loc = back_propagate( Y, target_subset);

        if(counter % error_collect_frequency ==0){  err.push_back( err_loc);       }
       

        // Update weights and biases
        // According to: http://page.mi.fu-berlin.de/rojas/neural/chapter/K8.pdf
        for(L = 0; L < Nlayers; L++){         

          dWold[L] = (learning_rate * (dW[L] - weight_decay[L]*W[L]) + momentum_term * dWold[L]);
          dBold[L] = (learning_rate * (dB[L] - bias_decay[L]*B[L]) + momentum_term * dBold[L]);

          W[L] += dWold[L];
          B[L] += dBold[L];

        }


        counter++;

    }// for i

    if(verbosity>=1){  cout<<"epoch = "<<epoch<<" (local) error = "<<err_loc<<"\n";  }

  }// for epoch

  return err;

}




/*

int NeuralNetwork::Propagate(boost::python::list input, boost::python::list& result){

  int i;
  if(len(input)!=sz_x){
    std::cout<<"Error: Size of the input "<<len(input)<<" does not match the ANN architecture "<<sz_x<<std::endl;
  }

  MATRIX* x;  x = new MATRIX(sz_x, 1);
  MATRIX* y;  y = new MATRIX(sz_y, 1);

  for(i=0; i<sz_x; i++){  x->M[i] = extract<double>(input[i]);   }

  this->Propagate(*x, *y);

  for(i=0; i<sz_y; i++){  result.append(y->M[i]); }


  delete x;
  delete y;

}

*/


int NeuralNetwork::Propagate(boost::python::list input, boost::python::list& result, boost::python::list& derivatives){

  int i;
  if(len(input)!=sz_x){
    std::cout<<"Error: Size of the input "<<len(input)<<" does not match the ANN architecture "<<sz_x<<std::endl;
  }

  MATRIX* x;  x = new MATRIX(sz_x, 1);
  MATRIX* y;  y = new MATRIX(sz_y, 1);
  MATRIX* d;  d = new MATRIX(sz_d, 1);

  for(i=0; i<sz_x; i++){  x->M[i] = extract<double>(input[i]);  }

  this->Propagate(*x, *y, *d);

  for(i=0; i<sz_y; i++){  result.append(y->M[i]); }
  for(i=0; i<sz_d; i++){  derivatives.append(d->M[i]); }


  delete x;
  delete y;
  delete d;

}




int NeuralNetwork::Propagate(const MATRIX& input, MATRIX& result, MATRIX& derivs){
/**  This function propagates a given input through the ANN and also computes 
   the derivaties of the outputs w.r.t. inputs

   input - input matrix
   result - the output of the ANN
   derivs - the derivatives

*/

//  sz_x = Inputs.size();
//  sz_y = Outputs.size();

  // Input check
  if(input.n_rows!=sz_x){
    std::cout<<"Error: Size of the input "<<input.n_rows
    <<" does not match the ANN architecture "<<Npe[0]<<std::endl;
    exit(0);
  }

  // Define and setup some variables
  int NL = Nlayers - 1;
  int i,L;
  double tmp;
  MATRIX x(sz_x,1);
  vector<MATRIX> Y;
  vector<MATRIX> gamma,ksi; // "Conjugate" variables, ksi[L] - is basically dx[L]/dx[0]

  // Setup some variables
  result = 0.0;

  MATRIX g(sz_x,sz_x); g.Init_Unit_Matrix(1.0);
  gamma.push_back(g);
  ksi.push_back(g); // thus ksi[0] = I - unity matrix of size sz_x - by - sz_x

  for(L=1;L<=NL;L++){
      MATRIX d1(Npe[L],sz_x); d1 = 0.0;

      gamma.push_back(d1);
      ksi.push_back(d1);
  }// for L


  for(i=0;i<sz_x;i++){

    //---------- Use the same linear transformation of input as during the training ----------
    tmp=input.M[i];

    tmp = Inputs[i].scale_factor * tmp + Inputs[i].shift_amount;
    x.M[i] = tmp;

  }// for i



  //-----------------------------------------------------
  //------- Forward propagation of the signal ------------
  vector<MATRIX> tmpF;
  Y.clear();
  Y.push_back(x); //0-th item;

  for (L=1;L<=NL;L++){

    MATRIX NET(Npe[L],1);
    MATRIX   y(Npe[L],1);

    NET = W[L] * Y[L-1] + B[L];
    D[L] = 0.0;


    for(int j=0;j<Npe[L];j++){

      y.M[j]=tanh(NET.M[j]);
      D[L].M[j*Npe[L]+j] = (1.0 - y.M[j]*y.M[j]);

    }
          /*
          MATRIX F(Npe[L],Npe[L-1]);

          F = D[L]*W[L];

          if(tmpF.size()==0){
          tmpF.push_back(F);
          }else{
            MATRIX tmp_F(Npe[L],tmpF[0].num_of_cols);
            tmp_F = F*tmpF[0];
            tmpF.clear();
            tmpF.push_back(tmp_F);
          }
          */
    Y.push_back(y);
 
    gamma[L] = W[L]*ksi[L-1];
    ksi[L]   = D[L]*gamma[L];

  }// for L

  derivs = ksi[NL]; // tmpF[0];
  //    tmpF.clear();


  //------- Linear transformation of the output ----------------

  for(i=0;i<Npe[NL];i++){
    tmp = Y[NL].M[i];

    // Inverse transform of output
    if(scale_method=="normalize_and_transform"){
      // Invert non-linear transformation
      //   tmp = tanh(tmp) = (exp(tmp)-exp(-tmp))/(exp(tmp)+exp(-tmp))   =>

      if(tmp>=0.99){ tmp = 0.99; }
      else if(tmp<=-0.99) { tmp = -0.99; }

      tmp = 0.5*log((1.0+tmp)/(1.0-tmp)); // <=
    }

    // Invert linear transformation of the output
    if(Outputs[i].scale_factor!=0.0){
         tmp = (1.0/Outputs[i].scale_factor)*(tmp - Outputs[i].shift_amount);
    }

    result.M[i] = tmp;

  }// for i

  //--------- Linear transform of the derivatives ----------------

  for(i=0;i<sz_y;i++){
    for(int j=0;j<sz_x;j++){
      // scale dYi/dXj
      derivs.M[sz_x*i+j] = (1.0/Derivs[sz_x*i+j].scale_factor)*(derivs.M[sz_x*i+j]);
    }// for j
  }// for i

  return 0;
}




void NeuralNetwork::ANNTrain(){
/**
  This function will train the ANN using the parameters set up and the methods 
  selected.
*/

  //=====================================================================
  //========== Parameters and auxiliary valiables are here ==============
  //=====================================================================



  //----------- Parameters --------------------------

  int NL= Nlayers-1;// Maximal index of matrixes and biases

  sz_x = Inputs.size();
  sz_y = Outputs.size();
  num_of_patterns = Inputs[0].Data.size();

  int nrows=sz_x; 
  int ncols=num_of_patterns;   

  int nrowsd=sz_y; 
  int ncolsd=num_of_patterns;  

  double dmax = 50;
  double dmin = 1e-8;   
  

  std::cout<<" sz_x (input space dimension)  = "<<sz_x<<"\n";
  std::cout<<" sz_y (output space dimension) = "<<sz_y<<"\n";
  std::cout<<" sz_d (derivatives space dimension) = "<<sz_d<<"\n";


  //----------- Auxiliary variables ----------------
  int i,j,L;

  vector<int> rperm;   // a random permutation

  MATRIX X(sz_x,1);
  MATRIX Target(sz_y,1);
  MATRIX Tangent(sz_d,1);
  MATRIX derivs(sz_y,sz_x);

  vector<MATRIX> D2; // Second derivatives
  vector<MATRIX> Tau;   // Delta(already have) and Tau 
  vector<MATRIX> gamma,ksi; // "Conjugate" variables, ksi[L] - is basically dx[L]/dx[0]
  MATRIX interm1(Npe[NL],Npe[NL]);

  MATRIX d(3,1); d = 0.0;
  D2.push_back(d);
  Tau.push_back(d);

  MATRIX g(sz_x,sz_x); g.Init_Unit_Matrix(1.0);
  gamma.push_back(g);
  ksi.push_back(g); // thus ksi[0] = I - unity matrix of size sz_x - by - sz_x

  for(L=1;L<=NL;L++){
    MATRIX d1(Npe[L],sz_x); d1 = 0.0;
    MATRIX d2(Npe[L],Npe[L]); d2 = 0.0;
 
    Tau.push_back(d1);
    D2.push_back(d2);
 
    gamma.push_back(d1);
    ksi.push_back(d1);
  }// for L  



  //===============================================================================
  //============== Checking for appropriate settings is here ======================
  //===============================================================================


  //------------ First stage: Need to check all settings and to set up them if necessary --------------
  
  if(!is_learning_method){
    std::cout<<"Error: Learning method is not defined!\n";
    std::cout<<"Now exiting...\n";
    exit(103);
  }else{ // Here is an initialization for different methods

    if(learning_method=="RProp"){      

         epoch_size = num_of_patterns;
         is_epoch_size = 1;

    }// if RProp


  }// else
  
  if(!is_learning_rate){ 
    learning_rate = 0.01/((double)num_of_patterns);
    is_learning_rate = 1;
	std::cout<<"Warning: Learning rate is not defined!\n";
	std::cout<<"Setting learning rate to = "<<learning_rate<<"\n";
  }
  if(!is_epoch_size){
     std::cout<<"Warning: The epoch_size parameter has not been defined\n";
     std::cout<<"Using default value = 1 (online training)\n";
         epoch_size = 1;
      is_epoch_size = 1; 
  }
  if(!is_momentum_term){
     std::cout<<"Warning: The momentum term will not be used\n";
     std::cout<<"Using default value of momentum_term parameter = 0\n";
         momentum_term = 0.0;
      is_momentum_term = 1;
  }
  if(!is_grad_weight){
     std::cout<<"Warning: The gradients of the ANN will not be fitted\n";
     std::cout<<"Using default value of grad_weight parameter = 0\n";
         grad_weight = 0.0;
      is_grad_weight = 1;
  }else{
     if(!derivs_flag){
        std::cout<<"Error: Can not use information about derivatives:  the derivatives are not defined\n";
        std::cout<<"Now exiting...\n";
        exit(103);
     }
  }
  if(!is_weight_decay){
     std::cout<<"Warning: The weight decay constants has not been defined\n";
     std::cout<<"Setting them all to default value = 0\n";
     for(i=0;i<NL;i++){   weight_decay.push_back(0.0);}
     is_weight_decay = 1;
  }else{
        if(weight_decay.size()!=NL){
           std::cout<<"Error: The number of decay constants should be equal to the number of weight matrices\n";
           std::cout<<"Now exiting...\n";
           exit(91);
        }// if
  }// else
  if(!is_norm_exp){
     std::cout<<"Warning: The norm exponent has not been defined\n";
     std::cout<<"Using default value = 0\n";
     std::cout<<"This corresponds to error as sum{ (t-o)^2 }\n";
         norm_exp = 0; // This corresponds to error as  sum{(t-o)^2}
      is_norm_exp = 1;
  }
  if(!is_iterations_in_cycle){
     std::cout<<"Warning: The number of iterations in one cycle has not been defined\n";
     std::cout<<"Using default value = 1000\n";
     iterations_in_cycle = 1000;
     is_iterations_in_cycle = 1;
  }
  if(!is_a_plus){
     std::cout<<"Warning: The a_plus parameter has not been defined\n";
     std::cout<<"Using default value = 1.2\n";
     a_plus = 1.2;
     is_a_plus = 1;
  }
  if(!is_a_minus){
     std::cout<<"Warning: The a_minus parameter has not been defined\n";
     std::cout<<"Using default value = 0.5\n";
     a_minus = 0.5;
     is_a_minus = 1;
  }  

  //============================================================================
  //============== Training algorithms are here ================================
  //============================================================================

  // Initialize random number generator
//  srand((unsigned)time(0)); 

  while (Iteration<iterations_in_cycle){

    // Clear container
    if(rperm.size()>0){ rperm.clear(); }

    // Choose a random sub-set of the training data set 
    // after this operation, the variable rperm will contain epoch_size integers
    // ranging from 0 to num_of_patterns
    // Essentially, the rperm will contain a subset of the training patterns used 
    // in this eapoch to train the ANN
    randperm(epoch_size,num_of_patterns,rperm);

    // Initialize total gradients (summ over epoch size)
    for (L=1;L<=NL;L++){
      dWcurr[L] = 0.0;
      dBcurr[L] = 0.0;
    }

    // Iterate over all subset patterns
    for(int ep=0;ep<epoch_size;ep++){

      // Pick a single pattern - inputs, outputs, and may be derivatives
      int indx = rperm[ep];
//      cout<<"Iteration = "<<Iteration<<" ep = "<<ep<<" indx = "<<indx<<endl;
      
      for(j=0;j<sz_x;j++){ X.M[j]       = Inputs[j].Data[indx];  }
      for(j=0;j<sz_y;j++){ Target.M[j]  = Outputs[j].Data[indx]; }

      if(derivs_flag){
        for(j=0;j<sz_d;j++){ Tangent.M[j] = Derivs[j].Data[indx];  }
      }else{
        for(j=0;j<sz_d;j++){ Tangent.M[j] = 0.0;}
      }

      //------------------------------------
      // Propagate forward - compute the inputs to all intermediate layers
      // and store them in the vector Y, layer by layer:
      // Y[0] - input layer, Y[1] - input to the first hidden layer, etc.
      // We will also compute the derivatives of these inputs w.r.t. to 
      // previous inputs (actually, to the previous NETs)
      vector<MATRIX> Y;  

      // 0-th item == input	
      Y.push_back(X);
      //Y.push_back(X);

      for(L=1;L<=NL;L++){

        MATRIX NET(Npe[L],1);
        MATRIX   y(Npe[L],1);

        NET = W[L]*Y[L-1] + B[L];
    
        /**
            (u/v)' =  (u'v - v'u)/v^2

          Y[j] = tanh( NET[j] )    tanh(x) = ( exp(x) - exp(-x) ) / (exp(x) + exp(-x) )

          dY[j]/dNET[j] =  (e(x) + e(-x) )^2 - (e(x) - e(-x) )^2 / (exp(x) + exp(-x) )^2 = 
 
         =  1 - tanh(x)^2

        */

        D[L] = 0.0;   // first derivatives of the transfer function ( tanh(net) ) 
        D2[L]= 0.0;   // second derivarives
       
        for(j=0;j<Npe[L];j++){

	  y.M[j]=tanh(NET.M[j]);  
          D[L].M[j*Npe[L]+j] = (1.0 - y.M[j]*y.M[j]); 
          D2[L].M[j*Npe[L]+j] = -2.0*y.M[j]*D[L].M[j*Npe[L]+j]; 

        }// for j
        
        Y.push_back(y);
          
        gamma[L] = W[L]*ksi[L-1];
        ksi[L]   = D[L]*gamma[L];
              
    }// for L
    
    // Calculate the error - at the output layer
    MATRIX e(Npe[NL],1);      e = 0.0;      // error in value
    MATRIX tau(Npe[NL],sz_x); tau = 0.0;    // error in the derivatives


    // Error for last (output) layer
    for (j=0;j<Npe[NL];j++){ 
      e.M[j]=pow((Target.M[j]-Y[NL].M[j]),(2.0*norm_exp+1));

      // Also, take the derivatives into account
      if(derivs_flag){
        for(i=0;i<sz_x;i++){
          tau.M[j*sz_x+i] = grad_weight*pow((Tangent.M[j*sz_x+i]-ksi[NL].M[j*sz_x+i]),(2.0*norm_exp+1));
        }//for i
      }              
    }// for j


        
    Delta[NL] = D[NL] * e;
    Tau[NL]   = D[NL] * tau;         
    interm1 = (tau * (gamma[NL].T()));
// NO GRADIENTS
//    for(i=0;i<Npe[NL];i++){Delta[NL].M[i] +=  interm1.M[i*Npe[NL]+i]*D2[NL].M[Npe[NL]*i+i]; }


    // Calculate deltas (derivatives)
    for(L=NL-1;L>=1;L--){

      MATRIX dEdx(Npe[L],1);
      MATRIX dEdksi(Npe[L],sz_x);
      MATRIX interm(Npe[L],Npe[L]);

      dEdx = (W[L+1].T() * Delta[L+1]);
      dEdksi = (W[L+1].T() * Tau[L+1]);

      Tau[L]   = D[L] * dEdksi;
      Delta[L] = D[L] * dEdx;

      interm = (dEdksi * (gamma[L].T()));

// NO GRADIENTS
//      for(i=0;i<Npe[L];i++){Delta[L].M[i] +=  interm.M[i*Npe[L]+i]*D2[L].M[Npe[L]*i+i] ; }

    }// for L

    // Calculate total gradient over all training patterns in epoch
    for (L=1;L<=NL;L++){

      // These are the negative gradients
//      dWcurr[L] += learning_rate*(Delta[L]*(Y[L-1].T()) + Tau[L]*(ksi[L-1].T()) - weight_decay[L]*W[L]);
// NO GRADIENTS
      dWcurr[L] += learning_rate*(Delta[L]*(Y[L-1].T()) - weight_decay[L]*W[L]);
      dBcurr[L] += learning_rate*Delta[L];

    }

  }// for ep - for all patterns in 1 epoch

  // Now add momentum term (if it is not zero)
  if(learning_method=="BackProp"){



    for(L=1;L<=NL;L++){
      /**   TESTING ON 4/3/2021!!!  : should be "+", but we try - */
      dW[L] = dWcurr[L] + momentum_term*dW[L];
      dB[L] = dBcurr[L] + momentum_term*dB[L];
    }

  }// if learning_method == BackProp

  // Or use QuickProp algorithm
  else if(learning_method=="QuickProp"){

      for (L=1;L<=NL;L++){

          if(Iteration==0){
              dW[L] = dWcurr[L];
              dB[L] = dBcurr[L];

          }else{

              double temp;
              for(int q1=0;q1<Npe[L];q1++){
                  for(int q2=0;q2<Npe[L-1];q2++){

                  //-------------- Weights -------------------
                      double St  = -dWcurr[L].M[q1*Npe[L-1]+q2];
                      double St1 = -dWold[L].M[q1*Npe[L-1]*q2];
                      double dw  = dW[L].M[q1*Npe[L-1]+q2];

                      temp = (St/(St1-St))*dw;

                     if(temp>dmax*dw){ temp = dmax*dw;  }
                     dW[L].M[q1*Npe[L-1]+q2] = temp;
                  //------------------------------------------

                  }// for q2
         
                  //-------------- Biases -------------------
                  double St  = -dBcurr[L].M[q1];
                  double St1 = -dBold[L].M[q1];
                  double dw  = dB[L].M[q1];

                  temp = (St/(St1-St))*dw;

                  if(temp>dmax*dw){ temp = dmax*dw;  }
                  dB[L].M[q1] = temp;

                  //-------------------------------------------

              }// for q1

          }// else: t>0

      }// for L

  }// if learning_method == QuickProp

  // or use RProp algorithm
  else if(learning_method=="RProp"){

      for(L=1;L<=NL;L++){

          if(Iteration==0){       
          // The very first initialization
              dWcurr[L] = 0.0;
              dBcurr[L] = 0.0;

              dW[L] = 0.1;
              dB[L] = 0.1;

          }else{

               // These are dE/dw and dE/db not -dE/dw and -dE/db, thus

              dWcurr[L] = -dWcurr[L]/learning_rate;
              dBcurr[L] = -dBcurr[L]/learning_rate;

              double pij;
              for(int q1=0;q1<Npe[L];q1++){       
                  for(int q2=0;q2<Npe[L-1];q2++){

                  //-------------- Weights -------------------
                      double Dercurr = dWcurr[L].M[q1*Npe[L-1]+q2];
                      double Derold  = dWold[L].M[q1*Npe[L-1]+q2];

                     if(Dercurr*Derold>0.0){
                         // Takes minimum
                         pij = (a_plus*dW[L].M[q1*Npe[L-1]+q2]<dmax)?(a_plus*dW[L].M[q1*Npe[L-1]+q2]):dmax;
                         dW[L].M[q1*Npe[L-1]+q2] = - SIGN(Dercurr)*pij;
                     }
                     else if(Dercurr*Derold<0){
                         // Takes maximum
                         pij = (a_minus*dW[L].M[q1*Npe[L-1]+q2]>dmin)?(a_minus*dW[L].M[q1*Npe[L-1]+q2]):dmin;
                         dW[L].M[q1*Npe[L-1]+q2] = - SIGN(Dercurr)*pij;
                     }
                     else{
                         pij = dW[L].M[q1*Npe[L-1]+q2];
                         dW[L].M[q1*Npe[L-1]+q2] = - SIGN(Dercurr)*pij;                        
                     }

                  //----------------------------------------

                  }// for q2           

                  //-------------- Biases -------------------
                  double Dercurr = dBcurr[L].M[q1];
                  double Derold  = dBold[L].M[q1];

                  if(Dercurr*Derold>0.0){
                      // Takes minimum
	             pij = (a_plus*dB[L].M[q1]<dmax)?(a_plus*dB[L].M[q1]):dmax;
                     dB[L].M[q1] = - SIGN(Dercurr)*pij;
                  }
                  else if(Dercurr*Derold<0){
                     // Takes maximum
                     pij = (a_minus*dB[L].M[q1]>dmin)?(a_minus*dB[L].M[q1]):dmin;
                     dB[L].M[q1] = - SIGN(Dercurr)*pij;
                  }
                  else{             
                     pij = dB[L].M[q1];
                     dB[L].M[q1] = - SIGN(Dercurr)*pij;                    
                  }
                  //----------------------------------------

              }// for q1
          }// else: t>0.0
      }// for L
  }// if learning_method == RProp

  // Or use conjugate gradients (Polak-Ribiere version)
  else if(learning_method=="ConjGradProp"){

      double betha = 0.0;
      double denom = 0.0;

      for(L=1;L<=NL;L++){
        
          for(int q1=0;q1<Npe[L];q1++){
              for(int q2=0;q2<Npe[L-1];q2++){

                  betha += dWcurr[L].M[q1*Npe[L-1]+q2]*(dWcurr[L].M[q1*Npe[L-1]+q2]-dWold[L].M[q1*Npe[L-1]+q2]);
                  denom += dWold[L].M[q1*Npe[L-1]+q2]*dWold[L].M[q1*Npe[L-1]+q2];
              }// for q2

              betha += dBcurr[L].M[q1]*(dBcurr[L].M[q1]-dBold[L].M[q1]);
              denom += dBold[L].M[q1]*dBold[L].M[q1];
          }// for q1
      }// for L

      if(fabs(denom)>1e-8){  betha = (betha/denom); }
      else{  betha = 0.0;}
   
      for (L=1;L<=NL;L++){

          dW[L] = dWcurr[L] + betha*dW[L];
          dB[L] = dBcurr[L] + betha*dB[L];

      }// for L
  }// if learning_method == ConjGradProp

  // Update weights and biases
  for ( L=1;L<=NL;L++){

      W[L] = W[L] + dW[L];  
      B[L] = B[L] + dB[L];

      dWold[L] = dWcurr[L];
      dBold[L] = dBcurr[L];

  }

  Iteration++;


  }// while

  // When done - set up iterations counter to 1, that is not to 0
  // this means if we repeat training cycle we will continue out
  // trainint rather then restart it from the beginning
  Iteration = 1;
  Cycle++;
 

}



void NeuralNetwork::LearningHistory(std::string filename,std::string data_flag){
// filename  - is a name of the file to which the information will be written
// data_flag - is a flag which take only 2 values: original, internal
//             it determines what the target values are - either original data
//             which we want to approximate or a transformed data on which 
//             the ANN is actually trained (scaled, internal representation)


//=====================================================================
//================= Here are some auxiliary variables =================
//=====================================================================

  int NL = Nlayers - 1;
  int i,j;
  sz_x = Inputs.size();
  sz_y = Outputs.size();
  sz_d = Derivs.size();

  ofstream output(filename.c_str(),ios::app);

  double error = 0.0;
  vector<double> Ave_Perc(Npe[NL]*(1+sz_x),0.0);      // Average relative error
  vector<double> Max_Perc(Npe[NL]*(1+sz_x),0.0);      // Maximal relative error
  vector<double> Min_Perc(Npe[NL]*(1+sz_x),1.0);      // Minimal relative error
  vector<double> Perc_less_half(Npe[NL]*(1+sz_x),0.0);// Percent of the pattern which relative error is less then 0.5
  vector<double> Error(Npe[NL]*(1+sz_x),0.0);         // Average error for each output (over all training patterns)

  MATRIX X(sz_x,1);
  MATRIX Target(sz_y,1);
  MATRIX Y(sz_y,1);
  MATRIX Tangent(sz_d,1);
  MATRIX derivs(sz_y,sz_x);
  MATRIX diff(sz_y,1);
                
  double perc;
  double max_perc,min_perc;
  double set_size = Inputs[0].Data.size();



  //==============================================================================
  //================ Calculate the errors and other parameters ===================
  //==============================================================================

  for(i=0;i<set_size;i++){

      if(data_flag=="internal"){ // For training testing

          // Transform to original inputs
          for(j=0;j<sz_x;j++){      X.M[j]=(1.0/Inputs[j].scale_factor)*(Inputs[j].Data[i]-Inputs[j].shift_amount);        }
          // Transform to original targets
          for(j=0;j<sz_y;j++){ Target.M[j]=(1.0/Outputs[j].scale_factor)*(Outputs[j].Data[i]-Outputs[j].shift_amount);     } 
          // Transform to original derivatives
          if(derivs_flag){
             for(j=0;j<sz_d;j++){ Tangent.M[j]=(1.0/Derivs[j].scale_factor)*(Derivs[j].Data[i]);     }
          }

      }// internal

      else if(data_flag=="original"){ // For recall testing

          for(j=0;j<sz_x;j++){      X.M[j]=Inputs[j].Data[i];        }
          // Transform to original targets
          for(j=0;j<sz_y;j++){ Target.M[j]=Outputs[j].Data[i];     }
          // Transform to original derivatives
          if(derivs_flag){
             for(j=0;j<sz_d;j++){ Tangent.M[j]=Derivs[j].Data[i];     }
          }

      }// original
      else{
          std::cout<<"Error: second argument may only take values: internal, original\n";
          std::cout<<"Now exiting...\n";
          exit(105);
      }

      this->Propagate(X,Y,derivs); // Y will be in external units

                                   					
      diff = (Target - Y);
  
                      
      // Functions
      for(j=0;j<Npe[NL];j++){
                          
          perc = fabs(diff.M[j])/fabs(Target.M[j]);

          if(perc<=0.5){ Perc_less_half[j] += (1.0/set_size); }
          Ave_Perc[j] += (perc/set_size);
          if(perc>=Max_Perc[j]){ Max_Perc[j] = perc;}
          if(perc<=Min_Perc[j]){ Min_Perc[j] = perc;}

          Error[j] += pow((diff.M[j]*diff.M[j]),(2*norm_exp+2));
      }// for j

      

      if(derivs_flag){
          // Their derivatives
          for(j=0;j<sz_d;j++){
              perc = fabs(derivs.M[j]-Tangent.M[j])/fabs(Tangent.M[j]);
              if(perc<=0.5){ Perc_less_half[sz_y+j] += (1.0/set_size); }
              Ave_Perc[sz_y+j] += (perc/set_size);
              if(perc>=Max_Perc[sz_y+j]){ Max_Perc[sz_y+j] = perc;}
              if(perc<=Min_Perc[sz_y+j]){ Min_Perc[sz_y+j] = perc;}

              Error[sz_y+j] += pow((derivs.M[j]-Tangent.M[j]),(2*norm_exp+2));

          }// for j

      }// derivs_flag



  }// for i - all training patterns




  //=====================================================================
  //============== Now just output parameters ===========================
  //=====================================================================

  output<<iterations_in_cycle*(Cycle)+(Iteration-1)<<"   ";
  for(i=0;i<Npe[NL];i++){
      output<<"[ ";
      output<<pow((Error[i]/set_size),(0.5/(norm_exp+1.0)))<<"  "<<Ave_Perc[i]<<"  "<<Max_Perc[i]<<"  "<<Min_Perc[i]<<"  "<<Perc_less_half[i]<<"  ";
      output<<" ]";
  }


  if(derivs_flag){
      for(i=0;i<sz_d;i++){
          output<<"[ ";
          output<<pow((Error[sz_y+i]/set_size),(0.5/(norm_exp+1.0)))<<"  "<<Ave_Perc[sz_y+i]<<"  "<<Max_Perc[sz_y+i]<<"  "<<Min_Perc[sz_y+i]<<"  "<<Perc_less_half[sz_y+i]<<"  ";
          output<<" ]";
      }// for i
  }// if derivs_flag


  output<<endl;


   
  output.close();


}


}// namespace libann
}// namespace liblibra

