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

#include "NeuralNetwork.h"

/// liblibra namespace
namespace liblibra{

using namespace boost;

/// libann namespace
namespace libann{




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
    grad_w[L] = 0.0;
    grad_b[L] = 0.0;

    for(j = 0; j<sz; j++){
      grad_w[L] -=  delta[L].col(j) * Y[L-1].col(j).T();
      grad_b[L] -=  delta[L].col(j);
    }// for j

    grad_w[L] /= double(sz);
    grad_b[L] /= double(sz);

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



vector<double> NeuralNetwork::train(Random& rnd, bp::dict params, MATRIX& inputs, MATRIX& targets){
/**
  References:
 
  [1] http://page.mi.fu-berlin.de/rojas/neural/chapter/K8.pdf
  [2] https://arxiv.org/pdf/1711.05101.pdf
  [3] http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.17.1332

*/

  int i, j, epoch, L, a1, a2;

  ///============ Get the parameters ==================
  /**
    Learning class in parenthesis

    1  (1) - Back Propagation (BProp) and options, no momentum, Algorithm 1 of [2], neither purple nor green [default]
    11 (1) - BProp with L2 regularization - Algorithm 1 of [2], purple option
    12 (1) - BProp with decoulpled decay  - Algorithm 1 of [2], green option
    13 (1) - Adam with L2 regularization - Algorithm 2 of [2], purple option
    14 (1) - Adam with decoupled decay - Algorithm 2 of [2], green option
    2  (2) - Resilient Propagation without weight-backtracking (RProp-) - Eq. 1 of [3] + section 2.2
    21 (2) - Resilient Propagation with weight-backtracking (RProp+) - Eq. 1 of [3] + section 2.1
    22 (2) - Modified RProp- (iRprop-)
    23 (2) - Modified RProp+ (iRprop+)


  */
  int learning_method = 1; 
  int learning_class = 1;

  /// `alpha` in algorithm (1) of [2]
  double learning_rate = 0.001;      

  /// `beta_1` in algorithm (1) of [2]
  double momentum_term = 0.0;        

  /// `lambda` in algorithm (1) of [2]
  /// L2 regularization factor  L_new = L_old + lambda * w^T w
  double weight_decay_lambda = 0.0;  
                                     
  /// `etha_t` in algorithm (1) of [2]
  double etha = 1.0;

  int num_epochs = 1;
  int steps_per_epoch = 1;
  int epoch_size = 1;
  int n_patterns = inputs.n_cols;
  int verbosity = 0;
  int is_error_collect_frequency = 0; 
  int error_collect_frequency = 1;

  /// RProp parameters
  /// RProp parameter - good for up to 1.2
  double a_plus = 1.1;   

  /// RProp parameter - good for up to 0.5
  double a_minus = 0.6;  

  double dB_min =  0.1*learning_rate;
  double dB_max =  learning_rate;
  double dW_min =  0.1*learning_rate;
  double dW_max =  learning_rate;


  std::string key;
  for(int i=0;i<len(params.values());i++){
    key = bp::extract<std::string>(params.keys()[i]);


    if(key=="learning_rate") { learning_rate = bp::extract<double>(params.values()[i]); }
    else if(key=="learning_method") { learning_method = bp::extract<int>(params.values()[i]);   }
    else if(key=="momentum_term") { momentum_term = bp::extract<double>(params.values()[i]); }
    else if(key=="weight_decay_lambda") { weight_decay_lambda = bp::extract<double>(params.values()[i]); }
    else if(key=="etha") { etha = bp::extract<double>(params.values()[i]); }
 
    else if(key=="num_epochs") { num_epochs = bp::extract<int>(params.values()[i]);   }
    else if(key=="steps_per_epoch") { steps_per_epoch = bp::extract<int>(params.values()[i]);   }
    else if(key=="epoch_size") { epoch_size = bp::extract<int>(params.values()[i]);   }
    else if(key=="verbosity") { verbosity = bp::extract<int>(params.values()[i]);   }
    else if(key=="error_collect_frequency") { 
      is_error_collect_frequency = 1;
      error_collect_frequency = bp::extract<int>(params.values()[i]);  
    }
    else if(key=="a_plus") { a_plus = bp::extract<double>(params.values()[i]); }
    else if(key=="a_minus") { a_minus = bp::extract<double>(params.values()[i]); }
    else if(key=="dB_min") { dB_min = bp::extract<double>(params.values()[i]); }
    else if(key=="dB_max") { dB_max = bp::extract<double>(params.values()[i]); }
    else if(key=="dW_min") { dW_min = bp::extract<double>(params.values()[i]); }
    else if(key=="dW_max") { dW_max = bp::extract<double>(params.values()[i]); }



  } // for i


  // Sanity check and other setups
  if(learning_method==1 || learning_method==11 || learning_method==12 || learning_method==13 || learning_method==14 ){
    learning_class = 1;
  }
  else if(learning_method == 2 || learning_method == 21 || learning_method == 22 || learning_method == 23 ){
    learning_class = 2;
  }

  // Implementation status
  if(learning_method==13){  cout<<"The method 13 is not yet implemented: using 11 instead"; learning_method = 11;  }
  if(learning_method==14){  cout<<"The method 14 is not yet implemented: using 12 instead"; learning_method = 12;  }

  if(learning_method==21){  cout<<"The method 22 is not yet implemented: using 2 instead"; learning_method = 2;  }
  if(learning_method==22){  cout<<"The method 22 is not yet implemented: using 2 instead"; learning_method = 2;  }
  if(learning_method==23){  cout<<"The method 23 is not yet implemented: using 2 instead"; learning_method = 2;  }


  if(learning_class==1){ // Backprop setups

    for(L = 0; L < Nlayers; L++){ 
      dBold[L] = dB[L]; 
      dWold[L] = dW[L]; 
    }// for L

  }
  if(learning_class==2){   // RProp setups
    if(epoch_size!=n_patterns){  
      cout<<"WARNING in ANN.train : epoch_size ("<<epoch_size
          <<") should be equal to the total number of patters ("<<n_patterns<<")\n";
      epoch_size = n_patterns;
      cout<<"Using new value of epoch_size = "<<epoch_size<<endl;
    }

    // Initialize deltas: 
    for(L = 0; L < Nlayers; L++){ 

      for(a1=0; a1<Npe[L]; a1++){  
        dB[L].set(a1, 0, learning_rate * rnd.uniform(0.0, 1.0) );

        for(a2=0; a2<Npe[L-1]; a2++){    
          dW[L].set(a1, a2, learning_rate * rnd.uniform(0.0, 1.0) );
        }// for a2
      }// for a1
    }// for L

  }// RProp



  if(!is_error_collect_frequency){ error_collect_frequency = steps_per_epoch;  }

  MATRIX input_subset(sz_x, epoch_size);
  MATRIX target_subset(sz_y, epoch_size);
  vector<int> subset(epoch_size);
  vector<int> inp_dim(sz_x); for(i=0; i<sz_x; i++){ inp_dim[i] = i; }
  vector<int> tar_dim(sz_y); for(i=0; i<sz_y; i++){ tar_dim[i] = i; }
  vector<MATRIX> Y;

  int counter = 0;
  double err_loc = 0.0;
  vector<double> err;

  if(verbosity>0){
    cout<<"Training with parameters:\n";
    cout<<"learning_method = "<<learning_method<<endl;
    cout<<"learning_class = "<<learning_class<<endl;
    cout<<"learning_rate = "<<learning_rate<<endl;
    cout<<"momentum_term = "<<momentum_term<<endl;
    cout<<"weight_decay_lambda = "<<weight_decay_lambda<<endl;
    cout<<"etha = "<<etha<<endl;
    cout<<"num_epochs = "<<num_epochs<<endl;
    cout<<"steps_per_epoch = "<<steps_per_epoch<<endl;
    cout<<"epoch_size = "<<epoch_size<<endl;
    cout<<"n_patterns = "<<n_patterns<<endl;
    cout<<"verbosity = "<<verbosity<<endl;
    cout<<"is_error_collect_frequency = "<<is_error_collect_frequency<<endl;
    cout<<"error_collect_frequency = "<<error_collect_frequency<<endl;
    cout<<"a_plus = "<<a_plus<<endl;
    cout<<"a_minus = "<<a_minus<<endl;
    cout<<"dB_min = "<<dB_min<<endl;
    cout<<"dB_max = "<<dB_max<<endl;
    cout<<"dW_min = "<<dW_min<<endl;
    cout<<"dW_max = "<<dW_max<<endl;


  }// verbosity>0


  //===================================================================
   
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
        for(L = 0; L < Nlayers; L++){         

          //***************************** BProp family ********************************
          if(learning_class==1){ 

            //======================= Compute deltas dW and dB ========================
            // Momentum term is according to [1]
            // Wight decay regularization + momentum: [2] 
            // Here, dW and dB are understood as the momentum vectors
            dW[L] = momentum_term * dWold[L] + etha * learning_rate * grad_w[L];
            dB[L] = momentum_term * dBold[L] + etha * learning_rate * grad_b[L];

            if(learning_method==11 || learning_method==13){
              dW[L] += etha * learning_rate * weight_decay_lambda * W[L];
              dB[L] += etha * learning_rate * weight_decay_lambda * B[L];
            }

            //=========== Update step:  dW and dB are the momentum terms ==============
            W[L] -= dW[L];
            B[L] -= dB[L];

            if(learning_method==12 || learning_method==14){
              W[L] -=  etha * weight_decay_lambda * W[L];
              B[L] -=  etha * weight_decay_lambda * B[L];
            }

          }// BProp


          //***************************** RProp family ********************************
          // RProp   basic version: 
          else if(learning_class==2){

            // According to: [3], Eq. 1
            for(a1=0; a1<Npe[L]; a1++){
              // Biases:
              double gr_prod = grad_b[L].get(a1, 0) * grad_b_old[L].get(a1,0);  // product of B gradients

              // Scaling
              if(gr_prod > 0.0 ){   dB[L].scale(a1, 0,  a_plus);  } 
              else if(gr_prod < 0.0 ){ dB[L].scale(a1, 0, a_minus); }
              else{  ;; } // don't change


              //************ RProp- ****************
              if(learning_method==2){

                // Obey the bounds
                if( dB[L].get(a1, 0) > dB_max ){  dB[L].set(a1, 0, dB_max); }
                else if( dB[L].get(a1, 0) < dB_min ){  dB[L].set(a1, 0, dB_min); }

                //Update the biases
                B[L].add(a1, 0, -SIGN( grad_b[L].get(a1, 0) ) * dB[L].get(a1, 0) ); 

              }// method == 2

              //************ RProp+ ****************
              //if(learning_method==21){ ;; }



              // Weights 
              for(a2=0; a2<Npe[L-1]; a2++){

                gr_prod = grad_w[L].get(a1, a2) * grad_w_old[L].get(a1, a2);  // product of B gradients

                // Scaling
                if(gr_prod > 0.0 ){   dW[L].scale(a1, a2,  a_plus);  } 
                else if(gr_prod < 0.0 ){ dW[L].scale(a1, a2, a_minus); }
                else{  ;; } // don't change

                //************ RProp- ****************
                if(learning_method==2){

                  // Obey the bounds
                  if( dW[L].get(a1, a2) > dW_max ){  dW[L].set(a1, a2, dW_max); }
                  else if( dW[L].get(a1, a2) < dW_min ){  dW[L].set(a1, a2, dW_min); }

                  W[L].add(a1, a2, -SIGN( grad_w[L].get(a1, a2) ) * dW[L].get(a1, a2) ); 

                }// method == 2
                //************ RProp+ ****************
                //if(learning_method==21){ ;; }


              }// for a2
            }// for a1

          }// if RProp


          dWold[L] = dW[L];
          dBold[L] = dB[L];

          grad_w_old[L] = grad_w[L];
          grad_b_old[L] = grad_b[L];


        }// for L - layers


        counter++;

    }// for i

    if(verbosity>=1){  
      cout<<"epoch = "<<epoch<<" (local) error = "<<err_loc<<"\n"; 
    }
    if(verbosity>=2){
      for(L = 0; L < Nlayers; L++){         
        cout<<"dW["<<L<<"] = "; dW[L].show_matrix(); cout<<endl;
      }// L

      for(L = 0; L < Nlayers; L++){         
        cout<<"dB["<<L<<"] = "; dB[L].show_matrix(); cout<<endl;
      }// L

    }

  }// for epoch


  return err;

}



}// namespace libann
}// namespace liblibra

