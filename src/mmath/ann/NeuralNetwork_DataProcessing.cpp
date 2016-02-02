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

#include "NeuralNetwork.h"

using namespace boost;


/// libmmath namespace
namespace libmmath{

/// libann namespace
namespace libann{


int NeuralNetwork::CropTrainingData(boost::python::list minlim,boost::python::list maxlim){
/*  minlim and maxlim - are the lists of minimal and maximal inputs and outputs allowed
    dimesions of minlim and maxlin should coinside the dimensions of the patterns
*/
  int i,j;
  sz_x = Inputs.size();
  sz_y = Outputs.size();
  num_of_patterns = Inputs[0].Data.size();

 if(len(minlim)!=(sz_x+sz_y)){
    std::cout<<"The size of minlim list does not match the total number of input and output vectors in each pattern\n";
    std::cout<<"Exiting..\n";
    exit(103);
 }
 if(len(maxlim)!=(sz_x+sz_y)){
    std::cout<<"The size of maxlim list does not match the total number of input and output vectors in each pattern\n";
    std::cout<<"Exiting..\n";
    exit(103);
 }

 for(i=0;i<num_of_patterns;i++){
     int do_crop = 0;
     for(j=0;j<sz_x;j++){
         if(Inputs[j].Data[i]<minlim[j]) { do_crop = 1; break; }
         else if(Inputs[j].Data[i]>maxlim[j]) { do_crop = 1; break; }
     }
     if(!do_crop){
         for(j=0;j<sz_y;j++){
         if(Outputs[j].Data[i]<minlim[sz_x+j]) { do_crop = 1; break; }
         else if(Outputs[j].Data[i]>maxlim[sz_x+j]) { do_crop = 1; break; }
         }
     }// if !do_crop

     
     if(do_crop){
        for(j=0;j<sz_x;j++){  Inputs[j].Data.erase(Inputs[j].Data.begin()+i);     }
        for(j=0;j<sz_y;j++){ Outputs[j].Data.erase(Outputs[j].Data.begin()+i);     }
        i--;
        num_of_patterns--;
        
     }// if do_crop

 }// for i

  return 0;
}

int NeuralNetwork::ScaleDerivatives(){
// This is an auxiliary function which should be called 
// after all other scaling functions(which first scale inputs
// and outputs) and therefore contain corresponding scaling factors

// Scale derivatives:
  sz_x = Inputs.size();
  sz_y = Outputs.size();
  int i,j;

  if(derivs_flag){

  for(i=0;i<sz_y;i++){
      for(j=0;j<sz_x;j++){
          // scale dYi/dXj
          Derivs[sz_x*i+j].ScaleData(Outputs[i].scale_factor/Inputs[j].scale_factor);
      }// for j
  }// for i
  }// is sz_d = sz_x*sz_y

  return 0;

}

int NeuralNetwork::ScaleTrainingData(int X_scale,int Y_scale){

/*-------------------------------------
 This is one of the overloaded functions
 It scales the input (if X_scale ==1 ) and 
 the output (if Y_scale == 1) data onto the 
 interval [-0.95,0.95]
-----------------------------------------*/

  int i,j;
  double scale = 0.95;

  sz_x = Inputs.size();
  sz_y = Outputs.size();

  if(X_scale){
    for(i=0;i<sz_x;i++){  Inputs[i].ScaleData(-scale,scale);   }
  }
  if(Y_scale){
    for(i=0;i<sz_y;i++){  Outputs[i].ScaleData(-scale,scale);  }
  }

  this->ScaleDerivatives();
  scale_method = "scale_default";
 
  return 0;  

}


int NeuralNetwork::ScaleTrainingData(int X_scale,int Y_scale,boost::python::list minlim,boost::python::list maxlim){
/*  minlim and maxlim - are the lists of minimal and maximal inputs and outputs allowed
    in other words all input and output vectors will be mapped into corresponding intervals
    [minlim,maxlim]
    dimesions of minlim and maxlin should coinside the dimensions of the patterns

 This is one of the overloaded functions
 It scales the input and the output data onto the
 interval defined be minlim and maxlim lists

 minlim = [input minimal limits, output minimal limits]
 maxlim = [input maximal limits, output maximal limits]

-----------------------------------------*/

        sz_x = Inputs.size();
        sz_y = Outputs.size();     

        if(len(minlim)!=(sz_x+sz_y)){
           std::cout<<"The size of minlim list does not match the total number of input and output vectors in each pattern\n";
           std::cout<<"Exiting..\n";
           exit(103);
        }
        if(len(maxlim)!=(sz_x+sz_y)){
           std::cout<<"The size of maxlim list does not match the total number of input and output vectors in each pattern\n";
           std::cout<<"Exiting..\n";
           exit(103);
        }


        int i;

        if(X_scale){
        for(i=0;i<sz_x;i++){
            double min_val = extract<double>(minlim[i]);
            double max_val = extract<double>(maxlim[i]);

            Inputs[i].ScaleData(min_val,max_val); 

            std::cout<<i<<"-th inputs has been scaled by linear transformation: ";
            std::cout<<"X = "<<Inputs[i].scale_factor<<"X + "<<Inputs[i].shift_amount<<endl;
        }  
        }
        if(Y_scale){
        for(i=0;i<sz_y;i++){
            double min_val = extract<double>(minlim[sz_x+i]);
            double max_val = extract<double>(maxlim[sz_x+i]);

            Outputs[i].ScaleData(min_val,max_val);  
            std::cout<<i<<"-th outputs has been scaled by linear transformation: ";
            std::cout<<"X = "<<Outputs[i].scale_factor<<"X + "<<Outputs[i].shift_amount<<endl;

        }
        }

   this->ScaleDerivatives();


   scale_method = "scale";


   return 0;

}


int NeuralNetwork::NormalizeTrainingData(int X_scale,int Y_scale){

        sz_x = Inputs.size();
        sz_y = Outputs.size();

        int i;
        if(X_scale){
        for(i=0;i<sz_x;i++){  Inputs[i].NormalizeData();   }
        }
        if(Y_scale){
        for(i=0;i<sz_y;i++){  Outputs[i].NormalizeData();  }
        }

        this->ScaleDerivatives();

        scale_method = "normalize";

  return 0;
}


int NeuralNetwork::NormalizeAndScaleTrainingData(int X_scale,int Y_scale,boost::python::list minlim,boost::python::list maxlim){

        sz_x = Inputs.size();
        sz_y = Outputs.size();

        int i;
        if(X_scale){
        for(i=0;i<sz_x;i++){
            double min_val = extract<double>(minlim[i]);
            double max_val = extract<double>(maxlim[i]);

            Inputs[i].NormalizeData();      
            Inputs[i].ScaleData(min_val,max_val);           
        }
        }
        if(Y_scale){
        for(i=0;i<sz_y;i++){
            double min_val = extract<double>(minlim[sz_x+i]);
            double max_val = extract<double>(maxlim[sz_x+i]);

            Outputs[i].NormalizeData();
            Outputs[i].ScaleData(min_val,max_val);
        }
        }

        this->ScaleDerivatives();

        scale_method = "normalize_and_scale";


  return 0;
}


int NeuralNetwork::NormalizeAndTransformTrainingData(int X_scale,int Y_scale){

        sz_x = Inputs.size();
        sz_y = Outputs.size();
        num_of_patterns = Inputs[0].Data.size();

        int i,j;
        if(X_scale){
        for(i=0;i<sz_x;i++){  Inputs[i].NormalizeData();   }
        }
        if(Y_scale){
        for(i=0;i<sz_y;i++){  Outputs[i].NormalizeData();  }
        }

        //============ And now make a non-linear transformation of output values : x'' = tanh(x') ==============
        // thus mapping x' in [min,max] into x'' in [-1,1]
               
        for(j=0;j<sz_y;j++){ // Transform output to [-1,1]

           for(i=0;i<num_of_patterns;i++){

                Outputs[j].Data[i] = tanh(Outputs[j].Data[i]);


           }// for i        
        }// for j

        // Here the scaling of the derivatives is not so straightforward!


        scale_method = "normalize_and_transform";



  return 0;
}



int NeuralNetwork::SetTrainingData(object obj,int der_flag){

/*

  Need to set the training_set object as follows:       

  obj - is a list of patterns, where each pattern is the 
        list of two other lists - one of them - is input
        another list - is output

  obj = [ pattern1,
          pattern2,
            ...
          patternN
        ]

  pattern.Input  = [ vector of size n]
  pattern.Output = [ vector of size m]
  pattern.Derivs = [ vector of size n*k], may be omitted if der_flag = 0

  Derivatives should be ordered as follows:
  Derivs = [dy0/dx0, dy0/dx1, dy0/dx2, ..., dy0/dxn, 
            dy1/dx0, dy1/dx1, dy1/dx2, ..., dy1/dxn,
                       ...........
            dym/dx0, dym/dx1, dym/dx2, ..., dym/dxn
           ]

  This is for the ANN architecture:

 Input                              Output
 layer                              layer

    1                                  1
    2                                  2
   ...   --->   Hidden layers --->    ...

    n                                  m

    der_flag - do we input derivatives as well, or we use only function values
    der_flag = 0 - set up only function values as targets
    der_flag = 1 - set up also derivatives of the functions with respect to inputs as additional targets

*/


  std::cout<<"Setting data set...\n";
  int i;
  derivs_flag = der_flag;
  //--- Extract patterns from the object --------

  if(TrainData.size()>0) { TrainData.clear(); }

  for(i=0;i<len(obj);i++){

   ANNData pattern;
   int is_input = 0;
   int is_output = 0;
   int is_derivs = 0;

   boost::python::object obi = extract<boost::python::object>(obj[i]);

   set_list(is_input,  pattern.Input,  obi,"Input");
   set_list(is_output, pattern.Output, obi,"Output");

   if(der_flag){
      set_list(is_derivs, pattern.Derivs, obi,"Derivs");
      if(!is_derivs){
         std::cout<<"Error: Derivatives are not derined in pattern\n";
         std::cout<<"You either should define 'Derivs' keyword in your python object or use derivs_flag = 0\n";
         std::cout<<"Now exiting...\n";
         exit(101);
      }
   }   

   if(is_input&&is_output){
      TrainData.push_back(pattern);
   }

  }// for i

  

  //---- Now check if all input vectors have the same size ------
  int is_the_same = 1;

  num_of_patterns = TrainData.size();

  for(i=0;i<(num_of_patterns-1);i++){
     if(TrainData[i].Input.size()!=TrainData[i+1].Input.size()){
        is_the_same = 0;
        std::cout<<"Warning: The size of input vectors "<<i<<" and "<<(i+1)<<" is different"<<std::endl;
     }     
  }
  if(!is_the_same){ exit(102); }
  else{sz_x = TrainData[i].Input.size(); }

  //---- Now check if all output vectors have the same size ------
  is_the_same = 1;

  for(i=0;i<(num_of_patterns-1);i++){
     if(TrainData[i].Output.size()!=TrainData[i+1].Output.size()){
        is_the_same = 0;
        std::cout<<"Warning: The size of output vectors "<<i<<" and "<<(i+1)<<" is different"<<std::endl;
     }
  }
  if(!is_the_same){ exit(102); }
  else{sz_y = TrainData[i].Output.size(); }

  if(derivs_flag){
  //---- Now check if all derivative vectors have the same size ------
  // and also if the size of the derivatives vectors is equal to the product
  // of sizes of inputs and outputs 
  is_the_same = 1;

  for(i=0;i<(num_of_patterns-1);i++){
     if(TrainData[i].Derivs.size()!=TrainData[i+1].Derivs.size()){
        is_the_same = 0;
        std::cout<<"Warning: The size of output vectors "<<i<<" and "<<(i+1)<<" is different"<<std::endl;
     }   
  }
  if(!is_the_same){ exit(102); }
  else{sz_d = TrainData[i].Derivs.size(); }

  if(sz_d!=sz_x*sz_y){
     std::cout<<"Error: The number of derivatives should be the number of outputs times the number of inputs\n";
     std::cout<<"Actual values are:\n sz_x (inputs)      = "<<sz_x<<"\n";
     std::cout<<" sz_y (outputs)     = "<<sz_y<<"\n";
     std::cout<<" sz_d (derivatives) = "<<sz_d<<"\n";
     std::cout<<"Now exiting...\n";
     exit(102);
  }

  }// if derivs_flag





  // Now set up a new data type
  int j;

  // Set up inputs
  if(Inputs.size()>0){  Inputs.clear();  }

  for(i=0;i<sz_x;i++){
      vector<double> tmp;
     

      for(j=0;j<num_of_patterns;j++){
          tmp.push_back(TrainData[j].Input[i]);
      }
      DATA d(tmp); 
      Inputs.push_back(d);

  }


  // Set up outputs
  if(Outputs.size()>0){ Outputs.clear(); }

  for(i=0;i<sz_y;i++){
 
      vector<double> tmp1;
  
      for(j=0;j<num_of_patterns;j++){
          tmp1.push_back(TrainData[j].Output[i]);
      }
      DATA d1(tmp1);
      Outputs.push_back(d1);

  }

 
  // Set up derivatives of the outputs with respect to inputs

  if(Derivs.size()>0) { Derivs.clear(); }

// The order of these two for-cycles is important: we want to order
// derivatives as:
// Derivs[0] = dy0/dx0, Derivs[1] = dy0/dx1, ... Derivs[n] = dy0/dxn
// ans so on

  for(i=0;i<sz_y;i++){
      for(j=0;j<sz_x;j++){
          // scale dYi/dXj
          vector<double> tmp2;

          for(int p=0;p<num_of_patterns;p++){
              tmp2.push_back(TrainData[p].Derivs[sz_x*i+j]);
          }

          DATA d2(tmp2);
          Derivs.push_back(d2);

      }// for j
  }// for i



  return num_of_patterns;
}


}// namespace libann
}// namespace libmmath

