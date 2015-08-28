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
//#include <ctime> 
//#include <cstdlib>

//using namespace std;
using namespace boost;
//using namespace libmmath;
/*
using libmmath::VECTOR;
using libmmath::MATRIX;
using libmmath::CMATRIX;
using libmmath::MATRIX3x3;
using libmmath::QUATERNION;
using libmmath::DATA;
*/

namespace libmmath{
namespace libann{


int NeuralNetwork::Propagate(boost::python::list input,boost::python::list& result){

 if(len(input)!=sz_x){
 std::cout<<"Error: Size of the input "<<len(input)<<" does not match the ANN architecture "<<sz_x<<std::endl;
 }else{

 int NL = Nlayers - 1;
 double tmp;

 vector<MATRIX> Y;
 MATRIX x(sz_x,1);


 for(int i=0;i<sz_x;i++){

//---------- Use the same linear transformation of input as during the training ----------
    tmp=extract<double>(input[i]);
    tmp = Inputs[i].scale_factor * tmp + Inputs[i].shift_amount;
    x.M[i] = tmp;

 }// for i
//-----------------------------------------------------

//------- Forward propagation of the signal ------------
     Y.clear();
     Y.push_back(x);//0-th item;

     for (int L=1;L<=NL;L++){
          MATRIX NET(Npe[L],1);
          MATRIX   y(Npe[L],1);

          NET=W[L]*Y[L-1];
          NET = NET + B[L];

          for(int j=0;j<Npe[L];j++){ y.M[j]=tanh(NET.M[j]);}
          Y.push_back(y);
    }


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

    result.append(tmp);

 }// for i

 }// else if sizes match

  return 0;
}


int NeuralNetwork::Propagate(MATRIX input,MATRIX& result){

// Clean result
 result = 0.0;

 if(input.num_of_rows!=sz_x){
 std::cout<<"Error: Size of the input "<<input.num_of_rows<<" does not match the ANN architecture "<<sz_x<<std::endl;
 }else{

 int NL = Nlayers - 1;
 double tmp;

 vector<MATRIX> Y;
 MATRIX x(sz_x,1);

 for(int i=0;i<sz_x;i++){

//---------- Use the same linear transformation of input as during the training ----------
    tmp=input[i];
    tmp = Inputs[i].scale_factor * tmp + Inputs[i].shift_amount;
    x.M[i] = tmp;

 }// for i


//-----------------------------------------------------
//------- Forward propagation of the signal ------------
     Y.clear();
     Y.push_back(x);//0-th item;

     for (int L=1;L<=NL;L++){
          MATRIX NET(Npe[L],1);
          MATRIX   y(Npe[L],1);

          NET=W[L]*Y[L-1];
          NET = NET + B[L];

          for(int j=0;j<Npe[L];j++){ y.M[j]=tanh(NET.M[j]);}
          Y.push_back(y);
    }


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



 }// else if sizes match

  return 0;
}

int NeuralNetwork::Propagate(boost::python::list input,boost::python::list& result,boost::python::list& derivatives){

// Clean result
 sz_x = Inputs.size();
 sz_y = Outputs.size();


 if(len(input)!=sz_x){
 std::cout<<"Error: Size of the input "<<len(input)<<" does not match the ANN architecture "<<sz_x<<std::endl;
 }else{

 int NL = Nlayers - 1;
 int L;
 double tmp;

 vector<MATRIX> Y;
 MATRIX x(sz_x,1);
 MATRIX derivs(sz_y,sz_x);

  vector<MATRIX> gamma,ksi; // "Conjugate" variables, ksi[L] - is basically dx[L]/dx[0]

  MATRIX g(sz_x,sz_x); g.Init_Unit_Matrix(1.0);
  gamma.push_back(g);
  ksi.push_back(g); // thus ksi[0] = I - unity matrix of size sz_x - by - sz_x

  for(L=1;L<=NL;L++){
      MATRIX d1(Npe[L],sz_x); d1 = 0.0;

      gamma.push_back(d1);
      ksi.push_back(d1);
  }// for L



 for(int i=0;i<sz_x;i++){

//---------- Use the same linear transformation of input as during the training ----------
    tmp=extract<double>(input[i]);
    tmp = Inputs[i].scale_factor * tmp + Inputs[i].shift_amount;
    x.M[i] = tmp;

 }// for i



//-----------------------------------------------------
//------- Forward propagation of the signal ------------
     vector<MATRIX> tmpF;
     Y.clear();
     Y.push_back(x);//0-th item;

     for (L=1;L<=NL;L++){
          MATRIX NET(Npe[L],1);
          MATRIX   y(Npe[L],1);

          NET=W[L]*Y[L-1];
          NET = NET + B[L];

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
          gamma[L] = W[L]*ksi[L-1];
          ksi[L]   = D[L]*gamma[L];

          Y.push_back(y);

    }
    derivs = ksi[NL]; //tmpF[0];
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

    result.append(tmp);

 }// for i

//--------- Linear transform of the derivatives ----------------

  for(i=0;i<sz_y;i++){
      for(int j=0;j<sz_x;j++){
          // scale dYi/dXj
          derivs.M[sz_x*i+j] = (1.0/Derivs[sz_x*i+j].scale_factor)*(derivs.M[sz_x*i+j]);
          derivatives.append(derivs.M[sz_x*i+j]);
      }// for j
  }// for i




 }// else if sizes match

  return 0;
}



int NeuralNetwork::Propagate(MATRIX input,MATRIX& result,MATRIX& derivs){

// Clean result
 result = 0.0;
 sz_x = Inputs.size();
 sz_y = Outputs.size();

 if(input.num_of_rows!=sz_x){
 std::cout<<"Error: Size of the input "<<input.num_of_rows<<" does not match the ANN architecture "<<sz_x<<std::endl;
 }else{

 int NL = Nlayers - 1;
 int L;
 double tmp;

 vector<MATRIX> Y;
 MATRIX x(sz_x,1);


  vector<MATRIX> gamma,ksi; // "Conjugate" variables, ksi[L] - is basically dx[L]/dx[0]

  MATRIX g(sz_x,sz_x); g.Init_Unit_Matrix(1.0);
  gamma.push_back(g);
  ksi.push_back(g); // thus ksi[0] = I - unity matrix of size sz_x - by - sz_x

  for(L=1;L<=NL;L++){
      MATRIX d1(Npe[L],sz_x); d1 = 0.0;

      gamma.push_back(d1);
      ksi.push_back(d1);
  }// for L





 for(int i=0;i<sz_x;i++){

//---------- Use the same linear transformation of input as during the training ----------
    tmp=input.M[i];

    tmp = Inputs[i].scale_factor * tmp + Inputs[i].shift_amount;
    x.M[i] = tmp;

 }// for i



//-----------------------------------------------------
//------- Forward propagation of the signal ------------
     vector<MATRIX> tmpF;
     Y.clear();
     Y.push_back(x);//0-th item;

     for (L=1;L<=NL;L++){
          MATRIX NET(Npe[L],1);
          MATRIX   y(Npe[L],1);

          NET=W[L]*Y[L-1];
          NET = NET + B[L];

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

    }
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




 }// else if sizes match

  return 0;
}




void NeuralNetwork::ANNTrain(){
  // This is a combination of different training methods
  // controlled by corresponding parameters

  //=====================================================================
  //========== Parameters and auxiliary valiables are here ==============
  //=====================================================================



  //----------- Patameters --------------------------

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

  vector<int> rperm;

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
  srand((unsigned)time(0)); 

  while (Iteration<iterations_in_cycle){

      // Clear container
      if(rperm.size()>0){ rperm.clear(); }

      // Choose a random sub-set of the training data set           
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
          for(j=0;j<sz_x;j++){ X.M[j]       = Inputs[j].Data[indx];  }
          for(j=0;j<sz_y;j++){ Target.M[j]  = Outputs[j].Data[indx]; }
          if(derivs_flag){
          for(j=0;j<sz_d;j++){ Tangent.M[j] = Derivs[j].Data[indx];  }
          }else{
          for(j=0;j<sz_d;j++){ Tangent.M[j] = 0.0;}
          }

          //------------------------------------
          // Propagate forward
          vector<MATRIX> Y;  
          Y.push_back(X);//0-th item; == is input	

          for (L=1;L<=NL;L++){

              MATRIX NET(Npe[L],1);
              MATRIX   y(Npe[L],1);
		
              NET=W[L]*Y[L-1];
              NET = NET + B[L];
              D[L] = 0.0;
              D2[L]= 0.0;
	
    	      for(j=0;j<Npe[L];j++){
	          y.M[j]=tanh(NET.M[j]);
 		  D[L].M[j*Npe[L]+j] = (1.0 - y.M[j]*y.M[j]);               
                  D2[L].M[j*Npe[L]+j] = -2.0*y.M[j]*D[L].M[j*Npe[L]+j];                   
              }// for j
        
              Y.push_back(y);
          
              gamma[L] = W[L]*ksi[L-1];
              ksi[L]   = D[L]*gamma[L];              
		
          }// for L
    
          // Calculate error
          MATRIX e(Npe[NL],1);      e = 0.0;
          MATRIX tau(Npe[NL],sz_x); tau = 0.0;

 	  // Error for last (output) layer
          for (j=0;j<Npe[NL];j++){ 
              e.M[j]=pow((Target.M[j]-Y[NL].M[j]),(2.0*norm_exp+1));

              if(derivs_flag){
              for(i=0;i<sz_x;i++){
                  tau.M[j*sz_x+i] = grad_weight*pow((Tangent.M[j*sz_x+i]-ksi[NL].M[j*sz_x+i]),(2.0*norm_exp+1));
              }//for i
              }              
          }// for j
        
	  Delta[NL] = D[NL] * e;
          Tau[NL]   = D[NL] * tau;         
          interm1 = (tau * (gamma[NL].T()));
          for(i=0;i<Npe[NL];i++){Delta[NL].M[i] +=  interm1.M[i*Npe[NL]+i]*D2[NL].M[Npe[NL]*i+i]; }


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

  	      for(i=0;i<Npe[L];i++){Delta[L].M[i] +=  interm.M[i*Npe[L]+i]*D2[L].M[Npe[L]*i+i] ; }

          }// for L

          // Calculate total gradient over all training patterns in epoch
          for (L=1;L<=NL;L++){

              // These are the negative gradients
              dWcurr[L] = dWcurr[L] + learning_rate*(Delta[L]*(Y[L-1].T()) + Tau[L]*(ksi[L-1].T()) - weight_decay[L]*W[L]);
              dBcurr[L] = dBcurr[L] + learning_rate*Delta[L];
          }
      }// for ep - for all patterns in 1 epoch

  // Now add momentum term (if it is not zero)
  if(learning_method=="BackProp"){

      for(L=1;L<=NL;L++){
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

               // These are dE/dw amd dE/db not -dE/dw and -dE/db, thus

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
}// namespace libmmath

