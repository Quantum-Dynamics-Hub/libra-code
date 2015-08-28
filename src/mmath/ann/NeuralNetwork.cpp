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

#include <ctime> 
#include <cstdlib>
#include <boost/regex.hpp>
#include "NeuralNetwork.h"

using namespace std;
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



ANNData& ANNData::operator=(const ANNData& pattern){

  Input  = pattern.Input;
  Output = pattern.Output;
  Derivs = pattern.Derivs;

  return *this;
}



void NeuralNetwork::set(object obj){

   set_value(is_learning_method,learning_method, obj, "learning_method");

   set_value(is_learning_rate, learning_rate, obj, "learning_rate");
   set_value(is_momentum_term, momentum_term, obj, "momentum_term");
   set_value(is_epoch_size,    epoch_size,    obj, "epoch_size");
   set_value(is_iterations_in_cycle,iterations_in_cycle,obj, "iterations_in_cycle");
   set_value(is_grad_weight,   grad_weight,   obj, "grad_weight");
   set_value(is_norm_exp,      norm_exp,      obj, "norm_exp");
   set_value(is_a_plus,        a_plus,        obj, "a_plus");
   set_value(is_a_minus,       a_minus,       obj, "a_minus");

   set_list(is_weight_decay,   weight_decay,  obj, "weight_decay");


}


NeuralNetwork::NeuralNetwork(const NeuralNetwork& ann){

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


  // Copy constructor

  if(ann.is_learning_method){ learning_method  = ann.learning_method; is_learning_method = 1; }  

  if(ann.is_learning_rate) { learning_rate  = ann.learning_rate; is_learning_rate = 1; }
  if(ann.is_momentum_term) { momentum_term  = ann.momentum_term; is_momentum_term = 1; }
  if(ann.is_epoch_size)    { epoch_size     = ann.epoch_size;    is_epoch_size    = 1; }
  if(ann.is_iterations_in_cycle){ iterations_in_cycle = ann.iterations_in_cycle;is_iterations_in_cycle= 1; }
  if(ann.is_grad_weight)   { grad_weight    = ann.grad_weight;   is_grad_weight   = 1; }
  if(ann.is_norm_exp)      { norm_exp       = ann.norm_exp;      is_norm_exp      = 1; }
  if(ann.is_a_plus)        { a_plus         = ann.a_plus;        is_a_plus        = 1; }
  if(ann.is_a_minus)       { a_minus        = ann.a_minus;       is_a_minus       = 1; }


  if(ann.is_weight_decay)  { weight_decay   = ann.weight_decay;  is_weight_decay  = 1; }

  Iteration = ann.Iteration;
  Cycle     = ann.Cycle;
 
  scale_method = ann.scale_method;
  derivs_flag  = ann.derivs_flag;

  num_of_patterns = ann.num_of_patterns;
  
  Recall       = ann.Recall;
  TrainData    = ann.TrainData;

  Inputs       = ann.Inputs;
  Outputs      = ann.Outputs;
  Derivs       = ann.Derivs;

  
  sz_x         = ann.sz_x;
  sz_y         = ann.sz_y;
  sz_d         = ann.sz_d;

  Nlayers      = ann.Nlayers;
  Npe          = ann.Npe;

  W            = ann.W;
  dW           = ann.dW;
  dWcurr       = ann.dWcurr;
  dWold        = ann.dWold;

  B            = ann.B;
  dB           = ann.dB;
  dBcurr       = ann.dBcurr;
  dBold        = ann.dBold;

  D            = ann.D;
  Delta        = ann.Delta;



}

NeuralNetwork::~NeuralNetwork(){
 ;; 
}


NeuralNetwork& NeuralNetwork::operator=(const NeuralNetwork& ann){

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


  if(ann.is_learning_method){ learning_method  = ann.learning_method; is_learning_method = 1; }  

  if(ann.is_learning_rate) { learning_rate  = ann.learning_rate; is_learning_rate = 1; }
  if(ann.is_momentum_term) { momentum_term  = ann.momentum_term; is_momentum_term = 1; }
  if(ann.is_epoch_size)    { epoch_size     = ann.epoch_size;    is_epoch_size    = 1; }
  if(ann.is_iterations_in_cycle){ iterations_in_cycle = ann.iterations_in_cycle;is_iterations_in_cycle= 1; }
  if(ann.is_grad_weight)   { grad_weight    = ann.grad_weight;   is_grad_weight   = 1; }
  if(ann.is_norm_exp)      { norm_exp       = ann.norm_exp;      is_norm_exp      = 1; }
  if(ann.is_a_plus)        { a_plus         = ann.a_plus;        is_a_plus        = 1; }
  if(ann.is_a_minus)       { a_minus        = ann.a_minus;       is_a_minus       = 1; }


  if(ann.is_weight_decay)  { weight_decay   = ann.weight_decay;  is_weight_decay  = 1; }

  Iteration = ann.Iteration;
  Cycle     = ann.Cycle;
 
  scale_method = ann.scale_method;
  derivs_flag  = ann.derivs_flag;

  num_of_patterns = ann.num_of_patterns;
  
  Recall       = ann.Recall;
  TrainData    = ann.TrainData;

  Inputs       = ann.Inputs;
  Outputs      = ann.Outputs;
  Derivs       = ann.Derivs;

  
  sz_x         = ann.sz_x;
  sz_y         = ann.sz_y;
  sz_d         = ann.sz_d;

  Nlayers      = ann.Nlayers;
  Npe          = ann.Npe;

  W            = ann.W;
  dW           = ann.dW;
  dWcurr       = ann.dWcurr;
  dWold        = ann.dWold;

  B            = ann.B;
  dB           = ann.dB;
  dBcurr       = ann.dBcurr;
  dBold        = ann.dBold;

  D            = ann.D;
  Delta        = ann.Delta;


  return *this;
}



void NeuralNetwork::CreateANN(boost::python::list ann){


   for(int i=0;i<len(ann);i++){
     int tmp = extract<int>(ann[i]);
     Npe.push_back(tmp);
   }

   Iteration = 0;

   Nlayers = Npe.size();

        srand((unsigned)time(0));

        MATRIX w; w.Init_Unit_Matrix(1.0);
        MATRIX b(3,1); b = 0.0;
        MATRIX d(3,1); d = 0.0;
        // 0-th matrixes are just ballast
        D.push_back(d);
        Delta.push_back(d);
        B.push_back(b);
        dB.push_back(b);
        dBcurr.push_back(b);
        dBold.push_back(b);

        W.push_back(w);
        w = 0.0;
        dW.push_back(w);
        dWcurr.push_back(w);        
        dWold.push_back(w);        


        // Init weights and biases (additional edges)
    for (int L=1;L<Nlayers;L++){
                MATRIX w(Npe[L],Npe[L-1]);
                MATRIX b(Npe[L],1);
                MATRIX d1(Npe[L],1);
                MATRIX d2(Npe[L],Npe[L]);

                        for(int i=0;i<Npe[L];i++){

                                for(int j=0;j<Npe[L-1];j++){

                                        w.M[i*Npe[L-1]+j] = 0.1*((rand()/(double)(RAND_MAX + 1.0))-0.5);

                                }
                                    b.M[i] = ((rand()/(double)(RAND_MAX + 1.0))-0.5);
                        }

                W.push_back(w);
                w = 0.0;
                dW.push_back(w);
                dWcurr.push_back(w);
                dWold.push_back(w);
                B.push_back(b);
                b = 0.0;
                dB.push_back(b);
                dBcurr.push_back(b);
                dBold.push_back(b);
                Delta.push_back(d1);
                D.push_back(d2);
    }
        // max index of W as well as B is Nlayers-1

//------------------ Debugging ------------------------
        std::cout<<"The MLP has been created"<<std::endl;
        std::cout<<"Number of layers = "<<(W.size()-1)<<std::endl;

        std::cout<<"Initial weight matrices:"<<std::endl;
        for(i=1;i<W.size();i++){              
                std::cout<<W[i]<<std::endl;
        }     
        std::cout<<"Initial biases:"<<std::endl;
        for(i=1;i<B.size();i++){
                std::cout<<B[i]<<std::endl;
        }
//------------------------------------------------------

}



int NeuralNetwork::ShowANN(){

// int NL= Nlayers-1;
// for(int L=1;L<=NL;L++){
 for(int L=1;L<W.size();L++){
     std::cout<<"W["<<L<<"] = "<<endl<<W[L];
 }
 for(L=1;L<B.size();L++){
     std::cout<<"B["<<L<<"] = "<<endl<<B[L];
 }

 return 0;

}

int NeuralNetwork::ImportANN(std::string filename){
 // Reads a 'filename' defining an ANN data
 // the 'filename' format is consistent with that of ExportANN format

 ifstream inp(filename.c_str());
 std::string str;
 std::string res;

 boost::regex re_int("\\d+",boost::regex::perl);
 boost::regex re_double("[-+]?(\\d+\\.\\d+|\\d+)([eE][-+]\\d+)?");
 boost::regex reW("W\\[\\d+]",boost::regex::perl);
 boost::regex reB("B\\[\\d+]",boost::regex::perl);

 smatch m;

   // Clean all arrays if they are not clean
   if(Npe.size()>0)     {Npe.clear();}
   if(W.size()>0)       {W.clear();  }
   if(dW.size()>0)      {dW.clear();  }
   if(dWcurr.size()>0)  {dWcurr.clear();  }
   if(dWold.size()>0)   {dWold.clear(); }
   if(B.size()>0)       {B.clear();  }
   if(dB.size()>0)      {dB.clear();  }
   if(dBcurr.size()>0)  {dBcurr.clear();  }
   if(dBold.size()>0)   {dBold.clear();}
   if(D.size()>0)       {D.clear();  }
   if(Delta.size()>0)   {Delta.clear();  }
   

 while(!inp.eof()){
 
 inp>>str;

 if(str=="Number_of_layers"){ inp>>str>>str; Nlayers = atoi(str.c_str());  }
 else if(str=="NN_architecture") { inp>>str>>str; 

   // Now decompose a string into integers
   boost::sregex_token_iterator I(str.begin(), str.end(), re_int, 0);
   boost::sregex_token_iterator J;
   while(I != J){    
        res = *I;
        int L = atoi(res.c_str());
        Npe.push_back(L);   I++;  }

   // And initialize weights and biases 
        srand((unsigned)time(0));

        MATRIX w; w.Init_Unit_Matrix(1.0);
        MATRIX b(3,1); b = 0.0;
        MATRIX d(3,1); d = 0.0;
        // 0-th matrixes are just ballast
        D.push_back(d);
        Delta.push_back(d);
        B.push_back(b);
        dB.push_back(b);
        dBcurr.push_back(b);
        dBold.push_back(b);
        W.push_back(w);
        w = 0.0;
        dW.push_back(w);
        dWcurr.push_back(w);
        dWold.push_back(w);

        // Init weights and biases (additional edges)
        for (int L=1;L<Nlayers;L++){
                MATRIX w(Npe[L],Npe[L-1]);
                MATRIX b(Npe[L],1);
                MATRIX d1(Npe[L],1);
                MATRIX d2(Npe[L],Npe[L]);

                        for(int i=0;i<Npe[L];i++){

                                for(int j=0;j<Npe[L-1];j++){

                                        w.M[i*Npe[L-1]+j] = 0.1*((rand()/(double)(RAND_MAX + 1.0))-0.5);

                                }
                                    b.M[i] = ((rand()/(double)(RAND_MAX + 1.0))-0.5);
                        }

                W.push_back(w);
                w = 0.0;
                dW.push_back(w);
                dWcurr.push_back(w);
                dWold.push_back(w);
                B.push_back(b);
                b = 0.0;
                dB.push_back(b);
                dBcurr.push_back(b);
                dBold.push_back(b);
                Delta.push_back(d1);
                D.push_back(d2);


        }// for L

        // Also create Inputs, Outputs and perhaps Derivs data arrays 
        // without actual data - with just scale and shift parameters
        sz_x = Npe[0];
        sz_y = Npe[Nlayers-1];
        sz_d = sz_x * sz_y;

        DATA d0;  
        if(Inputs.size()==0) {  for(L=0;L<sz_x;L++){ Inputs.push_back(d0);   }      }
        if(Outputs.size()==0){  for(L=0;L<sz_y;L++){ Outputs.push_back(d0);  }      }
        if(Derivs.size()==0) {  for(L=0;L<sz_d;L++){ Derivs.push_back(d0);   }      }
        

  }// if NN_architecture
 
 else if(regex_match(str,m,reW)){ // Look for W[1], W[2], etc.

   boost::sregex_token_iterator i(str.begin(), str.end(), re_int, 0); // extract i from W[i]
   res = *i;
   int L = atoi(res.c_str());
   inp>>str;
   for(int k=0;k<(Npe[L]*Npe[L-1]);k++){
       inp>>str; W[L].M[k]=atof(str.c_str());
   }
 }// if W

 else if(regex_match(str,m,reB)){ // Look for B[1], B[2], etc.

   boost::sregex_token_iterator i(str.begin(), str.end(), re_int, 0); // extract i from B[i]
   res = *i;
   int L = atoi(res.c_str());
   inp>>str;
   for(int k=0;k<(Npe[L]);k++){
       inp>>str; B[L].M[k]=atof(str.c_str());
   }
 }// if B
 
 else if(str=="x_scale"){ inp>>str>>str;
   // Now decompose a string into doubles
   boost::sregex_token_iterator I(str.begin(), str.end(), re_double, 0);
   boost::sregex_token_iterator J;
   int indx = 0;
   while(I != J){
        res = *I;
        double L = atof(res.c_str());
        Inputs[indx].scale_factor = L;
        Inputs[indx].is_scale_factor = 1;             
        indx++;
        I++;  }
 }// if x_scale

 else if(str=="x_shift"){ inp>>str>>str;
   // Now decompose a string into doubles
   boost::sregex_token_iterator I(str.begin(), str.end(), re_double, 0);
   boost::sregex_token_iterator J;
   int indx = 0;
   while(I != J){
        res = *I;
        double L = atof(res.c_str());
        Inputs[indx].shift_amount = L;
        Inputs[indx].is_shift_amount = 1;
        indx++;
        I++;  }
 }// if x_shift

 else if(str=="y_scale"){ inp>>str>>str;
   // Now decompose a string into doubles
   boost::sregex_token_iterator I(str.begin(), str.end(), re_double, 0);
   boost::sregex_token_iterator J;
   int indx = 0;
   while(I != J){
        res = *I;
        double L = atof(res.c_str());
        Outputs[indx].scale_factor = L;
        Outputs[indx].is_scale_factor = 1;
        indx++;
        I++;  }
 }// if y_scale

 else if(str=="y_shift"){ inp>>str>>str;
   // Now decompose a string into doubles
   boost::sregex_token_iterator I(str.begin(), str.end(), re_double, 0);
   boost::sregex_token_iterator J;
   int indx = 0;
   while(I != J){
        res = *I;
        double L = atof(res.c_str());
        Outputs[indx].shift_amount = L;
        Outputs[indx].is_shift_amount = 1;
        indx++;
        I++;  }
 }// if y_shift

 else if(str=="scale_method"){ inp>>str>>str;
    scale_method = str;
 }



 }// while !inp.eof()

//-------- Post-processing ---------------
// Update scaling factor for derivatives
  int i,j;

  for(i=0;i<sz_y;i++){
      for(j=0;j<sz_x;j++){
          // scale dYi/dXj
          Derivs[sz_x*i+j].scale_factor = (Outputs[i].scale_factor/Inputs[j].scale_factor);
      }// for j
  }// for i

 inp.close();

 return 0;

}


int NeuralNetwork::ExportANN(std::string filename){
 // Writes the ANN internals to the 'filename'

 ofstream out(filename.c_str()); 


 int NL= Nlayers-1;
 out<<"Number_of_layers = "<<Nlayers<<"\n";

 out<<"NN_architecture = [";
 for (int L=0;L<=NL;L++){
 out<<Npe[L]<<",";
 }
 out<<"]\n";

 for( L=1;L<=NL;L++){
     out<<"W["<<L<<"] = "<<endl<<W[L];
     out<<"B["<<L<<"] = "<<endl<<B[L];
 }

 int j;

 out<<"x_scale = [";
 for(j=0;j<Inputs.size();j++){  out<<Inputs[j].scale_factor<<",";}
 out<<"]\n";

 out<<"x_shift = [";
 for(j=0;j<Inputs.size();j++){  out<<Inputs[j].shift_amount<<",";}
 out<<"]\n";

 out<<"y_scale = [";
 for(j=0;j<Outputs.size();j++){  out<<Outputs[j].scale_factor<<",";}
 out<<"]\n";

 out<<"y_shift = [";
 for(j=0;j<Outputs.size();j++){  out<<Outputs[j].shift_amount<<",";}
 out<<"]\n";

 out<<"scale_method = "<<scale_method<<"\n";


 out.close();
 
 
 
 return 0;
}


void NeuralNetwork::save(boost::property_tree::ptree& pt,std::string path){

  if(is_learning_method){  ::save(pt,path+".learning_method",learning_method);    }
  if(is_learning_rate){  ::save(pt,path+".learning_rate",learning_rate);    }
  if(is_momentum_term){  ::save(pt,path+".momentum_term",momentum_term);    }
  if(is_epoch_size){  ::save(pt,path+".epoch_size",epoch_size);    }
  if(is_iterations_in_cycle){  ::save(pt,path+".iterations_in_cycle",iterations_in_cycle);    }
  if(is_grad_weight){  ::save(pt,path+".grad_weight",grad_weight);    }
  if(is_weight_decay){  ::save(pt,path+".weight_decay",weight_decay);    }
  if(is_norm_exp){  ::save(pt,path+".norm_exp",norm_exp);    }
  if(is_a_plus){  ::save(pt,path+".a_plus",a_plus);    }
  if(is_a_minus){  ::save(pt,path+".a_minus",a_minus);    }

  ::save(pt,path+".scale_method",scale_method);
  ::save(pt,path+".Iteration",Iteration);
  ::save(pt,path+".Cycle",Cycle);
  ::save(pt,path+".derivs_flag",derivs_flag);
//  ::save(pt,path+".TrainData",TrainData);
//  ::save(pt,path+".Recall",Recall);
//  ::save(pt,path+".Inputs",Inputs);
//  ::save(pt,path+".Outputs",Outputs);
//  ::save(pt,path+".Derivs",Derivs);
  ::save(pt,path+".num_of_patterns",num_of_patterns);
  ::save(pt,path+".sz_x",sz_x);
  ::save(pt,path+".sz_y",sz_y);
  ::save(pt,path+".sz_d",sz_d);
  ::save(pt,path+".Nlayers",Nlayers);
  ::save(pt,path+".Npe",Npe);
  ::save(pt,path+".B",B);
  ::save(pt,path+".dB",dB);
  ::save(pt,path+".dBcurr",dBcurr);
  ::save(pt,path+".dBold",dBold);
  ::save(pt,path+".W",W);
  ::save(pt,path+".dW",dW);
  ::save(pt,path+".dWcurr",dWcurr);
  ::save(pt,path+".dWold",dWold);
  ::save(pt,path+".D",D);
  ::save(pt,path+".Delta",Delta);


}

void save(boost::property_tree::ptree& pt,std::string path,vector<NeuralNetwork>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".NeuralNetwork"+rt);
  }
}


void NeuralNetwork::load(boost::property_tree::ptree& pt,std::string path,int& status){

  int st;
  status = 0;

  ::load(pt,path+".learning_method",learning_method,is_learning_method); if(is_learning_method==1) { status=1;}
  ::load(pt,path+".learning_rate",learning_rate,is_learning_rate); if(is_learning_rate==1) { status=1;}
  ::load(pt,path+".momentum_term",momentum_term,is_momentum_term); if(is_momentum_term==1) { status=1;}
  ::load(pt,path+".epoch_size",epoch_size,is_epoch_size); if(is_epoch_size==1) { status=1;}
  ::load(pt,path+".iterations_in_cycle",iterations_in_cycle,is_iterations_in_cycle); if(is_iterations_in_cycle==1) { status=1;}
  ::load(pt,path+".grad_weight",grad_weight,is_grad_weight); if(is_grad_weight==1) { status=1;}
  ::load(pt,path+".weight_decay",weight_decay,is_weight_decay); if(is_weight_decay==1) { status=1;}
  ::load(pt,path+".norm_exp",norm_exp,is_norm_exp); if(is_norm_exp==1) { status=1;}
  ::load(pt,path+".a_plus",a_plus,is_a_plus); if(is_a_plus==1) { status=1;}
  ::load(pt,path+".a_minus",a_minus,is_a_minus); if(is_a_minus==1) { status=1;}

  ::load(pt,path+".scale_method",scale_method,st); if(st==1) { status=1;}
  ::load(pt,path+".Iteration",Iteration,st); if(st==1) { status=1;}
  ::load(pt,path+".Cycle",Cycle,st); if(st==1) { status=1;}
  ::load(pt,path+".derivs_flag",derivs_flag,st); if(st==1) { status=1;}
//  ::load(pt,path+".TrainData",TrainData,st); if(st==1) { status=1;}
//  ::load(pt,path+".Recall",Recall,st); if(st==1) { status=1;}
//  ::load(pt,path+".Inputs",Inputs,st); if(st==1) { status=1;}
//  ::load(pt,path+".Outputs",Outputs,st); if(st==1) { status=1;}
//  ::load(pt,path+".Derivs",Derivs,st); if(st==1) { status=1;}
  ::load(pt,path+".num_of_patterns",num_of_patterns,st); if(st==1) { status=1;}
  ::load(pt,path+".sz_x",sz_x,st); if(st==1) { status=1;}
  ::load(pt,path+".sz_y",sz_y,st); if(st==1) { status=1;}
  ::load(pt,path+".sz_d",sz_d,st); if(st==1) { status=1;}
  ::load(pt,path+".Nlayers",Nlayers,st); if(st==1) { status=1;}
  ::load(pt,path+".Npe",Npe,st); if(st==1) { status=1;}
  ::load(pt,path+".B",B,st); if(st==1) { status=1;}
  ::load(pt,path+".dB",dB,st); if(st==1) { status=1;}
  ::load(pt,path+".dBcurr",dBcurr,st); if(st==1) { status=1;}
  ::load(pt,path+".dBold",dBold,st); if(st==1) { status=1;}
  ::load(pt,path+".W",W,st); if(st==1) { status=1;}
  ::load(pt,path+".dW",dW,st); if(st==1) { status=1;}
  ::load(pt,path+".dWcurr",dWcurr,st); if(st==1) { status=1;}
  ::load(pt,path+".dWold",dWold,st); if(st==1) { status=1;}
  ::load(pt,path+".D",D,st); if(st==1) { status=1;}
  ::load(pt,path+".Delta",Delta,st); if(st==1) { status=1;}



}

void load(boost::property_tree::ptree& pt,std::string path,vector<NeuralNetwork>& vt,int& status){
  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      NeuralNetwork x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}



}// namespace libann
}// namespace libmmath


