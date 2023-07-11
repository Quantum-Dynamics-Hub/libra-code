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

#if defined(USING_PCH)
#include "../pch.h"
#else

#include <ctime> 
#include <cstdlib>
#include <boost/regex.hpp>

#endif 

#include "NeuralNetwork.h"


/// liblibra namespace
namespace liblibra{

using namespace std;
using namespace boost;
using namespace libio;

namespace bp = boost::python;


// The signature of boost::property_tree::xml_parser::write_xml() changed in Boost 1.56
// See https://github.com/PointCloudLibrary/pcl/issues/864
#include <boost/version.hpp>
#if (BOOST_VERSION >= 105600)
  typedef boost::property_tree::xml_writer_settings<std::string> xml_writer_settings;
#else
  typedef boost::property_tree::xml_writer_settings<char> xml_writer_settings;
#endif



/// libann namespace
namespace libann{


void NeuralNetwork::allocate(vector<int>& arch){
/**  
  arch  - architecture of the ANN: the number of units in each layer

  The main variables we care to initialize here are the weights (W), biases (B), and their 
  derivatives (dW, dB)

  L                   0                   1                       ....             NL = Nlayers - 1

  W, dW, dWold      [junk]              W[1]                                        W[NL]

  B, dB, dBold      [junk]              B[1]                                        B[NL]

  Y             [Y[0]=input]      [ f(W[1]*Y[0] + B[1]) ]           [ output = f(W[NL]*Y[NL-1] + B[NL]) ]

  
*/

   int L;

   Nlayers = arch.size();
   Npe = arch;

   sz_x = arch[0];
   sz_y = arch[Nlayers-1];

   MATRIX w(Npe[0],Npe[0]); w.Init_Unit_Matrix(1.0);  // weights
   MATRIX b(Npe[0],1); b = 0.0;                       // biases
   MATRIX d(Npe[0],1); d = 0.0;                       // just general matrices

   // 0-th matrixes are just junk
   B.push_back(b);
   grad_b.push_back(b);
   grad_b_old.push_back(b);
   dB.push_back(b);
   dBold.push_back(b);

   W.push_back(w);
   w = 0.0;
   grad_w.push_back(w);
   grad_w_old.push_back(w);
   dW.push_back(w);
   dWold.push_back(w);

   // Init weights and biases (additional edges)
   for(L=1;L<Nlayers;L++){

      MATRIX w(Npe[L],Npe[L-1]);
      MATRIX b(Npe[L],1);

      B.push_back(b);
      grad_b.push_back(b);
      grad_b_old.push_back(b);
      dB.push_back(b);
      dBold.push_back(b);

      W.push_back(w);
      grad_w.push_back(w);
      grad_w_old.push_back(w);
      dW.push_back(w);
      dWold.push_back(w);


   } // max index of W as well as B is Nlayers-1

   //------------------ Debugging ------------------------
   std::cout<<"A multi-layer perceptron has been created\n";
   std::cout<<"Number of layers = "<<Nlayers<<std::endl;
   std::cout<<"The architecture of the MLP is: [";
   for(L=0;L<Nlayers-1;L++){  std::cout<<Npe[L]<<" , "; }
   std::cout<<Npe[Nlayers-1]<<"]\n";
   //------------------------------------------------------

}

NeuralNetwork::NeuralNetwork(vector<int>& arch){

  allocate(arch);
}

NeuralNetwork::NeuralNetwork(std::string filename){

  load(filename);
/*
  using boost::property_tree::ptree;
  ptree pt;

  std::ifstream ifs(filename, std::ifstream::in);
  read_xml(ifs, pt);
  ifs.close();

  int st;
  load(pt,"ANN", st);
*/
}



NeuralNetwork::NeuralNetwork(const NeuralNetwork& ann){

  Nlayers      = ann.Nlayers;
  Npe          = ann.Npe;

  W            = ann.W;
  grad_w       = ann.grad_w;
  grad_w_old   = ann.grad_w_old;
  dW           = ann.dW;
  dWold        = ann.dWold;

  B            = ann.B;
  grad_b       = ann.grad_b;
  grad_b_old   = ann.grad_b_old;
  dB           = ann.dB;
  dBold        = ann.dBold;

  D            = ann.D;
  Delta        = ann.Delta;

}


NeuralNetwork& NeuralNetwork::operator=(const NeuralNetwork& ann){

  Nlayers      = ann.Nlayers;
  Npe          = ann.Npe;

  W            = ann.W;
  grad_w       = ann.grad_w;
  grad_w_old   = ann.grad_w_old;
  dW           = ann.dW;
  dWold        = ann.dWold;

  B            = ann.B;
  grad_b       = ann.grad_b;
  grad_b_old   = ann.grad_b_old;
  dB           = ann.dB;
  dBold        = ann.dBold;

  D            = ann.D;
  Delta        = ann.Delta;


  return *this;
}





int NeuralNetwork::show(){
  int L;
  for(L=1;L<W.size();L++){ std::cout<<"W["<<L<<"] = "<<endl<<W[L]; }
  for(L=1;L<B.size();L++){ std::cout<<"B["<<L<<"] = "<<endl<<B[L]; }

  return 0;
}


void NeuralNetwork::save(boost::property_tree::ptree& pt,std::string path){

  libio::save(pt,path+".Nlayers",Nlayers);
  libio::save(pt,path+".Npe",Npe);
  libio::save(pt,path+".sz_x",sz_x);
  libio::save(pt,path+".sz_y",sz_y);

  liblinalg::save(pt,path+".B",B);
  liblinalg::save(pt,path+".grad_b",grad_b);
  liblinalg::save(pt,path+".grad_b_old",grad_b_old);
  liblinalg::save(pt,path+".dB",dB);
  liblinalg::save(pt,path+".dBold",dBold);

  liblinalg::save(pt,path+".W",W);
  liblinalg::save(pt,path+".grad_w",grad_w);
  liblinalg::save(pt,path+".grad_w_old",grad_w_old);
  liblinalg::save(pt,path+".dW",dW);
  liblinalg::save(pt,path+".dWold",dWold);

  liblinalg::save(pt,path+".D",D);
  liblinalg::save(pt,path+".Delta",Delta);

}

void NeuralNetwork::save(std::string filename){
/**
  \brief Save the ANN info into the XML file
  \param[in] The output file name 
  This is done by forming a property tree and then writing it as XML, using Boost write_xml function
*/

  using boost::property_tree::ptree;
  ptree pt;

  save(pt,"ANN");

//  boost::property_tree::xml_writer_settings<char> settings(' ', 4);
//  write_xml(filename, pt, std::locale(), settings);
  write_xml(filename, pt, std::locale(), xml_writer_settings('\t', 1));

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

  libio::load(pt,path+".Nlayers",Nlayers,st); if(st==1) { status=1;}
  libio::load(pt,path+".Npe",Npe,st); if(st==1) { status=1;}
  libio::load(pt,path+".sz_x",sz_x,st); if(st==1) { status=1;}
  libio::load(pt,path+".sz_y",sz_y,st); if(st==1) { status=1;}

  liblinalg::load(pt,path+".B",B,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".grad_b",grad_b,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".grad_b_old",grad_b_old,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".dB",dB,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".dBold",dBold,st); if(st==1) { status=1;}

  liblinalg::load(pt,path+".W",W,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".grad_w",grad_w,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".grad_w_old",grad_w_old,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".dW",dW,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".dWold",dWold,st); if(st==1) { status=1;}

  liblinalg::load(pt,path+".D",D,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".Delta",Delta,st); if(st==1) { status=1;}

}


int NeuralNetwork::load(std::string filename){
/**
  \brief Save the ANN info into the XML file
  \param[in] The output file name 
  This is done by forming a property tree and then writing it as XML, using Boost write_xml function
*/

  using boost::property_tree::ptree;
  ptree pt;

  std::ifstream ifs(filename, std::ifstream::in);
  read_xml(ifs, pt);
  ifs.close();

  int status;
  load(pt,"ANN", status);

  return status;

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
}// namespace liblibra


