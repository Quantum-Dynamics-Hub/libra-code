#include "AO.h"
#include "../molint/libmolint.h"


namespace libqchem{

using namespace libmolint;

namespace libqobjects{


//------------------ Members of  AO class -------------------------
AO::AO(){
  init();
}


void AO::init(){

  is_element = 0;
  is_ao_shell = 0;
  is_ao_shell_type = 0;
  is_ao_name = 0;

  is_x_exp = 0;
  is_y_exp = 0;
  is_z_exp = 0;

  expansion_size = 0;           is_expansion_size = 1;
  is_primitives = 0;
  is_coefficients = 0;

}

void AO::clear(){
  if(primitives.size()>0) { primitives.clear(); }
  if(coefficients.size()>0){ coefficients.clear(); }

  init();    
}

void AO::add_primitive(double c,PrimitiveG g){   
  coefficients.push_back(c);  is_coefficients = 1;
  primitives.push_back(g);    is_primitives   = 1;
  expansion_size++;
}


AO::AO(const AO& g){
  clear();

  if(g.is_x_exp)   {  x_exp = g.x_exp; is_x_exp = 1; }
  if(g.is_y_exp)   {  y_exp = g.y_exp; is_y_exp = 1; }
  if(g.is_z_exp)   {  z_exp = g.z_exp; is_z_exp = 1; }
  if(g.is_element) {  element = g.element; is_element = 1; }
  if(g.is_ao_shell){  ao_shell = g.ao_shell; is_ao_shell = 1; }
  if(g.is_ao_shell_type){  ao_shell_type = g.ao_shell_type; is_ao_shell_type = 1; }
  if(g.is_ao_name) {  ao_name = g.ao_name; is_ao_name = 1; }

  if(g.is_expansion_size){  expansion_size = g.expansion_size; is_expansion_size = 1;  }
  if(g.is_primitives) {
    for(int i=0;i<g.primitives.size();i++){ primitives.push_back(g.primitives[i]);  }
    is_primitives = 1;
  }
  if(g.is_coefficients){
    for(int i=0;i<g.coefficients.size();i++){ coefficients.push_back(g.coefficients[i]);  }
    is_coefficients = 1;
  }

}

// Assignment operator
AO& AO::operator=(const AO& g){

  clear();

  if(g.is_x_exp)   {  x_exp = g.x_exp; is_x_exp = 1; }
  if(g.is_y_exp)   {  y_exp = g.y_exp; is_y_exp = 1; }
  if(g.is_z_exp)   {  z_exp = g.z_exp; is_z_exp = 1; }
  if(g.is_element) {  element = g.element; is_element = 1; }
  if(g.is_ao_shell){  ao_shell = g.ao_shell; is_ao_shell = 1; }
  if(g.is_ao_shell_type){  ao_shell_type = g.ao_shell_type; is_ao_shell_type = 1; }
  if(g.is_ao_name) {  ao_name = g.ao_name; is_ao_name = 1; }
  if(g.is_expansion_size){  expansion_size = g.expansion_size; is_expansion_size = 1;  }

  if(g.is_primitives) { 
    for(int i=0;i<g.primitives.size();i++){ primitives.push_back(g.primitives[i]); }
    is_primitives = 1; 
  }
  if(g.is_coefficients){
    for(int i=0;i<g.coefficients.size();i++){ coefficients.push_back(g.coefficients[i]); }
    is_coefficients = 1; 
  }

  return *this;

}


void AO::show_info(){

  std::cout<<"AO properties:"<<std::endl;

  if(is_element){ std::cout<<"element = "<<element<<std::endl; }
  if(is_ao_shell){std::cout<<"ao_shell = "<<ao_shell<<std::endl;}
  if(is_ao_shell_type){std::cout<<"ao_shell_type = "<<ao_shell_type<<std::endl;}
  if(is_ao_name){ std::cout<<"ao_name = "<<ao_name<<std::endl; }
  if(is_x_exp)   {std::cout<<"x_exp = "<<x_exp<<" unitless\n";   }
  if(is_y_exp)   {std::cout<<"y_exp = "<<y_exp<<" unitless\n";   }
  if(is_z_exp)   {std::cout<<"z_exp = "<<z_exp<<" unitless\n";;   }
  if(is_expansion_size){ std::cout<<"expansion_size = "<<expansion_size<<" primitives"<<std::endl;
  
    for(int i=0;i<expansion_size;i++){ 
      std::cout<<"coefficients["<<i<<"] = "<<coefficients[i]<<std::endl;
      if(is_primitives){  std::cout<<"primitives["<<i<<"] = "<<std::endl;   primitives[i].show_info();   }
    }// for i
  } 

  std::cout<<std::endl;

}

double AO::compute(VECTOR& pos){

  double res = 0.0;
  for(int i=0;i<expansion_size;i++){
    res += coefficients[i] * primitives[i].compute(pos); 
  }

  return res;
}

void AO::normalize(){

  double res = normalization_factor();

  // Now scale contraction coefficents c[i] -> c'[i] = N*c[i], so AO' is normalized
  // int(AO'*AO'dr) = 1

  for(int i=0;i<expansion_size;i++){
    if(is_primitives){    coefficients[i] *= res; }
  }

}


double AO::norm2(){

  // Calculate the normalization of entire contraction
  // AO = N*summ(c_i * prim[i]) , where prim - are not normalized
  //         i
  // Calculate the norm:  <AO(A)|AO(A)>

  double res = 0.0;
  for(int i=0;i<expansion_size;i++){
    for(int j=0;j<expansion_size;j++){
      res += coefficients[i] * coefficients[j] * gaussian_overlap(primitives[i],primitives[j],1); // assume contraction of normalized Gaussians
    }
  }
  return res;
}

double AO::norm1(){

  return sqrt(norm2());
}

double AO::normalization_factor(){

  return (1.0/sqrt(norm2()));
}


void AO::shift_position(const VECTOR& dR){
// Move all primitives by vector dR
  for(int i=0;i<expansion_size;i++){  primitives[i].shift_position(dR);   }
}



///=======================================================================================================
///===================== Overload basic functions from libmolint to AO objects  ==========================

// Reference verions

double gaussian_overlap
( AO& AOa, AO& AOb,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
){

  dIdA = 0.0;
  dIdB = 0.0;
  
  VECTOR dida, didb;
  double w;
  double res = 0.0;
  for(int i=0;i<AOa.expansion_size;i++){
    for(int j=0;j<AOb.expansion_size;j++){

      w = AOa.coefficients[i] * AOb.coefficients[j];
      res += w * gaussian_overlap(AOa.primitives[i],AOb.primitives[j],is_normalize, is_derivs, dida, didb, auxd, n_aux); 
      dIdA += w * dida;
      dIdB += w * didb;

    }// j
  }// i

  return res;
}

double gaussian_overlap
( AO& AOa, AO& AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdB ){

  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(5);
  for(i=0;i<5;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_overlap(AOa, AOb, is_normalize, is_derivs, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<5;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}

boost::python::list gaussian_overlap
( AO& AOa, AO& AOb,int is_normalize, int is_derivs){

  VECTOR dIdA, dIdB;
  double I = gaussian_overlap(AOa, AOb, is_normalize, is_derivs, dIdA, dIdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dIdA);
    res.append(dIdB);
  }

  return res;
 
}


double gaussian_overlap(AO& AOa, AO& AOb,int is_normalize){

  VECTOR dIdA, dIdB;
  double res = gaussian_overlap(AOa, AOb, is_normalize, 0, dIdA, dIdB);
  return res;
}

double gaussian_overlap(AO& AOa, AO& AOb){

  double res = gaussian_overlap(AOa, AOb, 1);
  return res;

}



// Pointer versions

double gaussian_overlap
( AO* AOa, AO* AOb,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
){

  dIdA = 0.0;
  dIdB = 0.0;
  
  VECTOR dida, didb;
  double w;
  double res = 0.0;
  for(int i=0;i<AOa->expansion_size;i++){
    for(int j=0;j<AOb->expansion_size;j++){

      w = AOa->coefficients[i] * AOb->coefficients[j];
      res += w * gaussian_overlap(AOa->primitives[i],AOb->primitives[j],is_normalize, is_derivs, dida, didb, auxd, n_aux); 
      dIdA += w * dida;
      dIdB += w * didb;

    }// j
  }// i
  return res;
}

double gaussian_overlap
( AO* AOa, AO* AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdB ){

  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(5);
  for(i=0;i<5;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_overlap(AOa, AOb, is_normalize, is_derivs, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<5;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}


double gaussian_overlap(AO* AOa, AO* AOb,int is_normalize){

  VECTOR dIdA, dIdB;
  double res = gaussian_overlap(AOa, AOb, is_normalize, 0, dIdA, dIdB);
  return res;
}

double gaussian_overlap(AO* AOa, AO* AOb){

  double res = gaussian_overlap(AOa, AOb, 1);
  return res;

}




}// namespace libqobjects
}// namespace libqchem

