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
/**
  \file AO.cpp
  \brief The file implement: a) the AO class that represents atomic orbtals as a linear combination of 
  Gaussian primitives; b) related functions
    
*/

#include "AO.h"
#include "../molint/libmolint.h"

/// liblibra namespace
namespace liblibra{

using namespace libmolint;

/// libqobjects namespace
namespace libqobjects{


//------------------ Members of  AO class -------------------------
AO::AO(){
/**
  \brif The default constructor.
*/
  init();
}


void AO::init(){
/** 
  \brief Initialize the member data to the default values
*/

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
/**
  \brief Clears internal variables containing the primitive Gaussians and their coefficients in the contraction
*/
  if(primitives.size()>0) { primitives.clear(); }
  if(coefficients.size()>0){ coefficients.clear(); }

  init();    
}

void AO::add_primitive(double c,PrimitiveG g){   
/**
  \brief Add primitive Gaussian to the contraction.
  \param[in] c The coefficient with which the primitive Gaussian enters the contraction
  \param[in] g The primitive Gaussian object to add to AO
*/
  coefficients.push_back(c * g.normalization_factor());  is_coefficients = 1;
  primitives.push_back(g);    is_primitives   = 1;
  expansion_size++;
}


AO::AO(const AO& g){
/**
  \brief Copy constructor

  Only the properties defined in the source object will be copied to the destination object
*/

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
/**
  \brief Assignment operator

  Only the properties defined in the source object will be copied to the destination object  
*/


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
/**
  \brief Printing properties of the AO
*/

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
/**
  \brief Evaluate the AO function at given point 
  \param[in] pos The position at which the function is evaluated
  The result is returned.
*/

  double res = 0.0;
  for(int i=0;i<expansion_size;i++){
    res += coefficients[i] * primitives[i].compute(pos); 
  }

  return res;
}

void AO::normalize(){
/** 
  \brief Scale the coefficients of the contraction, so that the AO is normalized
*/

  double res = normalization_factor();

  // Now scale contraction coefficents c[i] -> c'[i] = N*c[i], so AO' is normalized
  // int(AO'*AO'dr) = 1

  for(int i=0;i<expansion_size;i++){
    if(is_primitives){    coefficients[i] *= res; }
  }

}


double AO::norm2(){
/**
  \brief The square of the norm of the entire contraction: <AO(A)|AO(A)>
  AO = N*summ(c_i * prim[i]) , where prim - are not normalized
           i
*/

  double res = 0.0;
  for(int i=0;i<expansion_size;i++){
    for(int j=0;j<expansion_size;j++){
      // coefficients already contain normalization factors for primitive Gaussians, so set the last parameter to 0
      res += coefficients[i] * coefficients[j] * gaussian_overlap(primitives[i],primitives[j],0); 
    }
  }
  return res;
}

double AO::norm1(){
/**
  \brief The magnitude (square root of the square of the norm) of the entire contraction: |<AO(A)|AO(A)>| = sqrt(<AO(A)|AO(A)>)
  AO = N*summ(c_i * prim[i]) , where prim - are not normalized
           i
*/

  return sqrt(norm2());
}

double AO::normalization_factor(){
/**
  \brief The normalization factor for the entire contraction: = 1.0/sqrt(<AO(A)|AO(A)>
*/

  return (1.0/sqrt(norm2()));
}


void AO::shift_position(VECTOR dR){
/**
  \brief Shift the center of all primitive by a given vector - Python-friendly
  \param[in] dR The translation vector
*/

  for(int i=0;i<expansion_size;i++){  primitives[i].shift_position(dR);   }
}

void AO::shift_position_const_ref(const VECTOR& dR){
/**
  \brief Shift the center of all primitive by a given vector
  \param[in] dR The translation vector
*/

  for(int i=0;i<expansion_size;i++){  primitives[i].shift_position_const_ref(dR);   }
}


void AO::set_position(VECTOR R_){
/**
  \brief Set the center of all primitives to a given vector - Python-friendly
  \param[in] R_ The new vector position
*/


  for(int i=0;i<expansion_size;i++){  primitives[i].set_position(R_);   }
}


void AO::set_position_const_ref(const VECTOR& R_){
/**
  \brief Set the center of all primitives to a given vector
  \param[in] R_ The new vector position
*/

  for(int i=0;i<expansion_size;i++){  primitives[i].set_position_const_ref(R_);   }
}


/*  This one is gonna MESS you a LOT!!!

 This is a good reminder that the Python uses references a lot !!!

  Alright!!! I keep the refernece version - for use in C++ only!!!
 
  Export the object version to Python to avoid messing up your calculations!!!

void AO::set_position(const VECTOR& R_){
// Move all primitives by vector dR
  for(int i=0;i<expansion_size;i++){  primitives[i].set_R(R_);   }
}
*/





//=======================================================================================================
//===================== Overload basic functions from libmolint to AO objects  ==========================

//=========================  Overlaps ================================

// Reference verions

double gaussian_overlap
( AO& AOa, AO& AOb,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
){
/**
  \brief Compute the overlap of two arbitrary AOs: <AO(A)|AO(B)>

  \param[in] AOa One primitive Gaussian
  \param[in] AOb Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized; is_normalize = 0 - no need for normalization (AOs are assumed to be normalized)
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each AO center
  \param[out] dIdA The derivative of the integral w.r.t. the coordinates of the first (A) AO (if is_derivs = 1)
  \param[out] dIdB The derivative of the integral w.r.t. the coordinates of the second (B) AO (if is_derivs = 1)
  \param[in,out] auxd The list of the pointers to pre-allocated pieces of memory (for variables of the double type)
  \param[in] n_aux The length of the array to which each of the auxd[i] pointers points.

*/


  dIdA = 0.0;
  dIdB = 0.0;
  
  VECTOR dida, didb;
  double w;
  double res = 0.0;
  for(int i=0;i<AOa.expansion_size;i++){
    for(int j=0;j<AOb.expansion_size;j++){

      w = AOa.coefficients[i] * AOb.coefficients[j];
      // use unnormalized primitives!
      res += w * gaussian_overlap(AOa.primitives[i],AOb.primitives[j], 0, is_derivs, dida, didb, auxd, n_aux); 
      dIdA += w * dida;
      dIdB += w * didb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa.normalization_factor() * AOb.normalization_factor();
    res *= nrm;
    dIdA *= nrm;
    dIdB *= nrm;
  }//

  return res;
}

double gaussian_overlap
( AO& AOa, AO& AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdB ){

  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_overlap(AOa, AOb, is_normalize, is_derivs, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
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
      // use unnormalized primitives!
      res += w * gaussian_overlap(AOa->primitives[i],AOb->primitives[j], 0, is_derivs, dida, didb, auxd, n_aux); 
      dIdA += w * dida;
      dIdB += w * didb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa->normalization_factor() * AOb->normalization_factor();
    res *= nrm;
    dIdA *= nrm;
    dIdB *= nrm;
  }//


  return res;
}

double gaussian_overlap
( AO* AOa, AO* AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdB ){

  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_overlap(AOa, AOb, is_normalize, is_derivs, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
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


//=========================  Moments ================================

// Reference verions

double gaussian_moment
( AO& AOa, PrimitiveG& G, AO& AOb,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB, vector<double*>& auxd,int n_aux
){

  dIdA = 0.0;
  dIdR = 0.0;
  dIdB = 0.0; 
  
  VECTOR dida, didr, didb;
  double w;
  double res = 0.0;
  for(int i=0;i<AOa.expansion_size;i++){
    for(int j=0;j<AOb.expansion_size;j++){

      w = AOa.coefficients[i] * AOb.coefficients[j];
      // use unnormalized primitives!
      res += w * gaussian_moment(AOa.primitives[i], G, AOb.primitives[j], 0, is_derivs, dida, didr, didb, auxd, n_aux); 
      dIdA += w * dida;
      dIdR += w * didr;
      dIdB += w * didb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa.normalization_factor() * AOb.normalization_factor();
    res *= nrm;
    dIdA *= nrm;
    dIdR *= nrm;
    dIdB *= nrm;
  }//
  return res;
}

double gaussian_moment
( AO& AOa, PrimitiveG& G, AO& AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB ){

  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_moment(AOa, G, AOb, is_normalize, is_derivs, dIdA, dIdR, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}

boost::python::list gaussian_moment
( AO& AOa, PrimitiveG& G, AO& AOb,int is_normalize, int is_derivs){

  VECTOR dIdA, dIdR, dIdB;
  double I = gaussian_moment(AOa, G, AOb, is_normalize, is_derivs, dIdA, dIdR, dIdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dIdA);
    res.append(dIdR);
    res.append(dIdB);
  }

  return res;
 
}


double gaussian_moment(AO& AOa, PrimitiveG& G, AO& AOb,int is_normalize){

  VECTOR dIdA, dIdR, dIdB;
  double res = gaussian_moment(AOa, G, AOb, is_normalize, 0, dIdA, dIdR, dIdB);
  return res;
}

double gaussian_moment(AO& AOa, PrimitiveG& G, AO& AOb){

  double res = gaussian_moment(AOa, G, AOb, 1);
  return res;

}



// Pointer versions

double gaussian_moment
( AO* AOa, PrimitiveG& G, AO* AOb,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB, vector<double*>& auxd,int n_aux
){

  dIdA = 0.0;
  dIdB = 0.0;
  
  VECTOR dida, didr, didb;
  double w;
  double res = 0.0;
  for(int i=0;i<AOa->expansion_size;i++){
    for(int j=0;j<AOb->expansion_size;j++){

      w = AOa->coefficients[i] * AOb->coefficients[j];
      // use unnormalized primitives!
      res += w * gaussian_moment(AOa->primitives[i], G, AOb->primitives[j], 0, is_derivs, dida, didr, didb, auxd, n_aux); 
      dIdA += w * dida;
      dIdR += w * didr;
      dIdB += w * didb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa->normalization_factor() * AOb->normalization_factor();
    res *= nrm;
    dIdA *= nrm;
    dIdR *= nrm;
    dIdB *= nrm;
  }//


  return res;
}

double gaussian_moment
( AO* AOa, PrimitiveG& G, AO* AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB ){

  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_moment(AOa, G, AOb, is_normalize, is_derivs, dIdA, dIdR, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}


double gaussian_moment(AO* AOa, PrimitiveG& G, AO* AOb,int is_normalize){

  VECTOR dIdA, dIdR, dIdB;
  double res = gaussian_moment(AOa, G, AOb, is_normalize, 0, dIdA, dIdR, dIdB);
  return res;
}

double gaussian_moment(AO* AOa, PrimitiveG& G, AO* AOb){

  double res = gaussian_moment(AOa, G, AOb, 1);
  return res;

}


//=========================  Pseudopotentials ================================

// Reference verions

double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO& AOa, AO& AOb,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB,
                   vector<double*>& auxd,int n_aux
                  ){
  dIdR = 0.0;
  dIdA = 0.0;
  dIdB = 0.0; 
  
  VECTOR dida, didr, didb;
  double w;
  double res = 0.0;
  for(int i=0;i<AOa.expansion_size;i++){
    for(int j=0;j<AOb.expansion_size;j++){

      w = AOa.coefficients[i] * AOb.coefficients[j];
      // use unnormalized primitives!
      res += w * pseudopot02(C0, C2, alp, R, AOa.primitives[i], AOb.primitives[j], 0, is_derivs, didr, dida, didb, auxd, n_aux); 
      dIdR += w * didr;
      dIdA += w * dida;
      dIdB += w * didb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa.normalization_factor() * AOb.normalization_factor();
    res *= nrm;
    dIdR *= nrm;
    dIdA *= nrm;
    dIdB *= nrm;
  }//
  return res;
}


double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO& AOa, AO& AOb,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB
                  ){

  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = pseudopot02(C0, C2, alp, R, AOa, AOb, is_normalize, is_derivs, dIdR, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<5;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}

boost::python::list pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO& AOa, AO& AOb, int is_normalize, int is_derivs
                  ){


  VECTOR dIdR, dIdA, dIdB;
  double I = pseudopot02(C0, C2, alp, R, AOa, AOb, is_normalize, is_derivs, dIdR, dIdA, dIdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dIdR);
    res.append(dIdA);
    res.append(dIdB);
  }

  return res;
 
}


double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO& AOa, AO& AOb, int is_normalize                   
                  ){

  VECTOR dIdR, dIdA, dIdB;
  double res = pseudopot02(C0, C2, alp, R, AOa, AOb, is_normalize, 0, dIdR, dIdA, dIdB);
  return res;
}

double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO& AOa, AO& AOb
                  ){

  double res = pseudopot02(C0, C2, alp, R, AOa, AOb, 1);
  return res;

}



// Pointer versions

double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO* AOa, AO* AOb,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB,
                   vector<double*>& auxd,int n_aux
                  ){
  dIdR = 0.0;
  dIdA = 0.0;
  dIdB = 0.0; 
  
  VECTOR dida, didr, didb;
  double w;
  double res = 0.0;
  for(int i=0;i<AOa->expansion_size;i++){
    for(int j=0;j<AOb->expansion_size;j++){

      w = AOa->coefficients[i] * AOb->coefficients[j];
      // use unnormalized primitives!
      res += w * pseudopot02(C0, C2, alp, R, AOa->primitives[i], AOb->primitives[j], 0, is_derivs, didr, dida, didb, auxd, n_aux); 
      dIdR += w * didr;
      dIdA += w * dida;
      dIdB += w * didb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa->normalization_factor() * AOb->normalization_factor();
    res *= nrm;
    dIdR *= nrm;
    dIdA *= nrm;
    dIdB *= nrm;
  }//
  return res;
}


double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO* AOa, AO* AOb,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB
                  ){
  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = pseudopot02(C0, C2, alp, R, AOa, AOb, is_normalize, is_derivs, dIdR, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}

boost::python::list pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO* AOa, AO* AOb, int is_normalize, int is_derivs
                  ){

  VECTOR dIdR, dIdA, dIdB;
  double I = pseudopot02(C0, C2, alp, R, AOa, AOb, is_normalize, is_derivs, dIdR, dIdA, dIdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dIdR);
    res.append(dIdA);
    res.append(dIdB);
  }

  return res;
 
}


double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO* AOa, AO* AOb, int is_normalize                   
                  ){
  VECTOR dIdR, dIdA, dIdB;
  double res = pseudopot02(C0, C2, alp, R, AOa, AOb, is_normalize, 0, dIdR, dIdA, dIdB);
  return res;
}

double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO* AOa, AO* AOb
                  ){
  double res = pseudopot02(C0, C2, alp, R, AOa, AOb, 1);
  return res;

}



//=========================  Multipoles ================================

// Reference verions
VECTOR transition_dipole_moment
( AO& AOa, AO& AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
){

  dMdA = 0.0;
  dMdB = 0.0;
  
  MATRIX3x3 dmda, dmdb;
  double w;
  VECTOR res; res = 0.0;
  for(int i=0;i<AOa.expansion_size;i++){
    for(int j=0;j<AOb.expansion_size;j++){

      w = AOa.coefficients[i] * AOb.coefficients[j];
      // use unnormalized primitives!
      res += w * transition_dipole_moment(AOa.primitives[i], AOb.primitives[j], 0, is_derivs, dmda, dmdb, auxd, n_aux); 
      dMdA += w * dmda;
      dMdB += w * dmdb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa.normalization_factor() * AOb.normalization_factor();
    res *= nrm;
    dMdA *= nrm;
    dMdB *= nrm;
  }//
  return res;
}

VECTOR transition_dipole_moment
( AO& AOa, AO& AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
){
  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  VECTOR res; res = transition_dipole_moment(AOa, AOb, is_normalize, is_derivs, dMdA, dMdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}

boost::python::list transition_dipole_moment
( AO& AOa, AO& AOb, int is_normalize,int is_derivs
){

  MATRIX3x3 dMdA, dMdB;
  VECTOR I; I = transition_dipole_moment(AOa, AOb, is_normalize, is_derivs, dMdA, dMdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dMdA);
    res.append(dMdB);
  }

  return res;
 
}

VECTOR transition_dipole_moment
( AO& AOa, AO& AOb, int is_normalize
){
  MATRIX3x3 dMdA, dMdB;
  VECTOR res; res = transition_dipole_moment(AOa, AOb, is_normalize, 0, dMdA, dMdB);
  return res;
}

VECTOR transition_dipole_moment
( AO& AOa, AO& AOb
){
  VECTOR res; res = transition_dipole_moment(AOa, AOb, 1);
  return res;
}


// Pointer versions
VECTOR transition_dipole_moment
( AO* AOa, AO* AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
){

  dMdA = 0.0;
  dMdB = 0.0;
  
  MATRIX3x3 dmda, dmdb;
  double w;
  VECTOR res; res = 0.0;
  for(int i=0;i<AOa->expansion_size;i++){
    for(int j=0;j<AOb->expansion_size;j++){

      w = AOa->coefficients[i] * AOb->coefficients[j];
      // use unnormalized primitives!
      res += w * transition_dipole_moment(AOa->primitives[i], AOb->primitives[j], 0, is_derivs, dmda, dmdb, auxd, n_aux); 
      dMdA += w * dmda;
      dMdB += w * dmdb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa->normalization_factor() * AOb->normalization_factor();
    res *= nrm;
    dMdA *= nrm;
    dMdB *= nrm;
  }//
  return res;
}

VECTOR transition_dipole_moment
( AO* AOa, AO* AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
){
  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  VECTOR res; res = transition_dipole_moment(AOa, AOb, is_normalize, is_derivs, dMdA, dMdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}

boost::python::list transition_dipole_moment
( AO* AOa, AO* AOb, int is_normalize,int is_derivs
){
  MATRIX3x3 dMdA, dMdB;
  VECTOR I; I = transition_dipole_moment(AOa, AOb, is_normalize, is_derivs, dMdA, dMdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dMdA);
    res.append(dMdB);
  }

  return res;
 
}

VECTOR transition_dipole_moment
( AO* AOa, AO* AOb, int is_normalize
){
  MATRIX3x3 dMdA, dMdB;
  VECTOR res; res = transition_dipole_moment(AOa, AOb, is_normalize, 0, dMdA, dMdB);
  return res;
}

VECTOR transition_dipole_moment
( AO* AOa, AO* AOb
){
  VECTOR res; res = transition_dipole_moment(AOa, AOb, 1);
  return res;
}



//=========================  Derivative coupling integrals ================================

// Reference verions
VECTOR derivative_coupling_integral
( AO& AOa, AO& AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
){

  dMdA = 0.0;
  dMdB = 0.0;
  
  MATRIX3x3 dmda, dmdb;
  double w;
  VECTOR res; res = 0.0;
  for(int i=0;i<AOa.expansion_size;i++){
    for(int j=0;j<AOb.expansion_size;j++){

      w = AOa.coefficients[i] * AOb.coefficients[j];
      // use unnormalized primitives!
      res += w * derivative_coupling_integral(AOa.primitives[i], AOb.primitives[j], 0, is_derivs, dmda, dmdb, auxd, n_aux); 
      dMdA += w * dmda;
      dMdB += w * dmdb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa.normalization_factor() * AOb.normalization_factor();
    res *= nrm;
    dMdA *= nrm;
    dMdB *= nrm;
  }//
  return res;
}

VECTOR derivative_coupling_integral
( AO& AOa, AO& AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
){
  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  VECTOR res; res = derivative_coupling_integral(AOa, AOb, is_normalize, is_derivs, dMdA, dMdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}

boost::python::list derivative_coupling_integral
( AO& AOa, AO& AOb, int is_normalize,int is_derivs
){
  MATRIX3x3 dMdA, dMdB;
  VECTOR I; I = derivative_coupling_integral(AOa, AOb, is_normalize, is_derivs, dMdA, dMdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dMdA);
    res.append(dMdB);
  }

  return res;
 
}

VECTOR derivative_coupling_integral
( AO& AOa, AO& AOb, int is_normalize
){
  MATRIX3x3 dMdA, dMdB;
  VECTOR res; res = derivative_coupling_integral(AOa, AOb, is_normalize, 0, dMdA, dMdB);
  return res;
}

VECTOR derivative_coupling_integral
( AO& AOa, AO& AOb
){
  VECTOR res; res = derivative_coupling_integral(AOa, AOb, 1);
  return res;
}


// Pointer versions
VECTOR derivative_coupling_integral
( AO* AOa, AO* AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
){

  dMdA = 0.0;
  dMdB = 0.0;
  
  MATRIX3x3 dmda, dmdb;
  double w;
  VECTOR res; res = 0.0;
  for(int i=0;i<AOa->expansion_size;i++){
    for(int j=0;j<AOb->expansion_size;j++){

      w = AOa->coefficients[i] * AOb->coefficients[j];
      // use unnormalized primitives!
      res += w * derivative_coupling_integral(AOa->primitives[i], AOb->primitives[j], 0, is_derivs, dmda, dmdb, auxd, n_aux); 
      dMdA += w * dmda;
      dMdB += w * dmdb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa->normalization_factor() * AOb->normalization_factor();
    res *= nrm;
    dMdA *= nrm;
    dMdB *= nrm;
  }//
  return res;
}

VECTOR derivative_coupling_integral
( AO* AOa, AO* AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
){
  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  VECTOR res; res = derivative_coupling_integral(AOa, AOb, is_normalize, is_derivs, dMdA, dMdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}

boost::python::list derivative_coupling_integral
( AO* AOa, AO* AOb, int is_normalize,int is_derivs
){
  MATRIX3x3 dMdA, dMdB;
  VECTOR I; I = derivative_coupling_integral(AOa, AOb, is_normalize, is_derivs, dMdA, dMdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dMdA);
    res.append(dMdB);
  }

  return res;
 
}

VECTOR derivative_coupling_integral
( AO* AOa, AO* AOb, int is_normalize
){
  MATRIX3x3 dMdA, dMdB;
  VECTOR res; res = derivative_coupling_integral(AOa, AOb, is_normalize, 0, dMdA, dMdB);
  return res;
}

VECTOR derivative_coupling_integral
( AO* AOa, AO* AOb
){
  VECTOR res; res = derivative_coupling_integral(AOa, AOb, 1);
  return res;
}






//=========================  Kinetic integrals ================================

// Reference verions

double kinetic_integral
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
      // use unnormalized primitives!
      res += w * kinetic_integral(AOa.primitives[i],AOb.primitives[j], 0, is_derivs, dida, didb, auxd, n_aux); 
      dIdA += w * dida;
      dIdB += w * didb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa.normalization_factor() * AOb.normalization_factor();
    res *= nrm;
    dIdA *= nrm;
    dIdB *= nrm;
  }//
  return res;
}

double kinetic_integral
( AO& AOa, AO& AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdB ){

  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = kinetic_integral(AOa, AOb, is_normalize, is_derivs, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}

boost::python::list kinetic_integral
( AO& AOa, AO& AOb,int is_normalize, int is_derivs){

  VECTOR dIdA, dIdB;
  double I = kinetic_integral(AOa, AOb, is_normalize, is_derivs, dIdA, dIdB);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dIdA);
    res.append(dIdB);
  }

  return res;
 
}


double kinetic_integral(AO& AOa, AO& AOb,int is_normalize){

  VECTOR dIdA, dIdB;
  double res = kinetic_integral(AOa, AOb, is_normalize, 0, dIdA, dIdB);
  return res;
}

double kinetic_integral(AO& AOa, AO& AOb){

  double res = kinetic_integral(AOa, AOb, 1);
  return res;

}



// Pointer versions

double kinetic_integral
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
      // use unnormalized primitives!
      res += w * kinetic_integral(AOa->primitives[i],AOb->primitives[j], 0, is_derivs, dida, didb, auxd, n_aux); 
      dIdA += w * dida;
      dIdB += w * didb;

    }// j
  }// i

  if(is_normalize){
    double nrm = AOa->normalization_factor() * AOb->normalization_factor();
    res *= nrm;
    dIdA *= nrm;
    dIdB *= nrm;
  }//


  return res;
}

double kinetic_integral
( AO* AOa, AO* AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdB ){

  // Allocate working memory
  int i;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = kinetic_integral(AOa, AOb, is_normalize, is_derivs, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;
}


double kinetic_integral(AO* AOa, AO* AOb,int is_normalize){

  VECTOR dIdA, dIdB;
  double res = kinetic_integral(AOa, AOb, is_normalize, 0, dIdA, dIdB);
  return res;
}

double kinetic_integral(AO* AOa, AO* AOb){

  double res = kinetic_integral(AOa, AOb, 1);
  return res;

}


//====================== Nuclear Attraction Integral ================================
// Versions with references

double nuclear_attraction_integral
( AO& AOa, AO& AOb, VECTOR& Rc, int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
){

  DA = 0.0;
  DB = 0.0;
  DC = 0.0;
  
  VECTOR dida, didb, didc;
  double w;
  double res = 0.0;
  for(int i=0;i<AOa.expansion_size;i++){
    for(int j=0;j<AOb.expansion_size;j++){

      w = AOa.coefficients[i] * AOb.coefficients[j];

      // use unnormalized primitives!
      res += w * nuclear_attraction_integral(AOa.primitives[i], AOb.primitives[j],
                 Rc , 0, is_derivs, dida, didb, didc, aux, n_aux, auxv, n_auxv); 

      DA += w * dida;
      DB += w * didb;
      DC += w * didc;
   
    }// j
  }// i

  if(is_normalize){
    double nrm = AOa.normalization_factor() * AOb.normalization_factor();
    res *= nrm;
    DA *= nrm;
    DB *= nrm;
    DC *= nrm;

  }//
  return res;
}

double nuclear_attraction_integral
( AO& AOa, AO& AOb, VECTOR& Rc, int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC
){

  // Allocate working memory
  int i;
  int n_aux = 20;
  int n_auxv = 10;
  vector<double*> auxd(20);
  for(i=0;i<20;i++){ auxd[i] = new double[n_aux]; }
  vector<VECTOR*> auxv(5);
  for(i=0;i<5;i++){ auxv[i] = new VECTOR[n_auxv]; }

  // Do computations
  double res = nuclear_attraction_integral(AOa, AOb, Rc,is_normalize, is_derivs, DA, DB, DC,
                                           auxd, n_aux, auxv, n_auxv);
  // Clean working memory
  for(i=0;i<20;i++){ delete [] auxd[i]; }  
  auxd.clear();
  for(i=0;i<5;i++){ delete [] auxv[i]; }  
  auxv.clear();
 
  return res;
}


boost::python::list nuclear_attraction_integral
( AO& AOa, AO& AOb, VECTOR& Rc, int is_normalize, int is_derivs
){

  VECTOR DA, DB, DC;
  double I = nuclear_attraction_integral(AOa, AOb, Rc, is_normalize, is_derivs, DA, DB, DC );

  boost::python::list res;
  res.append(I);
 
  if(is_derivs){
    res.append(DA);
    res.append(DB);
    res.append(DC);
  }

  return res;
 
}


double nuclear_attraction_integral
( AO& AOa, AO& AOb, VECTOR& Rc, int is_normalize
){
  VECTOR DA, DB, DC;
  double res = nuclear_attraction_integral(AOa, AOb, Rc, is_normalize, 0, DA, DB, DC);
  return res;
}

double nuclear_attraction_integral
( AO& AOa, AO& AOb, VECTOR& Rc){

  double res = nuclear_attraction_integral(AOa, AOb, Rc, 1);
  return res;
}


// Pointer versions
double nuclear_attraction_integral
( AO* AOa, AO* AOb, VECTOR& Rc, int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
){

  DA = 0.0;
  DB = 0.0;
  DC = 0.0;
  
  VECTOR dida, didb, didc;
  double w;
  double res = 0.0;
  for(int i=0;i<AOa->expansion_size;i++){
    for(int j=0;j<AOb->expansion_size;j++){

      w = AOa->coefficients[i] * AOb->coefficients[j];

      // use unnormalized primitives!
      res += w * nuclear_attraction_integral(AOa->primitives[i], AOb->primitives[j],
                 Rc , 0, is_derivs, dida, didb, didc, aux, n_aux, auxv, n_auxv); 

      DA += w * dida;
      DB += w * didb;
      DC += w * didc;
   
    }// j
  }// i

  if(is_normalize){
    double nrm = AOa->normalization_factor() * AOb->normalization_factor();
    res *= nrm;
    DA *= nrm;
    DB *= nrm;
    DC *= nrm;

  }//
  return res;
}

double nuclear_attraction_integral
( AO* AOa, AO* AOb, VECTOR& Rc, int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC
){

  // Allocate working memory
  int i;
  int n_aux = 20;
  int n_auxv = 10;
  vector<double*> auxd(20);
  for(i=0;i<20;i++){ auxd[i] = new double[n_aux]; }
  vector<VECTOR*> auxv(5);
  for(i=0;i<5;i++){ auxv[i] = new VECTOR[n_auxv]; }

  // Do computations
  double res = nuclear_attraction_integral(AOa, AOb, Rc,is_normalize, is_derivs, DA, DB, DC,
                                           auxd, n_aux, auxv, n_auxv);
  // Clean working memory
  for(i=0;i<20;i++){ delete [] auxd[i]; }  
  auxd.clear();
  for(i=0;i<5;i++){ delete [] auxv[i]; }  
  auxv.clear();
 
  return res;
}


boost::python::list nuclear_attraction_integral
( AO* AOa, AO* AOb, VECTOR& Rc, int is_normalize, int is_derivs
){
  VECTOR DA, DB, DC;
  double I = nuclear_attraction_integral(AOa, AOb, Rc, is_normalize, is_derivs, DA, DB, DC );

  boost::python::list res;
  res.append(I);
 
  if(is_derivs){
    res.append(DA);
    res.append(DB);
    res.append(DC);
  }

  return res;
 
}


double nuclear_attraction_integral
( AO* AOa, AO* AOb, VECTOR& Rc, int is_normalize
){
  VECTOR DA, DB, DC;
  double res = nuclear_attraction_integral(AOa, AOb, Rc, is_normalize, 0, DA, DB, DC);
  return res;
}

double nuclear_attraction_integral
( AO* AOa, AO* AOb, VECTOR& Rc){

  double res = nuclear_attraction_integral(AOa, AOb, Rc, 1);
  return res;
}



//====================== Electron Repulsion Integral ================================
// Versions with references
double electron_repulsion_integral
( AO& AOa, AO& AOb, AO& AOc, AO& AOd, int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
){

  DA = 0.0;
  DB = 0.0;
  DC = 0.0;
  DD = 0.0;

  
  VECTOR dida, didb, didc, didd;
  double w;
  double res = 0.0;
  for(int i=0;i<AOa.expansion_size;i++){
    for(int j=0;j<AOb.expansion_size;j++){
      for(int k=0;k<AOc.expansion_size;k++){
        for(int l=0;l<AOd.expansion_size;l++){

          w = AOa.coefficients[i] * AOb.coefficients[j] * AOc.coefficients[k] * AOd.coefficients[l];

          // use unnormalized primitives!
          res += w * electron_repulsion_integral(AOa.primitives[i], AOb.primitives[j],
                  AOc.primitives[k], AOd.primitives[l], 0, is_derivs, dida, didb, didc, didd, aux, n_aux, auxv, n_auxv); 

          DA += w * dida;
          DB += w * didb;
          DC += w * didc;
          DD += w * didd;
   
        }// l
      }// k
    }// j
  }// i

  if(is_normalize){
    double nrm = AOa.normalization_factor() * AOb.normalization_factor() * AOc.normalization_factor() * AOd.normalization_factor();
    res *= nrm;
    DA *= nrm;
    DB *= nrm;
    DC *= nrm;
    DD *= nrm;

  }//
  return res;
}


double electron_repulsion_integral
( AO& AOa, AO& AOb, AO& AOc, AO& AOd,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD
){

  // Allocate working memory
  int i;
  int n_aux = 40;
  int n_auxv = 40;
  vector<double*> auxd(30);
  for(i=0;i<30;i++){ auxd[i] = new double[n_aux]; }
  vector<VECTOR*> auxv(5);
  for(i=0;i<5;i++){ auxv[i] = new VECTOR[n_auxv]; }

  // Do computations
  double res = electron_repulsion_integral(AOa, AOb, AOc, AOd, 
                                           is_normalize, is_derivs, DA, DB, DC, DD,
                                           auxd, n_aux, auxv, n_auxv);
  // Clean working memory
  for(i=0;i<30;i++){ delete [] auxd[i]; }  
  auxd.clear();
  for(i=0;i<5;i++){ delete [] auxv[i]; }  
  auxv.clear();
 
  return res;
}


boost::python::list electron_repulsion_integral
( AO& AOa, AO& AOb, AO& AOc, AO& AOd,
  int is_normalize, int is_derivs
){


  VECTOR DA, DB, DC, DD;
  double I = electron_repulsion_integral(AOa, AOb, AOc, AOd, is_normalize, is_derivs, DA, DB, DC, DD);

  boost::python::list res;
  res.append(I);
 
  if(is_derivs){
    res.append(DA);
    res.append(DB);
    res.append(DC);
    res.append(DD);
  }

  return res;
 
}


double electron_repulsion_integral
( AO& AOa, AO& AOb, AO& AOc, AO& AOd, int is_normalize
){

  VECTOR DA, DB, DC, DD;
  double res = electron_repulsion_integral(AOa, AOb, AOc, AOd, is_normalize, 0, DA, DB, DC, DD);
  return res;
}

double electron_repulsion_integral
( AO& AOa, AO& AOb, AO& AOc, AO& AOd){

  double res = electron_repulsion_integral(AOa, AOb, AOc, AOd, 1);
  return res;
}


// Pointer versions

double electron_repulsion_integral
( AO* AOa, AO* AOb, AO* AOc, AO* AOd,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
){

  DA = 0.0;
  DB = 0.0;
  DC = 0.0;
  DD = 0.0;
  
  VECTOR dida, didb, didc, didd;
  double w;
  double res = 0.0;
  for(int i=0;i<AOa->expansion_size;i++){
    for(int j=0;j<AOb->expansion_size;j++){
      for(int k=0;k<AOc->expansion_size;k++){
        for(int l=0;l<AOd->expansion_size;l++){

          w = AOa->coefficients[i] * AOb->coefficients[j]*AOc->coefficients[k] * AOd->coefficients[l];

          // use unnormalized primitives!
          res += w * electron_repulsion_integral(AOa->primitives[i], AOb->primitives[j],
          AOc->primitives[k], AOd->primitives[l], 0, is_derivs, dida, didb, didc, didd, aux, n_aux, auxv, n_auxv); 

          DA += w * dida;
          DB += w * didb;
          DC += w * didc;
          DD += w * didd;

        }// l
      }// k
    }// j
  }// i

  if(is_normalize){
    double nrm = AOa->normalization_factor() * AOb->normalization_factor() * AOc->normalization_factor() * AOd->normalization_factor();
    res *= nrm;
    DA *= nrm;
    DB *= nrm;
    DC *= nrm;
    DD *= nrm;

  }//
  return res;
}


double electron_repulsion_integral
( AO* AOa, AO* AOb, AO* AOc, AO* AOd,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD
){

  // Allocate working memory
  int i;
  int n_aux = 40;
  int n_auxv = 40;
  vector<double*> auxd(30);
  for(i=0;i<30;i++){ auxd[i] = new double[n_aux]; }
  vector<VECTOR*> auxv(5);
  for(i=0;i<5;i++){ auxv[i] = new VECTOR[n_auxv]; }

  // Do computations
  double res = electron_repulsion_integral(AOa, AOb, AOc, AOd, 
                                           is_normalize, is_derivs, DA, DB, DC, DD,
                                           auxd, n_aux, auxv, n_auxv);
  // Clean working memory
  for(i=0;i<30;i++){ delete [] auxd[i]; }  
  auxd.clear();
  for(i=0;i<5;i++){ delete [] auxv[i]; }  
  auxv.clear();
 
  return res;
}


boost::python::list electron_repulsion_integral
( AO* AOa, AO* AOb, AO* AOc, AO* AOd,
  int is_normalize, int is_derivs
){

  VECTOR DA, DB, DC, DD;
  double I = electron_repulsion_integral(AOa, AOb, AOc, AOd, is_normalize, is_derivs, DA, DB, DC, DD);

  boost::python::list res;
  res.append(I);
 
  if(is_derivs){
    res.append(DA);
    res.append(DB);
    res.append(DC);
    res.append(DD);
  }

  return res;
 
}


double electron_repulsion_integral
( AO* AOa, AO* AOb, AO* AOc, AO* AOd, int is_normalize
){
  VECTOR DA, DB, DC, DD;
  double res = electron_repulsion_integral(AOa, AOb, AOc, AOd, is_normalize, 0, DA, DB, DC, DD);
  return res;
}

double electron_repulsion_integral
( AO* AOa, AO* AOb, AO* AOc, AO* AOd){

  double res = electron_repulsion_integral(AOa, AOb, AOc, AOd, 1);
  return res;
}




}// namespace libqobjects
}// namespace liblibra

