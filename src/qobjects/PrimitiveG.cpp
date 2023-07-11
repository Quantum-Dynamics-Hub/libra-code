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
  \file PrimitiveG.cpp
  \brief The file implements: a) the PrimitiveG class that represents Gaussian primitives; b) related functions
    
*/

#include "PrimitiveG.h"
#include "../molint/libmolint.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;
using namespace libmolint;

/// libqobjects namespace
namespace libqobjects{


//----------------------- Members of PrimitiveG class -----------------------

int PrimitiveG::get_x_exp(){ 
/** 
  \brief Returns x exponent (x-projection of angular momentum quantum number)
*/
  if(is_x_exp){ return x_exp; }
   else{ std::cout<<"Error: x_exp is not set\n"; exit(0); } 
}

int PrimitiveG::get_y_exp(){
/** 
  \brief Returns y exponent (y-projection of angular momentum quantum number)
*/

  if(is_y_exp){ return y_exp; }
  else{ std::cout<<"Error: y_exp is not set\n"; exit(0); } 
}

int PrimitiveG::get_z_exp(){
/** 
  \brief Returns z exponent (z-projection of angular momentum quantum number)
*/

  if(is_z_exp){ return z_exp; }
  else{ std::cout<<"Error: z_exp is not set\n"; exit(0); } 
}

double PrimitiveG::get_alpha(){
/** 
  \brief Returns Gaussian exponent 
*/

  if(is_alpha){ return alpha; } else{ std::cout<<"Error: alpha is not set\n"; exit(0); } 
}

VECTOR PrimitiveG::get_R(){ 
/** 
  \brief Returns the position of the primitive Gaussian
*/

  if(is_R){ return R; } else{ std::cout<<"Error: R is not set\n"; exit(0); } 
}

double PrimitiveG::get_value(){
/** 
  \brief Returns the value of the primitive Gaussian function, computed at some point
  Note: this functions doesn't do the calculation, only returns the pre-computed value
*/

  if(is_value){ return value; } else{ std::cout<<"Error: value is not set\n"; exit(0); } 
}


void PrimitiveG::set_x_exp(int _x){ 
/** 
  \brief Sets the x exponent (x-projection of angular momentum quantum number)
  \param[in] _x The x exponent input
*/

  x_exp = _x; is_x_exp = 1; 
}

void PrimitiveG::set_y_exp(int _y){
/** 
  \brief Sets the y exponent (y-projection of angular momentum quantum number)
  \param[in] _y The y exponent input
*/

  y_exp = _y; is_y_exp = 1; 
}

void PrimitiveG::set_z_exp(int _z){
/** 
  \brief Sets the z exponent (z-projection of angular momentum quantum number)
  \param[in] _z The z exponent input
*/

  z_exp = _z; is_z_exp = 1; 
}

void PrimitiveG::set_alpha(double _alp){ 
/** 
  \brief Sets the Gaussian exponent, alpha
  \param[in] _alp The Gaussian exponent, alpha, input
*/

  alpha = _alp; is_alpha = 1; 
}

void PrimitiveG::set_position_const_ref(const VECTOR& _R){ 
/** 
  \brief Sets the position of the primitive
  \param[in] _R The input position
  Note: here we use the setting by the const reference, so changing _R outside, will not affect the 
  internal variable R - anyways, this version is not exported to Python, so don't worry.
  
*/

  R = _R; is_R = 1; 
}

void PrimitiveG::set_position(VECTOR _R){ 
/** 
  \brief Sets the position of the primitive
  \param[in] _R The input position
  Note: here the assignment made by value, that is the copy constructory is used,
  so changing _R outside, will not affect the internal variable R.
  
*/

  R = _R; is_R = 1; }


// General initialization
void PrimitiveG::init(int l,int m,int n,double alp,VECTOR& center){
/** 
  \brief Initialization with parameters
  \param[in] l x-exponent
  \param[in] m y-exponent
  \param[in] n z-exponent
  \param[in] alp Gaussian exponent
  \param[in] center The position of the primitive Gaussian. Note: the assignment is made by reference, so changing the 
  center variable outside *WILL* change the position of the primitive
*/

  x_exp = l;           is_x_exp = 1;
  y_exp = m;           is_y_exp = 1;
  z_exp = n;           is_z_exp = 1;
  alpha = alp;         is_alpha = 1;
  R = center;          is_R = 1;
  value = 0.0;         is_value = 0;
}

// Default constructor
PrimitiveG::PrimitiveG(){  
/** 
  \brief Default ctor: Simple S-type function at the global origin: R = (0,0,0)^T, alpha = 1.0
*/
  VECTOR r;
  init(0,0,0,1.0, r );
}

// Arbitrary L function constructor
PrimitiveG::PrimitiveG(int l,int m,int n,double alp,VECTOR& center){ 
/** 
  \brief Constructor with initialization
  \param[in] l x-exponent
  \param[in] m y-exponent
  \param[in] n z-exponent
  \param[in] alp Gaussian exponent
  \param[in] center The position of the primitive Gaussian. Note: the assignment is made by reference, so changing the 
  center variable outside *WILL* change the position of the primitive
*/

  init(l,m,n,alp, center );
}

// Copy constructor
PrimitiveG::PrimitiveG(const PrimitiveG& g){
/**
  \brief Copy constructor

  Only the properties defined in the source object will be copied to the destination object
*/

  // Default parameters
  VECTOR r;
  init(0,0,0,1.0, r );

  if(g.is_x_exp)   {  x_exp = g.x_exp; is_x_exp = 1; }
  if(g.is_y_exp)   {  y_exp = g.y_exp; is_y_exp = 1; }
  if(g.is_z_exp)   {  z_exp = g.z_exp; is_z_exp = 1; }
  if(g.is_alpha)   {  alpha = g.alpha; is_alpha = 1; }
  if(g.is_R)       {  R     = g.R;     is_R     = 1; }
  if(g.is_value)   {  value = g.value; is_value = 1; }

}

PrimitiveG& PrimitiveG::operator=(const PrimitiveG& g){
/**
  \brief Assignment operator

  Only the properties defined in the source object will be copied to the destination object  
*/

  // Default parameters
  VECTOR r;
  init(0,0,0,1.0, r );
 
  if(g.is_x_exp) {  x_exp = g.x_exp; is_x_exp = 1; }
  if(g.is_y_exp) {  y_exp = g.y_exp; is_y_exp = 1; }
  if(g.is_z_exp) {  z_exp = g.z_exp; is_z_exp = 1; }
  if(g.is_alpha) {  alpha = g.alpha; is_alpha = 1; }
  if(g.is_R)     {  R     = g.R;     is_R     = 1; }
  if(g.is_value) {  value = g.value; is_value = 1; }

  return *this;

}

void PrimitiveG::show_info(){
/**
  \brief Printing properties of the primitive Gaussian
*/

  std::cout<<"Primitive Gaussian properties:"<<std::endl;

  if(is_x_exp) {std::cout<<"x_exp = "<<x_exp<<" unitless"<<std::endl;   }
  if(is_y_exp) {std::cout<<"y_exp = "<<y_exp<<" unitless"<<std::endl;   }
  if(is_z_exp) {std::cout<<"z_exp = "<<z_exp<<" unitless"<<std::endl;   }
  if(is_alpha) {std::cout<<"alpha = "<<alpha<<" Bohr^-1"<<std::endl;   }
  if(is_R)     {std::cout<<"R = "<<R<<" Bohr"<<std::endl;   }
  if(is_value) {std::cout<<"value = "<<value<<" unitless"<<std::endl;   }

  std::cout<<std::endl;

}


double PrimitiveG::compute(VECTOR& pos){
/**
  \brief Evaluate the primitive Gaussian function at given point 
  \param[in] pos The position at which the function is evaluated
  The result is stored in the member-function value
*/

  VECTOR r; r = pos - R;
  double r2 = r.length2();
  value = FAST_POW(r.x,x_exp)*FAST_POW(r.y,y_exp)*FAST_POW(r.z,z_exp)*exp(-alpha * r2);
  is_value = 1;

  return value;
}


double PrimitiveG::norm2(){
/**
  \brief Compute the square norm of the given primitive Gaussian, <G(A)|G(A)>

*/

  double res = gaussian_norm2(x_exp,alpha) * gaussian_norm2(y_exp,alpha) * gaussian_norm2(z_exp,alpha);
  return res;

}// norm2

double PrimitiveG::norm1(){
/**
  \brief Compute the square norm of the given primitive Gaussian, sqrt(<G(A)|G(A)>)

*/

  return sqrt(norm2()); 

}// norm1

double PrimitiveG::normalization_factor(){
/**
  \brief Compute normalization factor for given primitive Gaussian

  double nom = pow((8.0*G_alpha),(x_exp+y_exp+z_exp))*FACTORIAL(x_exp)*FACTORIAL(y_exp)*FACTORIAL(z_exp);
  double denom = FACTORIAL(2*x_exp)*FACTORIAL(2*y_exp)*FACTORIAL(2*z_exp);
  res = pow((2.0*G_alpha/M_PI),0.75)*sqrt(nom/denom);

  Meaning: if N - is a result of this function and
  G(l,m,n,alp) = x^l * y^m * z^n * exp(-alp*r^2)  a primitive Gaussian
  then
  g = N * G(l,m,n,alp) - is normaized: integral(g,g) = 1.0

  = 1/sqrt(<G(A)|G(A)>)
*/

  double res = (1.0/sqrt(norm2()));
  return res;

}// norm

void PrimitiveG::shift_position_const_ref(const VECTOR& dR){
/**
  \brief Shift the center of the primitive by a given vector
  \param[in] dR The translation vector
*/

  R += dR; 
}

void PrimitiveG::shift_position(VECTOR dR){  
/**
  \brief Shift the center of the primitive by a given vector - Python-friendly
  \param[in] dR The translation vector
*/

  R += dR; 
}


//=======================================================================================================
//===================== Overload basic functions from libmolint to Gaussian objects  ====================

//======================== Overlaps =================================

double gaussian_overlap
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
){
/**
  \brief Compute the overlap of two arbitrary primitive Gaussians: <G(A)|G(B)>

  \param[in] GA One primitive Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized)
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each primitive center
  \param[out] dIdA The derivative of the integral w.r.t. the coordinates of the first (A) primitive Gaussian (if is_derivs = 1)
  \param[out] dIdB The derivative of the integral w.r.t. the coordinates of the second (B) primitive Gaussian (if is_derivs = 1)
  \param[in,out] auxd The list of the pointers to pre-allocated pieces of memory (for variables of the double type)
  \param[in] n_aux The length of the array to which each of the auxd[i] pointers points.

*/

  double res = 
  libmolint::gaussian_overlap
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dIdA, dIdB, auxd, n_aux
  );

  return res;
}

double gaussian_overlap
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB
){
/**
  \brief Compute the overlap of two arbitrary primitive Gaussians: <G(A)|G(B)> - Python-friendly

  \param[in] GA One primitive Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized)
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each primitive center
  \param[out] dIdA The derivative of the integral w.r.t. the coordinates of the first (A) primitive Gaussian (if is_derivs = 1)
  \param[out] dIdB The derivative of the integral w.r.t. the coordinates of the second (B) primitive Gaussian (if is_derivs = 1)

  This version does not require externally pre-allocated memory pointers. This approach may reduce the performance (especially if called
  frequently), but is more convenient - some sufficient memory is allocated/deallocated temporarily, internally. If there are 
  some problems, you may increase the amount of internally-allocated memory.

*/


  double res = 
  libmolint::gaussian_overlap
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dIdA, dIdB
  );

  return res;
}

boost::python::list gaussian_overlap
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs){
/**
  \brief Compute the overlap of two arbitrary primitive Gaussians: <G(A)|G(B)> - Python-friendly

  \param[in] GA One primitive Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized)
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each primitive center

  This version does not require externally pre-allocated memory pointers. This approach may reduce the performance (especially if called
  frequently), but is more convenient - some sufficient memory is allocated/deallocated temporarily, internally. If there are 
  some problems, you may increase the amount of internally-allocated memory.

  This version returns a list of results, res, such that res[0] - the overlap, res[1] - the derivativ of the integral w.r.t. the coordinates
  of the A-center (dI/dA), res[2] - derivative w.r.t. B, dI/dB. The elements res[1] and res[2] are available only if is_derivs = 1.
  Otherwise, only the overlap is computed, so only res[0] is available.

*/


  boost::python::list res;
  res = 
  libmolint::gaussian_overlap
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs
  );

  return res;
}

double gaussian_overlap
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize){
/**
  \brief Compute the overlap of two arbitrary primitive Gaussians: <G(A)|G(B)> - Python-friendly

  \param[in] GA One primitive Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized)

  This version does not require externally pre-allocated memory pointers. 
  This version does not compute derivatives. Only the integral is computed.
*/


  double res = 
  libmolint::gaussian_overlap
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize
  );

  return res;
}

double gaussian_overlap
( PrimitiveG& GA, PrimitiveG& GB){
/**
  \brief Compute the overlap of two arbitrary primitive Gaussians: <G(A)|G(B)> - Python-friendly

  \param[in] GA One primitive Gaussian
  \param[in] GB Second primitive Gaussian

  This version does not require externally pre-allocated memory pointers. 
  This version does not compute derivatives. Only the integral is computed.
  This version does perform normalization by the default.
*/


  double res = 
  libmolint::gaussian_overlap
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R
  );

  return res;
}


//======================== Moments =================================

double gaussian_moment
( PrimitiveG& GA, PrimitiveG& G, PrimitiveG& GB,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB, vector<double*>& auxd,int n_aux
){
/**
  \brief Compute the overlap of 3 arbitrary primitive Gaussians: <G(A)|G|G(B)> 
  Note that specific cases of such integral include dipole, multipole moments, and components of pseudopotential:
  <G(A)|R_c|G(B)>,   <G(A)|R_c * R'_c|G(B)>,   <G(A)|exp(-alp*(r-R_c))|G(B)>, etc.

  \param[in] GA One primitive Gaussian
  \param[in] G The function that defines the momoent - this function is assumed in the form of generalized Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized). The normalization only uses the factors for GA and GB!
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each primitive center
  \param[out] dIdA The derivative of the integral w.r.t. the coordinates of the first (A) primitive Gaussian (if is_derivs = 1)
  \param[out] dIdR The derivative of the integral w.r.t. the coordinates of the moment-defining primitive Gaussian (if is_derivs = 1)
  \param[out] dIdB The derivative of the integral w.r.t. the coordinates of the second (B) primitive Gaussian (if is_derivs = 1)
  \param[in,out] auxd The list of the pointers to pre-allocated pieces of memory (for variables of the double type)
  \param[in] n_aux The length of the array to which each of the auxd[i] pointers points.

*/


  double res = 
  libmolint::gaussian_moment
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    G.x_exp, G.y_exp,  G.z_exp,  G.alpha,  G.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dIdA, dIdR, dIdB, auxd, n_aux
  );

  return res;
}

double gaussian_moment
( PrimitiveG& GA, PrimitiveG& G, PrimitiveG& GB,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB
){
/**
  \brief Compute the overlap of 3 arbitrary primitive Gaussians: <G(A)|G|G(B)> - Python-friendly
  Note that specific cases of such integral include dipole, multipole moments, and components of pseudopotential:
  <G(A)|R_c|G(B)>,   <G(A)|R_c * R'_c|G(B)>,   <G(A)|exp(-alp*(r-R_c))|G(B)>, etc.

  \param[in] GA One primitive Gaussian
  \param[in] G The function that defines the momoent - this function is assumed in the form of generalized Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized) The normalization only uses the factors for GA and GB!
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each primitive center
  \param[out] dIdA The derivative of the integral w.r.t. the coordinates of the first (A) primitive Gaussian (if is_derivs = 1)
  \param[out] dIdR The derivative of the integral w.r.t. the coordinates of the moment-defining primitive Gaussian (if is_derivs = 1)
  \param[out] dIdB The derivative of the integral w.r.t. the coordinates of the second (B) primitive Gaussian (if is_derivs = 1)

  This version does not require externally pre-allocated memory pointers. This approach may reduce the performance (especially if called
  frequently), but is more convenient - some sufficient memory is allocated/deallocated temporarily, internally. If there are 
  some problems, you may increase the amount of internally-allocated memory.

*/


  double res = 
  libmolint::gaussian_moment
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    G.x_exp, G.y_exp,  G.z_exp,  G.alpha,  G.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dIdA, dIdR, dIdB
  );

  return res;
}

boost::python::list gaussian_moment
( PrimitiveG& GA, PrimitiveG& G, PrimitiveG& GB,int is_normalize, int is_derivs){
/**
  \brief Compute the overlap of 3 arbitrary primitive Gaussians: <G(A)|G|G(B)> - Python-friendly
  Note that specific cases of such integral include dipole, multipole moments, and components of pseudopotential:
  <G(A)|R_c|G(B)>,   <G(A)|R_c * R'_c|G(B)>,   <G(A)|exp(-alp*(r-R_c))|G(B)>, etc.

  \param[in] GA One primitive Gaussian
  \param[in] G The function that defines the momoent - this function is assumed in the form of generalized Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized). The normalization only uses the factors for GA and GB!
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each primitive center

  This version does not require externally pre-allocated memory pointers. This approach may reduce the performance (especially if called
  frequently), but is more convenient - some sufficient memory is allocated/deallocated temporarily, internally. If there are 
  some problems, you may increase the amount of internally-allocated memory.

  This version returns a list of results, res, such that res[0] - the overlap, res[1] - the derivative of the integral w.r.t. the coordinates
  of the A-center (dI/dA), res[2] - derivative w.r.t. the moment-generating primitive, res[3] - derivative w.r.t. B, dI/dB. 
  The elements res[1], res[2], res[3] are available only if is_derivs = 1.
  Otherwise, only the overlap is computed, so only res[0] is available.

*/


  boost::python::list res;
  res = 
  libmolint::gaussian_moment
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    G.x_exp, G.y_exp,  G.z_exp,  G.alpha,  G.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs
  );

  return res;
}

double gaussian_moment
( PrimitiveG& GA, PrimitiveG& G, PrimitiveG& GB,int is_normalize){
/**
  \brief Compute the overlap of 3 arbitrary primitive Gaussians: <G(A)|G|G(B)> - Python-friendly
  Note that specific cases of such integral include dipole, multipole moments, and components of pseudopotential:
  <G(A)|R_c|G(B)>,   <G(A)|R_c * R'_c|G(B)>,   <G(A)|exp(-alp*(r-R_c))|G(B)>, etc.

  \param[in] GA One primitive Gaussian
  \param[in] G The function that defines the momoent - this function is assumed in the form of generalized Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized). The normalization only uses the factors for GA and GB!

  This version does not require externally pre-allocated memory pointers. 
  This version does not compute derivatives. Only the integral is computed.

*/


  double res = 
  libmolint::gaussian_moment
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    G.x_exp, G.y_exp,  G.z_exp,  G.alpha,  G.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize
  );

  return res;
}

double gaussian_moment
( PrimitiveG& GA, PrimitiveG& G, PrimitiveG& GB){
/**
  \brief Compute the overlap of 3 arbitrary primitive Gaussians: <G(A)|G|G(B)> - Python-friendly
  Note that specific cases of such integral include dipole, multipole moments, and components of pseudopotential:
  <G(A)|R_c|G(B)>,   <G(A)|R_c * R'_c|G(B)>,   <G(A)|exp(-alp*(r-R_c))|G(B)>, etc.

  \param[in] GA One primitive Gaussian
  \param[in] G The function that defines the momoent - this function is assumed in the form of generalized Gaussian
  \param[in] GB Second primitive Gaussian

  This version does not require externally pre-allocated memory pointers. 
  This version does not compute derivatives. Only the integral is computed.
  This version does perform normalization by the default (only using the factors for GA and GB!)

*/


  double res = 
  libmolint::gaussian_moment
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    G.x_exp, G.y_exp,  G.z_exp,  G.alpha,  G.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R
  );

  return res;
}



//======================== Pseudopotentials =================================

double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   PrimitiveG&  GA, PrimitiveG& GB,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB,
                   vector<double*>& auxd,int n_aux
                  ){
/**
  \brief Compute the matrix element of a 2-nd order pseudopotential with arbitrary primitive Gaussians: <G(A)|PP|G(B)>,
  where PP = [C0 + C2*(r-R)^2]*exp(-alp*(r-R)^2) 

  \param[in] C0 Pseudopotential linear parameter
  \param[in] C2 Pseudopotential linear parameter
  \param[in] alp Pseudopotential exponential parameter
  \param[in] R The center on which the pseudopotential is localized (e.g. could be one of the atomic centers)
  \param[in] GA One primitive Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized)
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each primitive center and w.r.t. the pseudopotential center
  \param[out] dIdR The derivative of the integral w.r.t. the coordinates of the pseudopotential center (if is_derivs = 1)
  \param[out] dIdA The derivative of the integral w.r.t. the coordinates of the first (A) primitive Gaussian (if is_derivs = 1)
  \param[out] dIdB The derivative of the integral w.r.t. the coordinates of the second (B) primitive Gaussian (if is_derivs = 1)
  \param[in,out] auxd The list of the pointers to pre-allocated pieces of memory (for variables of the double type)
  \param[in] n_aux The length of the array to which each of the auxd[i] pointers points.

*/


  double res = 
  libmolint::pseudopot02
  ( C0, C2, alp, R,
    GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dIdR, dIdA, dIdB, auxd, n_aux
  );

  return res;
}

double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   PrimitiveG&  GA, PrimitiveG& GB,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB
                  ){

  double res = 
  libmolint::pseudopot02
  ( C0, C2, alp, R,
    GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dIdR, dIdA, dIdB
  );

  return res;
}


boost::python::list pseudopot02
(double C0, double C2, double alp, const VECTOR& R,
 PrimitiveG&  GA, PrimitiveG& GB,
 int is_normalize, int is_derivs
){

  boost::python::list res;
  res = libmolint::pseudopot02
  ( C0, C2, alp, R,
    GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs
  );

  return res;
}

double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   PrimitiveG&  GA, PrimitiveG& GB,
                   int is_normalize
                  ){

  double res = 
  libmolint::pseudopot02
  ( C0, C2, alp, R,
    GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize
  );

  return res;
}


double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   PrimitiveG&  GA, PrimitiveG& GB
                  ){

  double res = 
  libmolint::pseudopot02
  ( C0, C2, alp, R,
    GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R
  );

  return res;
}




//======================== Multipoles =================================

VECTOR transition_dipole_moment
( PrimitiveG& GA, PrimitiveG& GB,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
){

  VECTOR res;
  res = libmolint::transition_dipole_moment
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dMdA, dMdB, auxd, n_aux
  );

  return res;
}

VECTOR transition_dipole_moment
( PrimitiveG& GA, PrimitiveG& GB,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
){

  VECTOR res;
  res = libmolint::transition_dipole_moment
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dMdA, dMdB
  );

  return res;
}

boost::python::list transition_dipole_moment
( PrimitiveG& GA, PrimitiveG& GB, int is_normalize,int is_derivs){

  boost::python::list res;
  res = libmolint::transition_dipole_moment
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs
  );

  return res;
}

VECTOR transition_dipole_moment
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize){

  VECTOR res;
  res = libmolint::transition_dipole_moment
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize
  );

  return res;
}

VECTOR transition_dipole_moment
( PrimitiveG& GA, PrimitiveG& GB){

  VECTOR res;
  res = libmolint::transition_dipole_moment
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R
  );

  return res;
}


//======================== Derivative coupling integrals =================================

VECTOR derivative_coupling_integral
( PrimitiveG& GA, PrimitiveG& GB,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
){
/**
  \brief Compute the matrix element of the derivative coupling <G(A)|d/dR_B|G(B)>
  note that the d/dR_B uses the derivative w.r.t. the coordinates of the primitive Gaussian B - otherwise it would be zero

  \param[in] GA One primitive Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized)
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each primitive center and w.r.t. the pseudopotential center
  \param[out] dMdA The derivative of the derivative coupling vector w.r.t. the coordinates of the primitive Gaussians A center (if is_derivs = 1)
  \param[out] dMdB The derivative of the derivative coupling vector w.r.t. the coordinates of the primitive Gaussians B center (if is_derivs = 1)
  \param[in,out] auxd The list of the pointers to pre-allocated pieces of memory (for variables of the double type)
  \param[in] n_aux The length of the array to which each of the auxd[i] pointers points.

*/


  VECTOR res;
  res = libmolint::derivative_coupling_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dMdA, dMdB, auxd, n_aux
  );

  return res;
}

VECTOR derivative_coupling_integral
( PrimitiveG& GA, PrimitiveG& GB,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
){

  VECTOR res;
  res = libmolint::derivative_coupling_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dMdA, dMdB
  );

  return res;
}

boost::python::list derivative_coupling_integral
( PrimitiveG& GA, PrimitiveG& GB, int is_normalize,int is_derivs){

  boost::python::list res;
  res = libmolint::derivative_coupling_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs
  );

  return res;
}

VECTOR derivative_coupling_integral
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize){

  VECTOR res;
  res = libmolint::derivative_coupling_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize
  );

  return res;
}

VECTOR derivative_coupling_integral
( PrimitiveG& GA, PrimitiveG& GB){

  VECTOR res;
  res = libmolint::derivative_coupling_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R
  );

  return res;
}



//======================== Kinetic integrals =================================

double kinetic_integral
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
){
/**
  \brief Compute the matrix element of kinetic energy for arbitrary primitive Gaussians: <G(A)|T|G(B)>, where
 -1/2 * ( d^2/dx^2 + d^2/dy^2 + d^2/dz^2 ), here the derivatives are w.r.t. electronic variables, not the centers of the primitives.

  \param[in] GA One primitive Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized)
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each primitive center
  \param[out] dIdA The derivative of the integral w.r.t. the coordinates of the first (A) primitive Gaussian (if is_derivs = 1)
  \param[out] dIdB The derivative of the integral w.r.t. the coordinates of the second (B) primitive Gaussian (if is_derivs = 1)
  \param[in,out] auxd The list of the pointers to pre-allocated pieces of memory (for variables of the double type)
  \param[in] n_aux The length of the array to which each of the auxd[i] pointers points.

*/


  double res = 
  libmolint::kinetic_integral
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dIdA, dIdB, auxd, n_aux
  );

  return res;
}

double kinetic_integral
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB
){

  double res = 
  libmolint::kinetic_integral
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs, dIdA, dIdB
  );

  return res;
}

boost::python::list kinetic_integral
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs){

  boost::python::list res;
  res = 
  libmolint::kinetic_integral
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize, is_derivs
  );

  return res;
}

double kinetic_integral
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize){

  double res = 
  libmolint::kinetic_integral
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    is_normalize
  );

  return res;
}

double kinetic_integral
( PrimitiveG& GA, PrimitiveG& GB){

  double res = 
  libmolint::kinetic_integral
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R
  );

  return res;
}


//======================== NAIs =================================

double nuclear_attraction_integral
( PrimitiveG& GA, PrimitiveG& GB, VECTOR& Rc,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
){
/**
  \brief Compute the matrix element of the nuclear attraction integral for arbitrary primitive Gaussians: <G(A)|V|G(B)>, where
  V(r) = 1/|r-Rc|

  \param[in] GA One primitive Gaussian
  \param[in] GB Second primitive Gaussian
  \param[in] Rc The coordinate of the nucleus with which the electron is interacting 
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized)
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each primitive center
  \param[out] DA The derivative of the integral w.r.t. the coordinates of the first (A) primitive Gaussian (if is_derivs = 1)
  \param[out] DB The derivative of the integral w.r.t. the coordinates of the second (B) primitive Gaussian (if is_derivs = 1)
  \param[out] DC The derivative of the integral w.r.t. the coordinates of the nucleus (C) (if is_derivs = 1)
  \param[in,out] auxd The list of the pointers to pre-allocated pieces of memory (for variables of the double type)
  \param[in] n_aux The length of the array to which each of the auxd[i] pointers points.
  \param[in,out] auxv The list of the pointers to pre-allocated pieces of memory (for variables of the VECTOR type)
  \param[in] n_auxv The length of the array to which each of the auxv[i] pointers points.

*/


  double res = 
  libmolint::nuclear_attraction_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    Rc, is_normalize, is_derivs, DA, DB, DC, aux, n_aux, auxv, n_auxv
  );

  return res;
}

double nuclear_attraction_integral
( PrimitiveG& GA, PrimitiveG& GB, VECTOR& Rc,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC
){

  double res = 
  libmolint::nuclear_attraction_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    Rc, is_normalize, is_derivs, DA, DB, DC  );

  return res;
}

boost::python::list nuclear_attraction_integral
( PrimitiveG& GA, PrimitiveG& GB, VECTOR& Rc,
  int is_normalize,  int is_derivs
){

  boost::python::list res = 
  libmolint::nuclear_attraction_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    Rc, is_normalize, is_derivs );

  return res;
}

double nuclear_attraction_integral
( PrimitiveG& GA, PrimitiveG& GB, VECTOR& Rc,
  int is_normalize
){

  double res = 
  libmolint::nuclear_attraction_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    Rc, is_normalize  );

  return res;
}

double nuclear_attraction_integral
( PrimitiveG& GA, PrimitiveG& GB, VECTOR& Rc
){

  double res = 
  libmolint::nuclear_attraction_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    Rc
  );

  return res;
}




//======================== ERIs =================================

double electron_repulsion_integral
( PrimitiveG& GA, PrimitiveG& GB, PrimitiveG& GC, PrimitiveG& GD,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
){
/**
  \brief Compute the matrix element of the electron repulsion integral for arbitrary primitive Gaussians: <G(A)G(B)|V|G(C)G(D)>, where
  V(r) = 1/|r1-r2|

  \param[in] GA Primitive Gaussian, #1
  \param[in] GB Primitive Gaussian, #2
  \param[in] GC Primitive Gaussian, #3
  \param[in] GD Primitive Gaussian, #4
  \param[in] is_normalize The flag telling whether we need to normalize the result: is_normalize = 1 - 
  the overlap will be normalized (Gaussians are assumed to be unnormalized), is_normalize = 0 - no need for
  normalization (Gaussians are assumed to be normalized)
  \param[in] is_derivs if = 1 - also compute the derivatives of the overlap w.r.t. coordinates of each primitive center
  \param[out] DA The derivative of the integral w.r.t. the coordinates of the first (A) primitive Gaussian (if is_derivs = 1)
  \param[out] DB The derivative of the integral w.r.t. the coordinates of the second (B) primitive Gaussian (if is_derivs = 1)
  \param[out] DC The derivative of the integral w.r.t. the coordinates of the third (C) primitive Gaussian (if is_derivs = 1)
  \param[out] DD The derivative of the integral w.r.t. the coordinates of the fourth (D) primitive Gaussian (if is_derivs = 1)
  \param[in,out] auxd The list of the pointers to pre-allocated pieces of memory (for variables of the double type)
  \param[in] n_aux The length of the array to which each of the auxd[i] pointers points.
  \param[in,out] auxv The list of the pointers to pre-allocated pieces of memory (for variables of the VECTOR type)
  \param[in] n_auxv The length of the array to which each of the auxv[i] pointers points.

*/


  double res = 
  libmolint::electron_repulsion_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    GC.x_exp, GC.y_exp, GC.z_exp, GC.alpha, GC.R,
    GD.x_exp, GD.y_exp, GD.z_exp, GD.alpha, GD.R,
    is_normalize, is_derivs, DA, DB, DC, DD, aux, n_aux, auxv, n_auxv
  );

  return res;
}

double electron_repulsion_integral
( PrimitiveG& GA, PrimitiveG& GB, PrimitiveG& GC, PrimitiveG& GD,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD
){

  double res = 
  libmolint::electron_repulsion_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    GC.x_exp, GC.y_exp, GC.z_exp, GC.alpha, GC.R,
    GD.x_exp, GD.y_exp, GD.z_exp, GD.alpha, GD.R,
    is_normalize, is_derivs, DA, DB, DC, DD  );

  return res;
}

boost::python::list electron_repulsion_integral
( PrimitiveG& GA, PrimitiveG& GB, PrimitiveG& GC, PrimitiveG& GD,
  int is_normalize,  int is_derivs
){

  boost::python::list res = 
  libmolint::electron_repulsion_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    GC.x_exp, GC.y_exp, GC.z_exp, GC.alpha, GC.R,
    GD.x_exp, GD.y_exp, GD.z_exp, GD.alpha, GD.R,
    is_normalize, is_derivs );

  return res;
}

double electron_repulsion_integral
( PrimitiveG& GA, PrimitiveG& GB, PrimitiveG& GC, PrimitiveG& GD,
  int is_normalize
){

  double res = 
  libmolint::electron_repulsion_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    GC.x_exp, GC.y_exp, GC.z_exp, GC.alpha, GC.R,
    GD.x_exp, GD.y_exp, GD.z_exp, GD.alpha, GD.R,
    is_normalize  );

  return res;
}

double electron_repulsion_integral
( PrimitiveG& GA, PrimitiveG& GB, PrimitiveG& GC, PrimitiveG& GD
){

  double res = 
  libmolint::electron_repulsion_integral
  ( GA.x_exp, GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp, GB.y_exp, GB.z_exp, GB.alpha, GB.R,
    GC.x_exp, GC.y_exp, GC.z_exp, GC.alpha, GC.R,
    GD.x_exp, GD.y_exp, GD.z_exp, GD.alpha, GD.R
  );

  return res;
}





}// namespace libqobjects
}// namespace liblibra



