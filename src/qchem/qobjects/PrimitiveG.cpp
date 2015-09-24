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

#include "PrimitiveG.h"
#include "../molint/libmolint.h"

namespace libqchem{

using namespace libmolint;

namespace libqobjects{


//----------------------- Members of PrimitiveG class -----------------------
int PrimitiveG::get_x_exp(){ if(is_x_exp){ return x_exp; } else{ std::cout<<"Error: x_exp is not set\n"; exit(0); } }
int PrimitiveG::get_y_exp(){ if(is_y_exp){ return y_exp; } else{ std::cout<<"Error: y_exp is not set\n"; exit(0); } }
int PrimitiveG::get_z_exp(){ if(is_z_exp){ return z_exp; } else{ std::cout<<"Error: z_exp is not set\n"; exit(0); } }
double PrimitiveG::get_alpha(){ if(is_alpha){ return alpha; } else{ std::cout<<"Error: alpha is not set\n"; exit(0); } }
VECTOR PrimitiveG::get_R(){ if(is_R){ return R; } else{ std::cout<<"Error: R is not set\n"; exit(0); } }
double PrimitiveG::get_value(){ if(is_value){ return value; } else{ std::cout<<"Error: value is not set\n"; exit(0); } }


void PrimitiveG::set_x_exp(int _x){ x_exp = _x; is_x_exp = 1; }
void PrimitiveG::set_y_exp(int _y){ y_exp = _y; is_y_exp = 1; }
void PrimitiveG::set_z_exp(int _z){ z_exp = _z; is_z_exp = 1; }
void PrimitiveG::set_alpha(double _alp){ alpha = _alp; is_alpha = 1; }
void PrimitiveG::set_R(const VECTOR& _R){ R = _R; is_R = 1; }


// General initialization
void PrimitiveG::init(int l,int m,int n,double alp,VECTOR& center){
  x_exp = l;           is_x_exp = 1;
  y_exp = m;           is_y_exp = 1;
  z_exp = n;           is_z_exp = 1;
  alpha = alp;         is_alpha = 1;
  R = center;          is_R = 1;
  value = 0.0;         is_value = 0;
}

// Default constructor
PrimitiveG::PrimitiveG(){  // Default ctor: Simple S-type function
  VECTOR r;
  init(0,0,0,1.0, r );
}

// Arbitrary L function constructor
PrimitiveG::PrimitiveG(int l,int m,int n,double alp,VECTOR& center){ 
  init(l,m,n,alp, center );
}

// Copy constructor
PrimitiveG::PrimitiveG(const PrimitiveG& g){

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

  VECTOR r; r = pos - R;
  double r2 = r.length2();
  value = FAST_POW(r.x,x_exp)*FAST_POW(r.y,y_exp)*FAST_POW(r.z,z_exp)*exp(-alpha * r2);
  is_value = 1;

  return value;
}


double PrimitiveG::norm2(){
/**
  <G(A)|G(A)>
**/
  double res = gaussian_norm2(x_exp,alpha) * gaussian_norm2(y_exp,alpha) * gaussian_norm2(z_exp,alpha);
  return res;

}// norm2

double PrimitiveG::norm1(){
/**
  sqrt(<G(A)|G(A)>)
**/
  return sqrt(norm2()); 

}// norm1

double PrimitiveG::normalization_factor(){
/**
  double nom = pow((8.0*G_alpha),(x_exp+y_exp+z_exp))*FACTORIAL(x_exp)*FACTORIAL(y_exp)*FACTORIAL(z_exp);
  double denom = FACTORIAL(2*x_exp)*FACTORIAL(2*y_exp)*FACTORIAL(2*z_exp);
  res = pow((2.0*G_alpha/M_PI),0.75)*sqrt(nom/denom);

  Meaning: if N - is a result of this function and
  G(l,m,n,alp) = x^l * y^m * z^n * exp(-alp*r^2)  a primitive Gaussian
  then
  g = N * G(l,m,n,alp) - is normaized: integral(g,g) = 1.0

  = 1/sqrt(<G(A)|G(A)>)
**/
  double res = (1.0/sqrt(norm2()));
  return res;

}// norm

void PrimitiveG::shift_position(const VECTOR& dR){  R += dR; }


///=======================================================================================================
///===================== Overload basic functions from libmolint to Gaussian objects  ====================

///======================== Overlaps =================================

double gaussian_overlap
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
){

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

  double res = 
  libmolint::gaussian_overlap
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R
  );

  return res;
}


///======================== Moments =================================

double gaussian_moment
( PrimitiveG& GA, PrimitiveG& G, PrimitiveG& GB,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB, vector<double*>& auxd,int n_aux
){

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

  double res = 
  libmolint::gaussian_moment
  ( GA.x_exp,GA.y_exp, GA.z_exp, GA.alpha, GA.R,
    G.x_exp, G.y_exp,  G.z_exp,  G.alpha,  G.R,
    GB.x_exp,GB.y_exp, GB.z_exp, GB.alpha, GB.R
  );

  return res;
}



///======================== Pseudopotentials =================================

double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   PrimitiveG&  GA, PrimitiveG& GB,
                   int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB,
                   vector<double*>& auxd,int n_aux
                  ){

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




///======================== Multipoles =================================
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


///======================== Derivative coupling integrals =================================
VECTOR derivative_coupling_integral
( PrimitiveG& GA, PrimitiveG& GB,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
){

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



///======================== Kinetic integrals =================================

double kinetic_integral
( PrimitiveG& GA, PrimitiveG& GB,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
){

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


///======================== NAIs =================================

double nuclear_attraction_integral
( PrimitiveG& GA, PrimitiveG& GB, VECTOR& Rc,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
){

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




///======================== ERIs =================================

double electron_repulsion_integral
( PrimitiveG& GA, PrimitiveG& GB, PrimitiveG& GC, PrimitiveG& GD,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
){

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
}// namespace libqchem



