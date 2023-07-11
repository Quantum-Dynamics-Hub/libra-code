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
/**
  \file Operators.cpp
  \brief The file implements the functions that perform some transformations being the operators found in many applications
    
*/

#include "Operators.h"

/// liblibra namespace
namespace liblibra{

/// liboperators namespace
namespace liboperators{

void rotate(double& x,double& y,double phi){
/** 
  This operator performs in-plane rotation - changes 2 components of an arbitrary vector
  
  \param[in,out] x The x component of 2D vector
  \param[in,out] y The y component of 2D vector
  \param[in] phi The angle of rotation (radians)

  Application of the operator 
  exp(iL * phi) = exp (phi * (x * d/dy - y * d/dx) ) to vector (x,y)^T:

                | x |     |  cos(phi)   -sin(phi)  |    | x |
  exp(iL * phi) |   |  =  |                        | *  |   |
                | y |     |  sin(phi)    cos(phi)  |    | y |
 
*/

  double c = std::cos(phi);
  double s = std::sin(phi);

  double tmpx =  c*x - s*y;
  double tmpy =  s*x + c*y;

  x = tmpx;
  y = tmpy;

}

boost::python::list expt_rotate(double x,double y,double phi){
/** 
  This is the Python version of the void rotate(double& x,double& y,double phi) function.
  The difference is that the input variables x and y are not changed. The new variables
  are returned as a 2-element Python list

  \param[in] x The x component of 2D vector
  \param[in] y The y component of 2D vector
  \param[in] phi The angle of rotation (radians)

*/
  boost::python::list res;
  double rx,ry; rx = x; ry = y;

  rotate(rx,ry,phi);

  res.append(rx);
  res.append(ry);
 
  return res;
}



void shift(double& x, double phi){
/** 
  This is the shift operator: 

  exp(iL * phi) = exp (phi * d/dx ) to vector variable x:
 
  exp(iL * phi) x = x + phi

  \param[in,out] x The variable on which the operator acts
  \param[in] phi The shift magnitude

*/

  x += phi;

}

double expt_shift(double x,double phi){
/** 
  This is the Python version of the void shift(double& x, double phi) function.
  The difference is that the input variable x is not changed. The new variable
  is returned as a double value

  \param[in] x The variable on which the operator acts
  \param[in] phi The shift magnitude

*/

  boost::python::list res;
  double rx; rx = x;

  shift(rx,phi);
 
  return rx;
}



void scale(double& x, double phi){
/** 
  Application of the operator 
  exp(iL * phi) = exp (phi * x*d/dx ) to vector variable x:

  exp(iL * phi) x = exp(phi) * x

  \param[in,out] x The variable on which the operator acts
  \param[in] phi The scaling factor

*/

  x *= std::exp(phi);

}

double expt_scale(double x,double phi){
/** 
  This is the Python version of the void scale(double& x, double phi) function.
  The difference is that the input variable x is not changed. The new variable
  is returned as a double value

  \param[in] x The variable on which the operator acts
  \param[in] phi The scaling factor

*/


  double rx; rx = x;

  scale(rx,phi);

  return rx;
}



}// namespace liboperators
}// namespace liblibra

