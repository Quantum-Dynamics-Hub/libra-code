#include "Operators.h"

namespace libmmath{
namespace liboperators{

void rotate(double& x,double& y,double phi){
// Application of the operator 
// exp(iL * phi) = exp (phi * (x * d/dy - y * d/dx) ) to vector (x,y)^T:
//
//               | x |     |  cos(phi)   -sin(phi)  |    | x |
// exp(iL * phi) |   |  =  |                        | *  |   |
//               | y |     |  sin(phi)    cos(phi)  |    | y |
//
  double c = std::cos(phi);
  double s = std::sin(phi);

  double tmpx =  c*x - s*y;
  double tmpy =  s*x + c*y;

  x = tmpx;
  y = tmpy;

}

boost::python::list expt_rotate(double x,double y,double phi){
  boost::python::list res;
  double rx,ry; rx = x; ry = y;

  rotate(rx,ry,phi);

  res.append(rx);
  res.append(ry);
 
  return res;
}



void shift(double& x, double phi){
// Application of the operator 
// exp(iL * phi) = exp (phi * d/dx ) to vector variable x:
//
// exp(iL * phi) x = x + phi

  x += phi;

}

double expt_shift(double x,double phi){

  boost::python::list res;
  double rx; rx = x;

  shift(rx,phi);
 
  return rx;
}



void scale(double& x, double phi){
// Application of the operator 
// exp(iL * phi) = exp (phi * x*d/dx ) to vector variable x:
//
// exp(iL * phi) x = exp(phi) * x

  x *= std::exp(phi);

}

double expt_scale(double x,double phi){

  double rx; rx = x;

  scale(rx,phi);

  return rx;
}



}// namespace liboperators
}// namespace libmmath

