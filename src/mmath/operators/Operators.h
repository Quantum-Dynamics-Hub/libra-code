#ifndef OPERATORS_H
#define OPERATORS_H

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;


namespace libmmath{
namespace liboperators{


void rotate(double&, double&, double);
boost::python::list expt_rotate(double x,double y,double phi);

void shift(double&, double);
double expt_shift(double x,double phi);

void scale(double&, double);
double expt_scale(double x,double phi);


}// namespace liboperators
}// namespace libmmath


#endif // OPERATORS_H
