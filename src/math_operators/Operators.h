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
  \file Operators.h
  \brief The file describes the functions that perform some transformations being the operators found in many applications
    
*/

#ifndef OPERATORS_H
#define OPERATORS_H

#if defined(USING_PCH)
#include "../pch.h"
#else

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// liboperators namespace
namespace liboperators{


void rotate(double&, double&, double);
boost::python::list expt_rotate(double x,double y,double phi);

void shift(double&, double);
double expt_shift(double x,double phi);

void scale(double&, double);
double expt_scale(double x,double phi);


}// namespace liboperators
}// namespace liblibra


#endif // OPERATORS_H
