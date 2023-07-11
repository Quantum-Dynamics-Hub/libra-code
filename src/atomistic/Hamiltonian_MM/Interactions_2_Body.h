/*********************************************************************************
* Copyright (C) 2017-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Interactions.h
  \brief The file describes functions and classes for molecular-mechanical Hamiltonian calculations 
*/

#ifndef INTERACTIONS_2_BODY_H
#define INTERACTIONS_2_BODY_H

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "Interactions.h"



/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_mm namespace
namespace libhamiltonian_mm{

using namespace libchemobjects;
using namespace libchemobjects::libchemsys;
using namespace libpot;
using namespace libcell;
using namespace libforcefield;
using namespace boost::python;



class Interaction_2_Body : public Interaction_N_Body{

public:

    // Constructor & Destructor
    Interaction_2_Body() : Interaction_N_Body(2) { int_type = 2; };

    void set_coords(VECTOR* r1_, VECTOR* r2_);
    void set_coords(VECTOR& r1_, VECTOR& r2_);

    void set_transl(VECTOR* t1_, VECTOR* t2_);
    void set_transl(VECTOR& t1_, VECTOR& t2_);

    void set_forces(VECTOR* f1_, VECTOR* f2_);
    void set_forces(VECTOR& f1_, VECTOR& f2_);

    void set_charges(double* q1_, double* q2_);
    void set_charges(double& q1_, double& q2_);


};



class Bond_Interaction : public Interaction_2_Body{
/**
  Functional type        Possible params keys      Meaning
   bond                    K                      Harmonic force constant
                           D                      Morse potential well depth
                           r0                     Equilibrium bond length
                           alpha                  Morse alpha parameter (exponent)

*/

public:

    double K,D,r0,alpha;

    Bond_Interaction() : Interaction_2_Body() { int_type = 20; };

    void set_functional(std::string f);
    void set_params(map<std::string,double>& params);
    void set_params(boost::python::dict);

    void compute();

};

class VdW_Interaction : public Interaction_2_Body{
/**
  vdw                      sigma                  Atomic vdw radius
                           epsilon                vdW binding energy in a pure substance
                           D                      Morse potential well depth - for vdw interactions
                           r0                     Equilibrium non-bonded length - for vdw interactions
                           alpha                  Morse alpha parameter (exponent) - for vdw interactions
                           scale                  The scaling factor for this specific interaction - usually 0 for bonded pairs and even for 1,3-pairs

*/

public:

    double sigma,epsilon,D,r0,alpha,scale;

    VdW_Interaction() : Interaction_2_Body() { int_type = 21; };

    void set_functional(std::string f);
    void set_params(map<std::string,double>& params);
    void set_params(boost::python::dict);

    void compute(double Ron, double Roff);
    void compute();
};


class Elec_Interaction : public Interaction_2_Body{

public:

    double J,xi1,xi2,eps,delta,scale;   

    Elec_Interaction() : Interaction_2_Body() { int_type = 22; };

    void set_functional(std::string f);
    void set_params(map<std::string,double>& params);
    void set_params(boost::python::dict);

    void compute(double Ron, double Roff);
    void compute();

};








}// namespace libhamiltonian_mm
}// namespace libatomistic
}// liblibra


#endif // INTERACTIONS_2_BODY_H
