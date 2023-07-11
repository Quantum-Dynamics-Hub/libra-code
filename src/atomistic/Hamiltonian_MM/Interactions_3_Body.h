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
  \file Interactions_3_Body.h
  \brief The file describes functions and classes for molecular-mechanical Hamiltonian calculations 
*/

#ifndef INTERACTIONS_3_BODY_H
#define INTERACTIONS_3_BODY_H

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



class Interaction_3_Body : public Interaction_N_Body{

public:


    // Constructor & Destructor
    Interaction_3_Body() : Interaction_N_Body(3){ int_type = 3; };

    void set_coords(VECTOR* r1_, VECTOR* r2_, VECTOR* r3_);
    void set_coords(VECTOR& r1_, VECTOR& r2_, VECTOR& r3_);

    void set_transl(VECTOR* t1_, VECTOR* t2_, VECTOR* t3_);
    void set_transl(VECTOR& t1_, VECTOR& t2_, VECTOR& t3_);

    void set_forces(VECTOR* f1_, VECTOR* f2_, VECTOR* f3_);
    void set_forces(VECTOR& f1_, VECTOR& f2_, VECTOR& f3_);

    void set_charges(double* q1_, double* q2_, double* q3_);
    void set_charges(double& q1_, double& q2_, double& q3_);


};


class Angle_Interaction : public Interaction_3_Body{

public:

    double k_theta;
    double theta_0;
    double cos_theta_0;
    double C0,C1,C2;
    int coordination;


    // Constructor & Destructor
    Angle_Interaction() : Interaction_3_Body(){ int_type = 30; };

    void set_functional(std::string f);
    void set_params(map<std::string,double>& params);

    void compute();

};






}// namespace libhamiltonian_mm
}// namespace libatomistic
}// liblibra


#endif // INTERACTIONS_3_BODY_H
