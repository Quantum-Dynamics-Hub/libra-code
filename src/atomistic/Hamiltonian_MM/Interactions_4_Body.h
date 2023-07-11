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
  \file Interactions_4_Body.h
  \brief The file describes functions and classes for molecular-mechanical Hamiltonian calculations 
*/

#ifndef INTERACTIONS_4_BODY_H
#define INTERACTIONS_4_BODY_H

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




class Interaction_4_Body : public Interaction_N_Body{

public:


    // Constructor & Destructor
    Interaction_4_Body() :  Interaction_N_Body(4){ int_type = 4; };

    void set_coords(VECTOR* r1_, VECTOR* r2_, VECTOR* r3_, VECTOR* r4_);
    void set_coords(VECTOR& r1_, VECTOR& r2_, VECTOR& r3_, VECTOR& r4_);

    void set_transl(VECTOR* t1_, VECTOR* t2_, VECTOR* t3_, VECTOR* t4_);
    void set_transl(VECTOR& t1_, VECTOR& t2_, VECTOR& t3_, VECTOR& t4_);

    void set_forces(VECTOR* f1_, VECTOR* f2_, VECTOR* f3_, VECTOR* f4_);
    void set_forces(VECTOR& f1_, VECTOR& f2_, VECTOR& f3_, VECTOR& f4_);

    void set_charges(double* q1_, double* q2_, double* q3_, double* q4_);
    void set_charges(double& q1_, double& q2_, double& q3_, double& q4_);

};


class Dihedral_Interaction : public Interaction_4_Body{

public:

    double Vphi,phi0;
    double Vphi1,Vphi2,Vphi3;
    int opt,n;

    // Constructor & Destructor
    Dihedral_Interaction() : Interaction_4_Body(){ int_type = 40; };

    void set_functional(std::string f);
    void set_params(map<std::string,double>& params);

    void compute();
};



class OOP_Interaction : public Interaction_4_Body{

public:


    double K,C0,C1,C2,xi_0;
    int opt;

    // Constructor & Destructor
    OOP_Interaction() : Interaction_4_Body(){ int_type = 41; };

    void set_functional(std::string f);
    void set_params(map<std::string,double>& params);

    void compute();
};










}// namespace libhamiltonian_mm
}// namespace libatomistic
}// liblibra


#endif // INTERACTIONS_4_BODY_H
