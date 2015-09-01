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

#ifndef BASIS_H
#define BASIS_H

#include <sstream>
using namespace std;

#include "../../mmath/libmmath.h"
using namespace libmmath;

#include "../qobjects/libqobjects.h"
using namespace libqchem::libqobjects;



namespace libqchem{
namespace libbasis{



// Basis.cpp
void basis_params_s(int, vector<double>&, vector<double>&);
void basis_params_p(int, vector<double>&, vector<double>&);
void basis_params_d(int, vector<double>&, vector<double>&);

void add_basis_ao(std::string Atom_name, VECTOR& R, std::string Atom_shell, int Nzeta, int Nquant,
                  double  IP, double exp1, double exp2, double coeff1, double coeff2,
                  vector<AO>& basis_ao);

void add_basis_ao(std::string Atom_name, VECTOR& R, std::string Atom_shell, int Nzeta, int Nquant,
                  double  IP, double exp1, double exp2, double coeff1, double coeff2,
                  boost::python::list basis_ao);



int num_valence_elec(int);


// Basis_ovlp.cpp
void update_overlap_matrix(int,int,int,const VECTOR&,const VECTOR&,const VECTOR&, vector<AO>&,MATRIX&);

// Basis_map.cpp
void show_mapping(const vector<vector<int> >&);



}//namespace libbasis
}//namespace libqchem

#endif // BASIS_H


