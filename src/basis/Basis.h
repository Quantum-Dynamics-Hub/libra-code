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
  \file Basis.h
  \brief The file describes functions for creating basis atomic orbitals from molecular structure information
    
*/

#ifndef BASIS_H
#define BASIS_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <sstream>
#endif

#include "../math_linalg/liblinalg.h"
#include "../qobjects/libqobjects.h"


/// liblibra namespace
namespace liblibra{

using namespace std;
using namespace liblinalg;
using namespace libqobjects;

/// libbasis namespace
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

void MO_overlap(MATRIX& Smo, vector<AO>& ao_i, vector<AO>& ao_j, MATRIX& Ci, MATRIX& Cj,
 vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2);

void MO_overlap(CMATRIX& Smo, vector<AO>& ao_i, vector<AO>& ao_j, CMATRIX& Ci, CMATRIX& Cj,
 vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2);

void MO_overlap(MATRIX& Smo, MATRIX& Ci, MATRIX& Cj, 
 vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2);

void MO_overlap(CMATRIX& Smo, CMATRIX& Ci, CMATRIX& Cj,
 vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2);


complex<double> SD_overlap(SD& sd_i, SD& sd_j);

CMATRIX SD_overlap(vector<SD>& sd_i, vector<SD>& sd_j);

void SD_overlap(CMATRIX& SD_ovlp, vector<SD>& sd_i, vector<SD>& sd_j);



// Basis_map.cpp
void show_mapping(const vector<vector<int> >&);


// Basis_nac.cpp
void update_derivative_coupling_matrix
(int x_period,int y_period,int z_period,const VECTOR& t1, const VECTOR& t2, const VECTOR& t3,
 vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
 vector<AO>& basis_ao, int c, MATRIX& Dao_x, MATRIX& Dao_y, MATRIX& Dao_z
);




}//namespace libbasis
}//namespace liblibra

#endif // BASIS_H


