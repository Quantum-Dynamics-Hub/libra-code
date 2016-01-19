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
/**
  \file Density_Matrix.h
  \brief The file defines functions for density matrix and Fock-to-density calculations
    
*/

#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

#include "../mmath/libmmath.h"
using namespace libmmath;

/// libcalculators namespace
namespace libcalculators{

void compute_density_matrix(vector< pair<int,double> >& occ, MATRIX* C, MATRIX* P);
MATRIX compute_density_matrix(boost::python::list occ, MATRIX C);

void Fock_to_P(int Norb,int Nocc, int degen, double Nel, std::string eigen_method, int pop_opt,
               MATRIX* Fao, MATRIX* Sao, MATRIX* C, MATRIX* E,
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,
               MATRIX* P, vector<Timer>& bench_t);

void Fock_to_P(MATRIX* Fao, MATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               MATRIX* E, MATRIX* C, MATRIX* P,                                              /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,           /*Outputs*/
               int BM, vector<Timer>& bench_t);                                              /*Benchmarking data*/
void Fock_to_P(MATRIX* Fao, MATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               MATRIX* E, MATRIX* C, MATRIX* P,                                                         /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ                       /*Outputs*/
              );
boost::python::list Fock_to_P(MATRIX Fao, MATRIX Sao, double Nel, double degen, double kT, double etol, int pop_opt);




}// namespace libcalculators

#endif // DENSITY_MATRIX_H
