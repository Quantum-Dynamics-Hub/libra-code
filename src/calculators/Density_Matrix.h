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

#include "../math_linalg/liblinalg.h"
#include "../math_meigen/libmeigen.h"
#include "../timer/libtimer.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libmeigen;


/// libcalculators namespace
namespace libcalculators{

void compute_density_matrix(vector< pair<int,double> >& occ, MATRIX* C, MATRIX* P);
void compute_density_matrix(vector< pair<int,double> >& occ, CMATRIX* C, CMATRIX* P);

MATRIX compute_density_matrix(boost::python::list occ, MATRIX C);
CMATRIX compute_density_matrix(boost::python::list occ, CMATRIX C);

void Fock_to_P(int Norb,int Nocc, int degen, double Nel, std::string eigen_method, int pop_opt,
               MATRIX* Fao, MATRIX* Sao, MATRIX* C, MATRIX* E,
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,
               MATRIX* P, vector<Timer>& bench_t);
void Fock_to_P(int Norb,int Nocc, int degen, double Nel, std::string eigen_method, int pop_opt,
               CMATRIX* Fao, CMATRIX* Sao, CMATRIX* C, CMATRIX* E,
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,
               CMATRIX* P, vector<Timer>& bench_t);


// Versions with the optional benchmark
void Fock_to_P(MATRIX* Fao, MATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               MATRIX* E, MATRIX* C, MATRIX* P,                                              /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,           /*Outputs*/
               int BM, vector<Timer>& bench_t);                                              /*Benchmarking data*/
void Fock_to_P(CMATRIX* Fao, CMATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               CMATRIX* E, CMATRIX* C, CMATRIX* P,                                              /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,           /*Outputs*/
               int BM, vector<Timer>& bench_t);                                              /*Benchmarking data*/

// Versions without the benchmark
void Fock_to_P(MATRIX* Fao, MATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               MATRIX* E, MATRIX* C, MATRIX* P,                                                         /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ                       /*Outputs*/
              );
void Fock_to_P(CMATRIX* Fao, CMATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               CMATRIX* E, CMATRIX* C, CMATRIX* P,                                                         /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ                       /*Outputs*/
              );

// Version for the Python
boost::python::list Fock_to_P(MATRIX Fao, MATRIX Sao, double Nel, double degen, double kT, double etol, int pop_opt);
boost::python::list Fock_to_P(CMATRIX Fao, CMATRIX Sao, double Nel, double degen, double kT, double etol, int pop_opt);




}// namespace libcalculators
}// liblibra

#endif // DENSITY_MATRIX_H
