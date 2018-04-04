/*********************************************************************************
* Copyright (C) 2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Ehrenfest.h
  \brief The file implements all about Ehrenfest calculations
    
*/

#ifndef EHRENFEST_H
#define EHRENFEST_H


// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../hamiltonian/libhamiltonian.h"

// Dynamics classes
#include "nuclear/libnuclear.h"
#include "electronic/libelectronic.h"
#include "ensemble/libensemble.h"

/// liblibra namespace
namespace liblibra{


using namespace libhamiltonian;

/// libdyn namespace
namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
using namespace libensemble;


double Ehrenfest_dia(CMATRIX& C, CMATRIX& H, vector<CMATRIX>& dHdR, vector<double>& f, int opt);
double Ehrenfest_dia(CMATRIX* C, CMATRIX* H, vector<CMATRIX*>& dHdR, vector<double*>& f, int opt);
double Ehrenfest_dia(CMATRIX& C, Hamiltonian& ham, vector<double>& f, int opt);

double Ehrenfest_adi(CMATRIX& C, CMATRIX& E, vector<CMATRIX>& dEdR, vector<CMATRIX>& D, vector<double>& f, int opt);
double Ehrenfest_adi(CMATRIX* C, CMATRIX* E, vector<CMATRIX*>& dEdR, vector<CMATRIX*>& D, vector<double*>& f, int opt);
double Ehrenfest_adi(CMATRIX& C, Hamiltonian& ham, vector<double>& f, int opt);



}// namespace libdyn
}// liblibra

#endif // EHRENFEST
