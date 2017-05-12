/*********************************************************************************
* Copyright (C) 2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
* This file is originally a part of the ErgoSCF code (see the info below). 
* Now, it is meant to use in the Libra code.
* WOW!: This is the FIRST file to integrate ErgoSCF and Libra (on May 11, 2017)
*
*********************************************************************************/


/* Ergo, version 3.3, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2013 Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */

#ifndef REALTYPEHEADER
#define REALTYPEHEADER

/// We incorporate all the definitions into the liblibra namespace
/// liblibra namespace
namespace liblibra{

/// Also, to keep everything organized, lets also create a libergoescf namespace
/// to keep all the original developments in there
/// libergoescf namespace 
namespace libergoscf{


/// AVA: We don't want that for now. Mostly, because
/// I have no idea on how to use it (seems that all the
/// definitions can be done in the CMakeLists.txt file)
//#include "config.h"

#ifdef PRECISION_SINGLE
typedef float ergo_real;
typedef double ergo_long_real;
#define REALTYPE_DEFINED_OK
#endif

#ifdef PRECISION_DOUBLE
typedef double ergo_real;
typedef double ergo_long_real;
#define REALTYPE_DEFINED_OK
#endif

#ifdef PRECISION_LONG_DOUBLE
typedef long double ergo_real;
typedef long double ergo_long_real;
#define REALTYPE_DEFINED_OK
#endif

/* if precision not specified, use double as default */
#ifndef REALTYPE_DEFINED_OK
typedef double ergo_real;
typedef double ergo_long_real;
#endif


}// libergoescf
}// liblibra


#endif
