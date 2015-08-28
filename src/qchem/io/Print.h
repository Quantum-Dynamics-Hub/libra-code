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

#ifndef PRINT_H
#define PRINT_H

#include "Nuclear.h"
#include "Electronic.h"

void print_xyz(std::string, Nuclear&);
void print_xyz(std::string, Nuclear&, double, double, double);
void print_xyzq(std::string, Nuclear&);
void print_xyzq(std::string, Nuclear&, double, double, double);
void print_el_struct(std::string, Electronic*,double,double);
void print_Mulliken_charges(std::string,Nuclear&);
void print_dipole(std::string,Electronic*,double,MATRIX*, MATRIX*, MATRIX*,VECTOR&);
void print_excitations(std::string filename, Electronic*, Control_Parameters& prms);

#endif // PRINT_H
