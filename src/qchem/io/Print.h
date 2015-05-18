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
