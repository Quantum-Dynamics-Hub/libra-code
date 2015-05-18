#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "Engine.h"
#include "Control_Parameters.h"
#include "Model_Parameters.h"
#include "Nuclear.h"
#include "Electronic.h"
#include "MOAO.h"
#include "Memory.h"

// Hamiltonian.cpp
void Hamiltonian_core(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,MATRIX*,MATRIX*,Memory*);
void Hamiltonian_Fock(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,Electronic*,Electronic*,Memory*);


// Hamiltonian_EHT.cpp
void Hamiltonian_core_eht(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,MATRIX*,MATRIX*,Memory*);
void Hamiltonian_Fock_eht(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,Electronic*,Electronic*,Memory*);

// Hamiltonian_GEHT.cpp
void Hamiltonian_core_geht(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,MATRIX*,MATRIX*,Memory*);
void Hamiltonian_Fock_geht(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,Electronic*,Electronic*,Memory*);


// Hamiltonian_HF.cpp
void Hamiltonian_core_hf(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,MATRIX*,MATRIX*,Memory*);
void Hamiltonian_Fock_hf(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,Electronic*,Electronic*,Memory*);


// Hamiltonian_INDO.cpp
void indo_core_parameters(vector<int>&, vector<int>&, vector<AO>&, Nuclear&, vector<double>&, vector<double>&, Memory*, int);

void Hamiltonian_core_indo(Control_Parameters&, Model_Parameters&, Nuclear&,
                           vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&, 
                           MATRIX*, MATRIX*, Memory*, vector<double>&, vector<double>&);


void get_integrals(int, int, vector<AO>&, double, double, double, double&, double&);

void Hamiltonian_Fock_indo(Control_Parameters&, Model_Parameters&, Nuclear&,
                           vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&,
                           Electronic*, Memory*, vector<double>&, vector<double>&);


// Hamiltonian_GEHT1.cpp
void Hamiltonian_core_geht1(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,MATRIX*,MATRIX*,Memory*);
void Hamiltonian_Fock_geht1_v1(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,Electronic*,Electronic*,Memory*);
void Hamiltonian_Fock_geht1(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,Electronic*,Electronic*,Memory*);

// Hamiltonian_GEHT1.cpp
void geht1_core_parameters(vector<int>&, vector<int>&, vector<AO>&, Nuclear&, vector<double>&, vector<double>&, Memory*, int);

//void Hamiltonian_core_geht1(Control_Parameters&, Model_Parameters&, Nuclear&,
//                            vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&, 
//                            MATRIX*, MATRIX*, Memory*, vector<double>&, vector<double>&);


void get_geht1_integrals(int, int, vector<AO>&, double, double, double, double&, double&);

//void Hamiltonian_Fock_geht1(Control_Parameters&, Model_Parameters&, Nuclear&,
//                            vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&,
//                            Electronic*, Memory*, vector<double>&, vector<double>&);


// Hamiltonian_GEHT2.cpp
void Hamiltonian_core_geht2(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,MATRIX*,MATRIX*,Memory*);
void Hamiltonian_Fock_geht2(Control_Parameters&,Model_Parameters&,Nuclear&,vector<int>&,vector<int>&,vector<AO>&,vector<vector<int> >&,Electronic*,Electronic*,Memory*);




#endif // HAMILTONIAN_H
