#ifndef BASIS_H
#define BASIS_H

#include <sstream>
#include "PrimitiveG.h"
#include "PrimitiveS.h"
#include "MOAO.h"
#include "Nuclear.h"
#include "Model_Parameters.h"
using namespace std;


// Basis.cpp
void basis_params_s(int, vector<double>&, vector<double>&);
void basis_params_p(int, vector<double>&, vector<double>&);
void basis_params_d(int, vector<double>&, vector<double>&);

void add_basis_ao(std::string, int, VECTOR&, std::string, int, int, double, double, double, double, double, vector<AO>&, int&);

int num_valence_elec(int);
int set_basis_STO_3G_DZ(Nuclear&, Model_Parameters&, vector<AO>&, int&, int&);


// Basis_ovlp.cpp
void update_overlap_matrix(int,int,int,const VECTOR&,const VECTOR&,const VECTOR&,
                           vector<int>&, vector<AO>&,MATRIX*,vector<double*>&, int, Nuclear&);

void update_overlap_matrix_new(int,int,int,const VECTOR&,const VECTOR&,const VECTOR&,
                              vector<int>&, vector<AO>&,MATRIX*,vector<double*>&, int, Nuclear&);



// Basis_map.cpp
void map_atoms_and_orbitals(int, const vector<AO>&, vector<vector<int> >&);
void show_mapping(const vector<vector<int> >&);


#endif // BASIS_H
