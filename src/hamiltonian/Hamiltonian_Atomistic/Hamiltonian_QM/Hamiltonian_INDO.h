#ifndef HAMILTONIAN_INDO_H
#define HAMILTONIAN_INDO_H

#include "../../../calculators/libcalculators.h"
using namespace libcalculators;

#include "../../../qchem/libqchem.h"
using namespace libqchem;
using namespace libqchem::libqobjects;

#include "../../../chemobjects/libchemobjects.h"
using namespace libchemobjects;
using namespace libchemobjects::libchemsys;

#include "Control_Parameters/libcontrol_parameters.h"
using namespace libhamiltonian::libhamiltonian_atomistic::libhamiltonian_qm::libcontrol_parameters;

#include "Model_Parameters/libmodel_parameters.h"
using namespace libhamiltonian::libhamiltonian_atomistic::libhamiltonian_qm::libmodel_parameters;

#include "Electronic.h"


namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_qm{



// Hamiltonian_INDO.cpp
void indo_core_parameters
( System& syst, vector<AO>& basis_ao, 
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  vector<double>& eri, vector<double>& V_AB, int opt);

/*
void indo_core_parameters(vector<int>&, vector<int>&, vector<AO>&, Nuclear&, vector<double>&, vector<double>&, Memory*, int);

void Hamiltonian_core_indo(Control_Parameters&, Model_Parameters&, Nuclear&,
                           vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&, 
                           MATRIX*, MATRIX*, Memory*, vector<double>&, vector<double>&);


void get_integrals(int, int, vector<AO>&, double, double, double, double&, double&);

void Hamiltonian_Fock_indo(Control_Parameters&, Model_Parameters&, Nuclear&,
                           vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&,
                           Electronic*, Memory*, vector<double>&, vector<double>&);
*/




}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


#endif // HAMILTONIAN_INDO_H
