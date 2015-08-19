#ifndef HAMILTONIAN_INDO_H
#define HAMILTONIAN_INDO_H

#include "../../../qchem/libqchem.h"
using namespace libqchem;

#include "../../../chemobjects/libchemobjects.h"
using namespace libchemobjects;
using namespace libchemobjects::libchemsys;




namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_indo{




class Hamiltonian_INDO{

  public: 
  Hamiltonian_INDO(){ ;; }

};

class listHamiltonian_INDO{

  public: 
  list_Hamiltonian_INDO(){ ;; }

};


// Hamiltonian_INDO.cpp
void indo_core_parameters(vector<int>&, vector<int>&, vector<AO>&, Nuclear&, vector<double>&, vector<double>&, Memory*, int);

void Hamiltonian_core_indo(Control_Parameters&, Model_Parameters&, Nuclear&,
                           vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&, 
                           MATRIX*, MATRIX*, Memory*, vector<double>&, vector<double>&);


void get_integrals(int, int, vector<AO>&, double, double, double, double&, double&);

void Hamiltonian_Fock_indo(Control_Parameters&, Model_Parameters&, Nuclear&,
                           vector<int>&, vector<int>&, vector<AO>&, vector<vector<int> >&,
                           Electronic*, Memory*, vector<double>&, vector<double>&);





}// namespace libhamiltonian_indo
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


#endif // HAMILTONIAN_INDO_H
