/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Hamiltonian_QM.h
  \brief The file describes functions for quantum-mechanical Hamiltonian calculations as well
  as the classes for organizing such computations in an object-oriented way. 
*/

#ifndef HAMILTONIAN_QM_H
#define HAMILTONIAN_QM_H

#include "../../basis_setups/libbasis_setups.h"
#include "../../control_parameters/libcontrol_parameters.h"
#include "../../common_types/libcommon_types.h"
#include "../../model_parameters/libmodel_parameters.h"


#include "Hamiltonian_EHT.h"
#include "Hamiltonian_INDO.h"
#include "Hamiltonian_HF.h"
#include "Electronic_Structure.h"

/// liblibra namespace
namespace liblibra{

using namespace libbasis_setups;
using namespace libcontrol_parameters;
using namespace libmodel_parameters;
using namespace libcommon_types;


namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{


// Hamiltonian_QM.cpp
void Hamiltonian_core(
  System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF
);

void Hamiltonian_Fock(
  Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
);

void Hamiltonian_Fock(
  Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
);



double energy_and_forces
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
);


void derivative_couplings
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
  int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3,
  MATRIX& Dmo_a_x, MATRIX& Dmo_a_y, MATRIX& Dmo_a_z,
  MATRIX& Dmo_b_x, MATRIX& Dmo_b_y, MATRIX& Dmo_b_z
);

void derivative_couplings1
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
  int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3,
  MATRIX& Dmo_a_x, MATRIX& Dmo_a_y, MATRIX& Dmo_a_z,
  MATRIX& Dmo_b_x, MATRIX& Dmo_b_y, MATRIX& Dmo_b_z,
  MATRIX& dEa_dx,  MATRIX& dEa_dy,  MATRIX& dEa_dz,
  MATRIX& dEb_dx,  MATRIX& dEb_dy,  MATRIX& dEb_dz
);

void derivative_couplings1
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
  int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3,
  MATRIX& Dmo_a_x, MATRIX& Dmo_a_y, MATRIX& Dmo_a_z,
  MATRIX& Dmo_b_x, MATRIX& Dmo_b_y, MATRIX& Dmo_b_z
);




VECTOR force
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int Norb, int at_indx, 
  int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3
);



// Excitations.cpp
void excite(int Norb, excitation& ex, 
            int Nocc_alp, vector< pair<int,double> >& occ_alp,
            int Nocc_bet, vector< pair<int,double> >& occ_bet);




class listHamiltonian_QM{
/**
  This class contains everything that is needed for quantum calculations on a 
  specified system. For this system we define its own charge, basis, nuclei, and electronic
  structure variables, as well as an excitation. Consider this class encoding 
  a fragment (sub-system). Such fragments can constitute a basis for large-scale calculations
*/

public:

    // Constructors
    listHamiltonian_QM(){ add_excitation(0,1,0,1); }
    listHamiltonian_QM(std::string ctrl_filename,System& syst);
    listHamiltonian_QM(const listHamiltonian_QM&);   ///< Copy constructor;
    ~listHamiltonian_QM();

    void operator=(const listHamiltonian_QM&);       ///< Copying one listHamiltonian_QM into the other one


    int Nelec;  ///< The number of electrons in this sub-system
    int Norb;   ///< The total number of orbitals in this sub-system
    Control_Parameters prms;  ///< Control parameters defining how to perform calculations for this sub-system
                              ///< Note: different levels of treatment are possible for different sub-systems
    Model_Parameters modprms; ///< Model parameters for this sub-system

    //============= Ground state =================
    std::vector<AO> basis_ao;       ///< Basis for this sub-system
    std::vector<vector<int> > atom_to_ao_map; ///< Nuclei --> AO mapping for this sub-system
    std::vector<int> ao_to_atom_map;          ///< AO --> Nuclei mapping for this sub-system
    Electronic_Structure* el;     ///< Electronic structure for this sub-system

    //============ Excited states ================
    vector<excitation> basis_ex;  ///< Excitations for this sub-system -  may be the same as in prms, but may be different

    

    // Methods
    void init(std::string ctrl_filename,System& syst);
    double compute_scf(System& syst);

    void get_parameters_from_file(std::string filename){ libcontrol_parameters::get_parameters_from_file(filename, prms); }
    Electronic_Structure get_electronic_structure(){ return *el; }
    void set_electronic_structure(Electronic_Structure& el_);

    void compute_overlap(System& syst);
    void compute_core_Hamiltonian(System& syst);
    double energy_and_forces(System& syst);

    void add_excitation(int f_o, int f_s, int t_o, int t_s);

    //void set_excitonic_basis(boost::python::list basis_ex);
    void excite_alp(int I,int J);
    void excite_bet(int I,int J);


/*
    friend int operator != (const listHamiltonian_QM& v1, const listHamiltonian_QM& v2){
      return (v1.Nelec!=v2.Nelec||v1.Norb!=v2.Norb);
    }
    friend int operator == (const listHamiltonian_QM& v1, const listHamiltonian_QM& v2){
      return !(v1!=v2);
    }
*/

};

typedef std::vector<listHamiltonian_QM> listHamiltonian_QMList;  ///< Data type that holds a vector of listHamiltonian_QM objects
typedef std::vector<vector<listHamiltonian_QM> > listHamiltonian_QMMap;  ///< Data type that holds the table (grid) of listHamiltonian_QM objects






}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra

#endif // HAMILTONIAN_QM_H
