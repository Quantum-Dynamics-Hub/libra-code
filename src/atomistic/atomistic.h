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
  \file Hamiltonian_Atomistic.h
  \brief The file describes the atomistic Hamiltonian class (derived from the generic one)
  for calculations of atomistic systems
    
*/

#ifndef ATOMISTIC_H
#define ATOMISTIC_H


#include "../chemobjects/libchemobjects.h"
#include "../calculators/libcalculators.h"
#include "Hamiltonian_Generic/libhamiltonian_generic.h"
#include "Hamiltonian_MM/libhamiltonian_mm.h"
#include "Hamiltonian_QM/libhamiltonian_qm.h"

/// liblibra namespace
namespace liblibra{


using namespace libcalculators;
using namespace libchemobjects;
using namespace libchemobjects::libchemsys;

/// libhamiltonian namespace
namespace libatomistic{

//using namespace liblinalg;
using namespace libhamiltonian_generic;
using namespace libhamiltonian_mm;
using namespace libhamiltonian_qm;


class Hamiltonian_Atomistic : public Hamiltonian{
/**
  The class for atomistic Hamiltonian calculations
*/

  System* _syst;          ///< reference to the system object (which can be allocated internally, or can be an external object)
  vector<int> ham_types;  ///< we can have several types of Hamiltonians. Their status is kept in this array:
                          ///< ham_types[0] - is the status of MM Hamiltonian (0 -not set, 1 - ready to use)
                          ///< ham_types[1] - is the status of QM Hamiltonian (0 -not set, 1 - ready to use)
  
public:

  // Data members:
  listHamiltonian_MM*  mm_ham;   ///< mm part: type = 0 This is the list of classical (MM) Hamiltonians for the sub-systems (e.g. bonds, angles, etc.)
  listHamiltonian_QM*  qm_ham;   ///< qm part: type = 1 This is the list of quantum (QM) Hamiltonians for the sub-systems (e.g. fragments, blocks, etc)


  /// Constructor: only allocates memory and sets up related variables
  Hamiltonian_Atomistic(int, int);

  /// Setups the type of Hamiltonian for this system
  void set_Hamiltonian_type(std::string); 


  //--------- MM Hamiltonians -----------
  // Interaction related functions:
  //int is_new_interaction(Hamiltonian_MM&);
  void show_interactions_statistics();
  void show_interactions();

  void set_atom_types(System& syst,vector<int>& lst,ForceField& ff);
  void set_fragment_types(System& syst, vector<int>& lst,ForceField& ff);

  bool is_active(Atom&,Atom&);
  bool is_active(Atom&,Atom&,Atom&);
  bool is_active(Atom&,Atom&,Atom&,Atom&);

  void set_atom_interactions_for_atoms(System& syst,string int_type,vector<Atom>& top_elt,vector<int>& lst1,vector<int>& lst2,ForceField& ff,int verb);
  void set_group_interactions_for_atoms(System& syst,string int_type,vector<Group>& top_elt,vector<int>& lst1,vector<int>& lst2,ForceField& ff);

  void set_interactions_for_atoms(System& syst, boost::python::list,boost::python::list,ForceField&,int verb, int assign_rings);
  void set_interactions_for_fragments(System& syst, boost::python::list,boost::python::list,ForceField&);

  void apply_pbc_to_interactions(System& syst, int int_type,int nx,int ny,int nz);
  void set_respa_types(std::string inter_type,std::string respa_type);

  MATRIX3x3 get_stress(std::string);


  //--------- QM Hamiltonians -----------   
  void init_qm_Hamiltonian(std::string ctrl_filename);
  void add_excitation(int f_o, int f_s, int t_o, int t_s){ qm_ham->add_excitation(f_o, f_s, t_o, t_s); }
//  void set_excitonic_basis(boost::python::list basis_ex) { qm_ham->set_excitonic_basis(basis_ex); }
  void excite_alp(int I,int J);
  void excite_bet(int I,int J);


  // Destructor
  ~Hamiltonian_Atomistic();

  void set_q(vector<double>& q_);
  void set_v(vector<double>& v_);
  void set_q(boost::python::list q_);
  void set_v(boost::python::list v_);


/*
  // Set properties
  void set_rep(int rep_);

  // Set parameters
  void set_params(vector<double>& params_);
  void set_params(boost::python::list params_);

  void set_q(boost::python::list q_);
  void set_v(vector<double>& v_);
  void set_v(boost::python::list v_);

  // Perform actual computations - this will construct the internals of the object of this type
  void compute();
*/
//  void compute();
  void set_system(System& syst){ _syst = &syst; }  ///< Making the internal _syst variable to point to the external object (binding)
  void compute_diabatic();
  void compute_adiabatic();




};

typedef std::vector<Hamiltonian_Atomistic> Hamiltonian_AtomisticList; ///< data type for keeping a list of atomistic Hamiltonians


}// namespace libhatomistic
}// liblibra

#endif // ATOMISTIC_H
