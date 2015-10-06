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

#ifndef HAMILTONIAN_ATOMISTIC_H
#define HAMILTONIAN_ATOMISTIC_H

#include "../Hamiltonian_Generic/Hamiltonian.h"
//#include "Hamiltonian_MM/Hamiltonian_MM.h"
#include "Hamiltonian_MM/libhamiltonian_mm.h"
#include "Hamiltonian_QM/libhamiltonian_qm.h"

#include "../../chemobjects/libchemobjects.h"
//#include "../../converters/libconverters.h"


using namespace libchemobjects;
using namespace libchemobjects::libchemsys;


namespace libhamiltonian{
namespace libhamiltonian_atomistic{

using namespace libmmath;
using namespace libhamiltonian_generic;
using namespace libhamiltonian_mm;
using namespace libhamiltonian_qm;


class Hamiltonian_Atomistic : public Hamiltonian{

  System* _syst;          // reference to the system object
  vector<int> ham_types;  // we can have several types

  
public:

  // Data members:
  listHamiltonian_MM*  mm_ham;   // mm part: type = 0
  listHamiltonian_QM*  qm_ham;   // qm part: type = 1


  // Constructor: only allocates memory and sets up related variables
  Hamiltonian_Atomistic(int, int);

  // Setups
  void set_Hamiltonian_type(std::string); 


  //--------- MM Hamiltonians -----------
  // Interaction related functions:
  //int is_new_interaction(Hamiltonian_MM&);
  void show_interactions_statistics();

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
  void set_system(System& syst){ _syst = &syst; }  // by reference
  void compute_diabatic();
  void compute_adiabatic();




};

typedef std::vector<Hamiltonian_Atomistic> Hamiltonian_AtomisticList;


}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian

#endif // HAMILTONIAN_ATOMISTIC_H
