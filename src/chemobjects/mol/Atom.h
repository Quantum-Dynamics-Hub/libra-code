/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Atom.h
  \brief The file describes the Atom class for keeping the atomic information
    
*/

#ifndef ATOM_H
#define ATOM_H

#include "../../dyn_rigidbody/librigidbody.h"
#include "../universe/libuniverse.h"
#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace librigidbody;
using namespace liblinalg;

/// libchemobjects namespace
namespace libchemobjects{

using namespace libuniverse;

/// libmol namespace
namespace libmol{


class Atom{
/**
  \brief The class for keeping atomic information: properties and some dynamic variables 
  Note the dynamic variables here (e.g. coordinates, etc) are still considered as properties (parameters) rather 
  that true (propagated) dynamical variables. This is only a conceptual consideration, in practice these 
  data members can, of course, be used for atomic propagation.
*/

  
  Universe* universe;

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const Atom&); // Copies the content which is defined

public:
  //--------- Atom identification ----------------------------
  int         Atom_id;             int is_Atom_id;   ///< the atomic integer-valued ID - can be any number

  //--------- Topology = Data model -------------------------
  int globAtom_Index;                   ///< the index of this atom in the global array of atoms (indexing starts with 0)
  int locAtom_Index;                    ///< the index of this atom in the local array of atoms - e.g. atoms of a given group (indexing starts with 0)
  int globGroup_Index;                  ///< the index of the group to which the atom belongs. This is the index in the global array of groups (indexing starts with 0)
  int globMolecule_Index;               ///< the index of the molecule to which the atom belongs. This is the index in the global array of molecules (indexing starts with 0)
  vector<int> globAtom_Adjacent_Atoms;  ///< The vector of globar atomic indices for the atoms directly connected (in MM) to the present atom

  //--------- Dynamical properties of the atom --------------
  RigidBody   Atom_RB;             int is_Atom_RB;  ///< this is the object containing the spatial information about the atom (position, velocity, etc, even the orientation, if needed)
  RigidBody   Atom_RB_old;         int is_Atom_RB_old; ///< previous position: from some previous calculations - needed to comparison-based calculations
  double      Atom_displ2;         int is_Atom_displ2; ///< square of atomic displacements (present - old)

  //----------Basic physical properties of the atom ---------
  int         Atom_Z;              int is_Atom_Z;  ///< the nucleus charge of the given atom (e.g. 3 for Li)
  std::string Atom_element;        int is_Atom_element; ///< name of the element of the given atom
  double      Atom_atomic_radius;  int is_Atom_atomic_radius; ///< atomic radius
  double      Atom_charge;         int is_Atom_charge;  ///< The partial charge of the atom - e.g. Mulliken charge
  double      Atom_electronegativity;   int is_Atom_electronegativity; ///< Atomic electronegativity

  //--------- Atom-in-Molecule topological properties --------------
  double      Atom_formal_charge;  int is_Atom_formal_charge; ///< This is the formal charge - e.g. the valence state charge - e.g. +1 on N in protonated amino group
  int         Atom_coordination;   int is_Atom_coordination; ///< defines an equilibrium angle
  std::string Atom_functional_group; int is_Atom_functional_group;///< Name of the functional group to which the atom belongs
  vector<int> Atom_ring_sizes;     int is_Atom_ring_sizes;   ///< Sizes of all rings to which this atom belongs
  int         Atom_min_ring_size;  int is_Atom_min_ring_size;  ///< The size of the minimal ring to which the atom belongs

  //----------- Force-field related properties ---------------
  std::string Atom_ff_type;         int is_Atom_ff_type;  ///< Force-field atom name
  double Atom_Zeff;                 int is_Atom_Zeff;     ///< Effective atomic core charge - e.g. +1 for Li, +2 for Ca, etc.

  //----------- Electronic structure-related properties ---------------
  double Atom_mull_charge_gross;         int is_Atom_mull_charge_gross;  ///< Specifically, the Mulliken gross charge
  double Atom_mull_charge_net;           int is_Atom_mull_charge_net;    ///< Specifically, the Mulliken net charge
 
  //---------- Keep this for a while --------------------
//  int         Atom_ff_int_type;     int is_Atom_ff_int_type;
//  int         Atom_is_surface_atom; int is_Atom_is_surface_atom;
//  int         Atom_surface_index;   int is_Atom_surface_index;
//  int         Atom_is_basis_atom;   int is_Atom_is_basis_atom;
//  int         Atom_is_C60_CT;       int is_Atom_is_C60_CT;


  //--------- Methods -----------------
  Atom(Universe&);
  Atom(Universe&,boost::python::dict);
  Atom(const Atom&); ///< Copy constructor
 ~Atom();            ///< Destructor
  Atom& operator=(const Atom&);
  bool operator==(const Atom& a) const { return globAtom_Index == a.globAtom_Index; }
  bool operator!=(const Atom& a) const { return globAtom_Index != a.globAtom_Index; }
  friend ostream& operator<<(ostream &strm,Atom& ob){ ob.show_info(); return strm; }
  void set(object);
  void show_info();

  void save(boost::property_tree::ptree& pt,std::string path);
  void save(std::string path);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);

};

void save(boost::property_tree::ptree& pt,std::string path,vector<Atom>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Atom>& vt,Universe& u,int& status);


}// namespace libmol
}// namespace libchemobjects
}// liblibra


#endif // ATOM_H
