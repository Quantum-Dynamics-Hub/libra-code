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

#ifndef GROUP_H
#define GROUP_H

#include "../../dyn_rigidbody/librigidbody.h"
#include "../universe/libuniverse.h"
#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace librigidbody;
using namespace liblinalg;


namespace libchemobjects{
namespace libmol{



class Group{

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const Group&); // Copies the content which is defined

public:

//--------- Topology = Data Model -------------------
  int  globGroup_Size;
  int  locGroup_Size;
  int  Group_Size;
  vector<int> globAtom_Index;
  vector<int> locAtom_Index;
  int globGroup_Index;
  int locGroup_Index;
  int globMolecule_Index;

//--------- Properties and their Descriptors --------------

// General properties

  std::string Group_name;         int is_Group_name;
  int         Group_id;           int is_Group_id;
//  double      Group_mass;         int is_Group_mass;
  double      Group_radius;       int is_Group_radius;
  RigidBody   Group_RB;           int is_Group_RB;

  //----------- Force-field related properties ---------------
  std::string Group_ff_type;      int is_Group_ff_type;


//~~~~~~~~~~~~~~~~ Force field related properties ~~~~~~~~~~~~~~~~~~

  // For bond this is length in Angstroems, for angles and dihedrals - angle in radians
  // for pairs - vdw this is sigma_ij
  // Group_equilibrium_geometry and Group_force_constant 1 and 2 - are additional variables
  // for storing informations for such constructs as stretch-bend and torsion angle Fourier expansion

//  double      Group_equilibrium_geometry;    int is_Group_equilibrium_geometry;
//  double      Group_equilibrium_geometry1;   int is_Group_equilibrium_geometry1;
//  double      Group_equilibrium_geometry2;   int is_Group_equilibrium_geometry2;

//  double      Group_current_geometry;        int is_Group_current_geometry;

//  double      Group_force_constant;          int is_Group_force_constant;
//  double      Group_force_constant1;         int is_Group_force_constant1;
//  double      Group_force_constant2;         int is_Group_force_constant2;

//  double      Group_force_constant_half;     int is_Group_force_constant_half;


// Group-specific properties
// Bonds
  double      Group_bond_order;              int is_Group_bond_order;
  double      Group_bond_alpha;              int is_Group_bond_alpha;

// Angles
//  double      Group_angle_Cijk;              int is_Group_angle_Cijk;
//  double      Group_angle_C0;                int is_Group_angle_C0;
//  double      Group_angle_C1;                int is_Group_angle_C1;
//  double      Group_angle_C2;                int is_Group_angle_C2;
//  double      Group_angle_coordination;      int is_Group_angle_coordination;

// Dihedrals
//  double      Group_dihedral_multiplicity;   int is_Group_dihedral_multiplicity;
//  double      Group_dihedral_periodicity;    int is_Group_dihedral_periodicity;

// Pairs
//  double      Group_pair_A_vdw;              int is_Group_pair_A_vdw;
//  double      Group_pair_B_vdw;              int is_Group_pair_B_vdw;
//  double      Group_pair_sigma_vdw;          int is_Group_pair_sigma_vdw;
//  double      Group_pair_epsil_vdw;          int is_Group_pair_epsil_vdw;
//  double      Group_pair_sigma7_vdw;         int is_Group_pair_sigma7_vdw;
//  double      Group_pair_alpha_vdw;          int is_Group_pair_alpha_vdw;
//  double      Group_pair_q1q2;               int is_Group_pair_q1q2;
//  double      Group_pair_shift_elec;         int is_Group_pair_shift_elec;
//  double      Group_pair_dielectric;         int is_Group_pair_dielectric;
  // Parameters for charge transfer force field (MALINA)
//  double      Group_pair_wij;                int is_Group_pair_wij;
//  double      Group_pair_wij_1;              int is_Group_pair_wij_1;
//  double      Group_pair_wij_2;              int is_Group_pair_wij_2;
//  double      Group_pair_alpij;              int is_Group_pair_alpij;
//  double      Group_pair_alpij_1;            int is_Group_pair_alpij_1;
//  double      Group_pair_alpij_2;            int is_Group_pair_alpij_2;
//  double      Group_pair_rij_1;              int is_Group_pair_rij_1;
//  double      Group_pair_rij_2;              int is_Group_pair_rij_2;

//  int         Group_pair_excluded;           int is_Group_pair_excluded;

//--------- Methods -----------------
  Group();
  Group(const Group&); // Copy constructor
  Group& operator=(const Group&);
  bool operator==(const Group& a) const { return globGroup_Index == a.globGroup_Index; }
  bool operator!=(const Group& a) const { return globGroup_Index != a.globGroup_Index; }
  friend ostream& operator<<(ostream &strm,Group& ob){ ob.show_info(); return strm;  }

  void set(object);
  void show_info();
  void save(boost::property_tree::ptree& pt,std::string path);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);


};

void save(boost::property_tree::ptree& pt,std::string path,vector<Group>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Group>& vt,int& status);



}// namespace libmol
}// namespace libchemobjects
}// liblibra


#endif // GROUP_H
