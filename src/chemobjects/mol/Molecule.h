#ifndef MOLECULE_H
#define MOLECULE_H

#include "../../dyn/rigidbody/librigidbody.h"
using namespace libdyn::librigidbody;

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libchemobjects{
namespace libmol{


class Molecule{

//  Group Group_data; // Group data for molecule as a group

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const Molecule&); // Copies the content which is defined

public:

  //---------- Topology = Data Model -------------------
  int globMolecule_Size;
  int locMolecule_Size;
  int Molecule_Size;
// Since one group(fragment) may be composed
// from different molecules or pieces of different
// molecules the bijection is violated. Thus we
// now adopting the idea of atoms in a molecule,
// rather than groups in a molecule!
//  vector<int> globGroup_Index;
//  vector<int> locGroup_Index;
  vector<int> globAtom_Index;
  vector<int> locAtom_Index;
  
  int globMolecule_Index;
  int locMolecule_Index;

  int Molecule_Number_of_bonds;
  int Molecule_Number_of_angles;
  int Molecule_Number_of_dihedrals;
  int Molecule_Number_of_impropers;

  //--------- General properties ----------------------------
  std::string Molecule_name;         int is_Molecule_name;
  int         Molecule_id;           int is_Molecule_id;

  //--------- Dynamic variables -----------------------------
  RigidBody   Molecule_RB;           int is_Molecule_RB;


//----------- Basic class operations ---------------------------
// Defined in Molecule.cpp
  Molecule();                // constructor
  Molecule(const Molecule&); // copy-constructor
 ~Molecule();                // destructor

  Molecule& operator=(const Molecule&); // assignment operator
  bool operator==(const Molecule& a) const { return globMolecule_Index == a.globMolecule_Index; }
  bool operator!=(const Molecule& a) const { return globMolecule_Index != a.globMolecule_Index; }

  void show_info();
  void set(object);

};

}// namespace libmol
}// namespace libchemobjects


#endif // MOLECULE_H
