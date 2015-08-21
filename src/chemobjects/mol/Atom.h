#ifndef ATOM_H
#define ATOM_H

#include "../../dyn/rigidbody/librigidbody.h"
#include "../universe/libuniverse.h"
using namespace libdyn::librigidbody;

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libchemobjects{

using namespace libuniverse;

namespace libmol{


class Atom{
  
  Universe* universe;

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const Atom&); // Copies the content which is defined

public:
  //--------- Atom identification ----------------------------
  int         Atom_id;             int is_Atom_id;

  //--------- Topology = Data model -------------------------
  int globAtom_Index;
  int locAtom_Index;
  int globGroup_Index;
  int globMolecule_Index;
  vector<int> globAtom_Adjacent_Atoms;

  //--------- Dynamical properties of the atom --------------
  RigidBody   Atom_RB;             int is_Atom_RB;
  RigidBody   Atom_RB_old;         int is_Atom_RB_old; // previous position
  double      Atom_displ2;         int is_Atom_displ2; // square of atomic displacements 

  //----------Basic physical properties of the atom ---------
  int         Atom_Z;              int is_Atom_Z;
  std::string Atom_element;        int is_Atom_element;
  double      Atom_atomic_radius;  int is_Atom_atomic_radius;
  double      Atom_charge;         int is_Atom_charge;
  double      Atom_electronegativity;   int is_Atom_electronegativity;

  //--------- Atom-in-Molecule topological properties --------------
  double      Atom_formal_charge;  int is_Atom_formal_charge;
  int         Atom_coordination;   int is_Atom_coordination; // defines an equilibrium angle
  std::string Atom_functional_group; int is_Atom_functional_group;// Name of the functional group
                                                                  // to which atom belongs
  vector<int> Atom_ring_sizes;     int is_Atom_ring_sizes;   // Sizes of rings to which this atom belong
  int         Atom_min_ring_size;  int is_Atom_min_ring_size;

  //----------- Force-field related properties ---------------
  std::string Atom_ff_type;         int is_Atom_ff_type;
  double Atom_Zeff;                 int is_Atom_Zeff;

  //----------- Electronic structure-related properties ---------------
  double Atom_mull_charge_gross;         int is_Atom_mull_charge_gross;
  double Atom_mull_charge_net;           int is_Atom_mull_charge_net;
 
  //---------- Keep this for a while --------------------
//  int         Atom_ff_int_type;     int is_Atom_ff_int_type;
//  int         Atom_is_surface_atom; int is_Atom_is_surface_atom;
//  int         Atom_surface_index;   int is_Atom_surface_index;
//  int         Atom_is_basis_atom;   int is_Atom_is_basis_atom;
//  int         Atom_is_C60_CT;       int is_Atom_is_C60_CT;


  //--------- Methods -----------------
  Atom(Universe&);
  Atom(Universe&,boost::python::dict);
  Atom(const Atom&); // Copy constructor
 ~Atom();            // Destructor
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



#endif // ATOM_H
