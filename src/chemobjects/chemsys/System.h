#ifndef SYSTEM_H
#define SYSTEM_H

#include "../../mmath/libmmath.h"
using namespace libmmath;
using namespace libmmath::libgraph;

#include "../mol/libmol.h"


//#include "Interaction.h"
//#include "ForceField.h"
//#include "RigidBody.h"
//#include "Universe.h"


namespace libchemobjects{
using namespace libmol;

namespace libchemsys{


struct connect{
public:
  int self;
  vector<int> others;
};


class System{

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const System&); // Copies the content which is defined


//---------- Internal variables --------------------
  int max_atom_id;     int is_max_atom_id;
  int max_fragment_id; int is_max_fragment_id;
  int max_molecule_id; int is_max_molecule_id;

  //----------- Defined in System_methods4.cpp ------------------
  // Chemistry related functions:
  int is(std::string,int,int,int,Atom**,vector<Atom*>&);
  void store(vector<Atom*>&,vector<Atom*>&);
  void set(vector<Atom*>&,std::string);

  //----------- Defined in System_methods2.cpp -----------------
  int is_in_vector(int indx,vector<int>& vect);


  //---------- Defined in System_aux.cpp ----------------------
  // Topology and builder related functions:
  void create_bond(int,int,int);
  void create_angle(int,int,int,int);
  void create_dihedral(int,int,int,int,int);
  void create_improper(int,int,int,int);


  //---------- Interactions parameters ------------------------
//  vector<Interaction> interactions;
//  vector<int> active_interactions;

  std::string stress_opt;int is_stress_opt;
  MATRIX3x3 stress_at;   int is_stress_at;
  MATRIX3x3 stress_fr;   int is_stress_fr;
  MATRIX3x3 stress_ml;   int is_stress_ml;
  MATRIX3x3 hessian;     int is_hessian;

  // RESPA auxiliary variables
  vector<VECTOR> respa_f_fast,respa_f_medium;
  vector<VECTOR> respa_t_fast,respa_t_medium;
  MATRIX3x3 respa_s_fast,respa_s_medium;
  double respa_E_fast,respa_E_medium;


public:

//---------- All objects comprising the system -------------------

  MATRIX*  GroupConnMatrix;     int is_GroupConnMatrix;
  GRAPH<Atom*,Group*> AtomGraph;int is_AtomGraph;

  MATRIX3x3 Box;              int is_Box;
  MATRIX3x3 Boxold;           int is_Boxold;
  VECTOR Box_origin;          int is_Box_origin;
  double dT_2;                int is_dT_2;  // cell vectors displacements 

  vector<Atom>  Atoms;        int Number_of_atoms;
  vector<Group> Bonds;        int Number_of_bonds;
  vector<Group> Angles;       int Number_of_angles;
  vector<Group> Dihedrals;    int Number_of_dihedrals;
  vector<Group> Impropers;    int Number_of_impropers;
  vector<Group> Pairs;        int Number_of_pairs;     // in future we need to get rid of this!
  vector<Group> Fragments;    int Number_of_fragments;
  vector<Group> Rings;        int Number_of_rings; // only smallest rings!
  vector<Molecule> Molecules; int Number_of_molecules;
//  vector<Surface> Surfaces;   int Number_of_surfaces;

  vector<int> Frag_bonds;    int Number_of_frag_bonds;
  vector<int> Frag_angles;   int Number_of_frag_angles;
  vector<int> Frag_dihedrals;int Number_of_frag_dihedrals;
  vector<int> Frag_impropers;int Number_of_frag_impropers;
  vector<int> Frag_pairs;    int Number_of_frag_pairs;
  vector<int> Surface_atoms; int Number_of_surface_atoms;

//---------- Properties ---------------
  std::string name;         int is_name;
  int         id;           int is_id;
  double      mass;         int is_mass;
  int         Nf_t;         int is_Nf_t;
  int         Nf_r;         int is_Nf_r;

  //----------- Basic class operations ---------------------------
  // Defined in System.cpp
  System();                // constructor
  System(const System&);   // copy-constructor
 ~System();                // destructor

  System& operator=(const System&); // assignment operator
  void show_info();
  void set(object);

  //--------- Defined in System_aux1.cpp ----------------------
  void move_atom_by_index(VECTOR&,int);
  void move_fragment_by_index(VECTOR&,int);
  void move_molecule_by_index(VECTOR&,int);

  //---------- Defined in System_methods.cpp  ------------------
  // Search and show functions
  void show_atoms();
  void show_bonds();
  void show_angles();
  void show_dihedrals();
  void show_impropers();
  void show_pairs();
  void show_frag_bonds();
  void show_frag_angles();
  void show_frag_dihedrals();
  void show_frag_impropers();
  void show_frag_pairs();
  void show_fragments();
  void show_rings();
  void show_molecules();

  void show_interactions();
  void show_interactions(std::string);


  int get_atom_index_by_atom_id(int);
  int get_fragment_index_by_fragment_id(int);
  int get_molecule_index_by_molecule_id(int);
  int Find_Bond(int,int);
  int Find_Frag_Pair(int,int);
  int Find_Angle(int,int);
  int Find_Angle(int,int,int);
  int Find_Dihedral(int,int,int,int);
  int Find_Improper(int);
  int is_12pair(int,int);
  int is_13pair(int,int);
  int is_14pair(int,int);
  int is_group_pair(int,int);

  //----------- Defined in System_methods1.cpp -----------
  // Topological functions
  void Generate_Connectivity_Matrix();
  void Assign_Rings();
  void DIVIDE_GRAPH(int,int, vector<int>&);  

  //----------- Defined in System_methods2.cpp ------------------
  // Builder functions
  void update_max_id();
  void CREATE_ATOM(Atom);
//  void CREATE_ATOM();
  void LINK_ATOMS(Atom&,Atom&);
  void LINK_ATOMS(int,int);
  void UPDATE_FRAG_TOPOLOGY();
  void ADD_ATOM_TO_FRAGMENT(int,int);
  void GROUP_ATOMS(boost::python::list,int);
  void CREATE_BONDS(boost::python::list,boost::python::dict);
  void CLONE_MOLECULE(int);

  //----------- Defined in System_methods3.cpp ------------------
  // Manipulation functions: Translations, rotations and updates
  void update_atoms_for_fragment(int);
  void update_fragments_for_molecule(int);
  void update_atoms_for_molecule(int);
  void rotate_atoms_of_fragment(int,MATRIX3x3&);
  void rotate_fragments_of_molecule(int,MATRIX3x3&);
  void rotate_atoms_of_molecule(int,MATRIX3x3&);
  void TRANSLATE_ATOM(double,VECTOR,int);
  void TRANSLATE_FRAGMENT(double,VECTOR,int);
  void TRANSLATE_MOLECULE(double,VECTOR,int);
  void ROTATE_FRAGMENT(double, VECTOR,int);
  void ROTATE_MOLECULE(double, VECTOR,int);

  //----------- Defined in System_methods4.cpp ------------------
  // Chemistry related functions:
  void determine_functional_groups(int assign_ring);


  //----------- Defined in System_methods5.cpp (extractors/converters) -------
  void extract_atomic_q(vector<double>& q);
  void set_atomic_q(vector<double>& q);

  void extract_atomic_p(vector<double>& p);
  void set_atomic_p(vector<double>& p);

  void extract_atomic_v(vector<double>& v);
  void set_atomic_v(vector<double>& v);

  void extract_atomic_f(vector<double>& f);
  void set_atomic_f(vector<double>& f);

  void extract_atomic_mass(vector<double>& mass);
  void set_atomic_mass(vector<double>& mass);


  void extract_fragment_q(vector<double>& q);
  void set_fragment_q(vector<double>& q);

  void extract_fragment_p(vector<double>& p);
  void set_fragment_p(vector<double>& p);

  void extract_fragment_v(vector<double>& v);
  void set_fragment_v(vector<double>& v);

  void extract_fragment_f(vector<double>& f);
  void set_fragment_f(vector<double>& f);

  void extract_fragment_mass(vector<double>& mass);
  void set_fragment_mass(vector<double>& mass);
  

  //------------- Defined in System_methods6.cpp ------------------
  void zero_atom_forces();
  void zero_fragment_forces();
  void zero_fragment_torques();
  void zero_forces();
  void zero_forces_and_torques();
  void update_fragment_forces();
  void update_fragment_torques();
  void update_fragment_forces_and_torques();

  void save_forces(vector<VECTOR>& frcs);
  void save_torques(vector<VECTOR>& trcs);
  void load_forces(vector<VECTOR>& frcs);
  void load_torques(vector<VECTOR>& trcs);
  void save_stress(MATRIX3x3& strs);  
  void increment_stress(MATRIX3x3& strs);
  void save_respa_state(std::string);
  void load_respa_state(std::string);


  void init_fragments();
  void init_molecules();

  void init_box_origin();
  void init_box();
  void init_box(double,double,double);
  void init_box(VECTOR,VECTOR,VECTOR);
  void apply_frag_pbc(std::string);

  void fix_fragment_translation(int);
  void fix_fragment_rotation(int);
  void fix_fragment(int);

//  double energy(std::string);
//  double energy();
//  double energy_respa(std::string respa_type);
  double ekin_tr();
  double ekin_tr_int(); // internal translational kinetic energy
  double ekin_rot();
  double volume();
  MATRIX3x3 pressure_tensor();


  //---------------- Defined in System_methods7.cpp -----------------
  void print_ent(std::string);
  void print_ent(std::string,int,std::string);
  void print_ent(std::string,boost::python::list);
  void print_ent(std::string,boost::python::list,int,std::string);

  void print_xyz(std::string,int);
  void print_xyz(std::string,int,std::string,int);


  
};




}// namespace libchemsys
}// namespace libchemobjects


#endif // SYSTEM_H
