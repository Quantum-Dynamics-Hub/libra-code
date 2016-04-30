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

#ifndef FORCE_FIELDS_H
#define FORCE_FIELDS_H


#include "../../../../mmath/libmmath.h"
using namespace libmmath;

#include "../../../../io/libio.h"
using namespace libio;



namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_mm{
namespace libforcefield{



class Atom_Record;
class Bond_Record;
class Angle_Record;
class Dihedral_Record;
class Fragment_Record;
class ForceField;


class Atom_Record{

public:

//-------------------- Type markers ----------------------
   std::string Atom_ff_type;       int is_Atom_ff_type;
   std::string Atom_ff_type_H;     int is_Atom_ff_type_H;
   int         Atom_ff_int_type;   int is_Atom_ff_int_type;
   std::string Atom_element;       int is_Atom_element;
   int         Atom_atomic_number; int is_Atom_atomic_number;
//-------------------- Data -----------------------------
   // General atom properties
   double      Atom_electronegativity; int is_Atom_electronegativity;
   double      Atom_partial_charge;    int is_Atom_partial_charge;
   // Atom type equivalences (used in many force fields)
   int         Atom_ff_eq_int_type2; int is_Atom_ff_eq_int_type2;
   int         Atom_ff_eq_int_type3; int is_Atom_ff_eq_int_type3;
   int         Atom_ff_eq_int_type4; int is_Atom_ff_eq_int_type4;
   int         Atom_ff_eq_int_type5; int is_Atom_ff_eq_int_type5;

   // UFF and DREIDING force field parameters
   double      Atom_radius;        int is_Atom_radius;
   double      Atom_Z_star;        int is_Atom_Z_star; // Effective charge of atom(UFF)
   double      Atom_theta;         int is_Atom_theta;  // Valence angle (in degrees)
   double      Atom_sigma;         int is_Atom_sigma;  // vdw radius of atom
   double      Atom_epsilon;       int is_Atom_epsilon;// vdw well depth for atom (D in terms of UFF)
   double      Atom_GMP;           int is_Atom_GMP;    // Generalized electronegativity

   // MMFF94 force field parameters
   int         Atom_crd;           int is_Atom_crd;
   int         Atom_val;           int is_Atom_val;
   int         Atom_pilp;          int is_Atom_pilp;
   int         Atom_mltb;          int is_Atom_mltb;
   int         Atom_arom;          int is_Atom_arom;
   int         Atom_lin;           int is_Atom_lin;
   int         Atom_sbmb;          int is_Atom_sbmb;
   double      Atom_alpha;         int is_Atom_alpha; 
   double      Atom_N_eff;         int is_Atom_N_eff;
   double      Atom_A_scale;       int is_Atom_A_scale;
   double      Atom_G_scale;       int is_Atom_G_scale;
   std::string Atom_DAN;           int is_Atom_DAN;
   double      Atom_pbci;          int is_Atom_pbci;
   double      Atom_fcadj;         int is_Atom_fcadj;


   // ESFF force field parameters
   double      Atom_dative;        int is_Atom_dative;
   double      Atom_brdr1;         int is_Atom_brdr1;
   double      Atom_brdr2;         int is_Atom_brdr2;
   double      Atom_brdr3;         int is_Atom_brdr3;

   // Sutton-Chen parameters 
   double      Atom_such_n;        int is_Atom_such_n;
   double      Atom_such_m;        int is_Atom_such_m;
   double      Atom_such_a;        int is_Atom_such_a;
   double      Atom_such_D;        int is_Atom_such_D;
   double      Atom_such_c;        int is_Atom_such_c;


   Atom_Record(){
 
   is_Atom_ff_type     = 0;
   is_Atom_ff_type_H   = 0;
   is_Atom_ff_int_type = 0;
   is_Atom_element     = 0;
   is_Atom_atomic_number = 0;

   is_Atom_electronegativity = 0;
   is_Atom_partial_charge = 0;

   is_Atom_ff_eq_int_type2 = 0;
   is_Atom_ff_eq_int_type3 = 0;
   is_Atom_ff_eq_int_type4 = 0;
   is_Atom_ff_eq_int_type5 = 0;
   
   is_Atom_radius      = 0;
   is_Atom_Z_star      = 0;
   is_Atom_theta       = 0;
   is_Atom_sigma       = 0;
   is_Atom_epsilon     = 0;
   is_Atom_GMP         = 0;

   is_Atom_crd         = 0; 
   is_Atom_val         = 0;
   is_Atom_pilp        = 0;
   is_Atom_mltb        = 0;
   is_Atom_arom        = 0;
   is_Atom_lin         = 0;
   is_Atom_sbmb        = 0;
   is_Atom_alpha       = 0;   
   is_Atom_N_eff       = 0;
   is_Atom_A_scale     = 0;
   is_Atom_G_scale     = 0;
   is_Atom_DAN         = 0; 
   is_Atom_pbci        = 0;
   is_Atom_fcadj       = 0;

   is_Atom_dative      = 0;
   is_Atom_brdr1       = 0;
   is_Atom_brdr2       = 0;
   is_Atom_brdr3       = 0;

   is_Atom_such_n      = 0;
   is_Atom_such_m      = 0;
   is_Atom_such_a      = 0;
   is_Atom_such_D      = 0;
   is_Atom_such_c      = 0;

   }

   int set(object at);
   int show_info();
   void merge(const Atom_Record&);
   Atom_Record& operator=(const Atom_Record&);
   void save(boost::property_tree::ptree& pt,std::string path);
   void load(boost::property_tree::ptree& pt,std::string path,int& status);


};

void save(boost::property_tree::ptree& pt,std::string path,vector<Atom_Record>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Atom_Record>& vt,int& status);


class Bond_Record{

public:

//-------------------- Type markers ----------------------
   std::string Atom1_ff_type;       int is_Atom1_ff_type;
   std::string Atom2_ff_type;       int is_Atom2_ff_type;
   int         Atom1_ff_int_type;   int is_Atom1_ff_int_type;
   int         Atom2_ff_int_type;   int is_Atom2_ff_int_type;
   int         Atom1_atomic_number; int is_Atom1_atomic_number;
   int         Atom2_atomic_number; int is_Atom2_atomic_number;
   std::string Atom1_element;       int is_Atom1_element;
   std::string Atom2_element;       int is_Atom2_element;

   

   int         Bond_type_index;     int is_Bond_type_index;
//-------------------- Data -----------------------------

   double      Bond_r_eq;           int is_Bond_r_eq;
   double      Bond_k_bond;         int is_Bond_k_bond;
   double      Bond_D_bond;         int is_Bond_D_bond;
   double      Bond_alpha;          int is_Bond_alpha;
   double      Bond_r_eq_ref;       int is_Bond_r_eq_ref;
   double      Bond_k_bond_ref;     int is_Bond_k_bond_ref;
   double      Bond_shift_elec;     int is_Bond_shift_elec;
   double      Bond_bci;            int is_Bond_bci;  // bond charge increment
 
   // For MALINA force field
   double      Bond_wij;            int is_Bond_wij;
   double      Bond_wij_1;          int is_Bond_wij_1;
   double      Bond_wij_2;          int is_Bond_wij_2;
   double      Bond_alpij;          int is_Bond_alpij;
   double      Bond_alpij_1;        int is_Bond_alpij_1;
   double      Bond_alpij_2;        int is_Bond_alpij_2;
   double      Bond_rij_1;          int is_Bond_rij_1;
   double      Bond_rij_2;          int is_Bond_rij_2;


   Bond_Record(){
  
   is_Atom1_ff_type = 0;
   is_Atom2_ff_type = 0;
   is_Atom1_ff_int_type = 0;
   is_Atom2_ff_int_type = 0;
   is_Atom1_atomic_number = 0;
   is_Atom2_atomic_number = 0;
   is_Atom1_element = 0;
   is_Atom2_element = 0;
   is_Bond_type_index = 0;


   is_Bond_r_eq = 0;
   is_Bond_k_bond = 0;
   is_Bond_D_bond = 0;
   is_Bond_alpha = 0;
   is_Bond_r_eq_ref = 0;
   is_Bond_k_bond_ref = 0;
   is_Bond_shift_elec = 0;
   is_Bond_bci = 0;
   is_Bond_wij = 0;
   is_Bond_wij_1 = 0;
   is_Bond_wij_2 = 0;
   is_Bond_alpij = 0;
   is_Bond_alpij_1 = 0;
   is_Bond_alpij_2 = 0;
   is_Bond_rij_1 = 0;
   is_Bond_rij_2 = 0;

    
   }

   int set(object at);
   int show_info();
   void merge(const Bond_Record&); // similar to = operator but does keeps old properties 
                                   // if they are defined
   Bond_Record& operator=(const Bond_Record&);
   void save(boost::property_tree::ptree& pt,std::string path);
   void load(boost::property_tree::ptree& pt,std::string path,int& status);



};

void save(boost::property_tree::ptree& pt,std::string path,vector<Bond_Record>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Bond_Record>& vt,int& status);


class Angle_Record{

public:

//-------------------- Type markers ----------------------
   std::string Atom1_ff_type;       int is_Atom1_ff_type;
   std::string Atom2_ff_type;       int is_Atom2_ff_type;
   std::string Atom3_ff_type;       int is_Atom3_ff_type;
   
   int         Atom1_ff_int_type;   int is_Atom1_ff_int_type;
   int         Atom2_ff_int_type;   int is_Atom2_ff_int_type;
   int         Atom3_ff_int_type;   int is_Atom3_ff_int_type;

   int         Angle_type_index;    int is_Angle_type_index;
//-------------------- Data -----------------------------

   double Angle_theta_eq;           int is_Angle_theta_eq;
   double Angle_k_angle;            int is_Angle_k_angle;
   double Angle_r_eq;               int is_Angle_r_eq;
   double Angle_k_ub;               int is_Angle_k_ub;
   double Angle_kijk_sb;            int is_Angle_kijk_sb; // _sb = stretch-bend
   double Angle_kkji_sb;            int is_Angle_kkji_sb;

   Angle_Record(){

   is_Atom1_ff_type = 0;
   is_Atom2_ff_type = 0;
   is_Atom3_ff_type = 0;


   is_Atom1_ff_int_type = 0;
   is_Atom2_ff_int_type = 0;
   is_Atom3_ff_int_type = 0;

   is_Angle_type_index  = 0;

   is_Angle_theta_eq = 0;
   is_Angle_k_angle  = 0;
   is_Angle_r_eq     = 0;
   is_Angle_k_ub     = 0;
   is_Angle_kijk_sb  = 0;
   is_Angle_kkji_sb  = 0;

   }

   int set(object at);
   int show_info();
   void merge(const Angle_Record&);
   Angle_Record& operator=(const Angle_Record&);
   void save(boost::property_tree::ptree& pt,std::string path);
   void load(boost::property_tree::ptree& pt,std::string path,int& status);


};

void save(boost::property_tree::ptree& pt,std::string path,vector<Angle_Record>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Angle_Record>& vt,int& status);


class Dihedral_Record{

public:

// This data class is also used for impropers (oops)
// In that case the atom2 will be treated as central atom!

//-------------------- Type markers ----------------------
   std::string Atom1_ff_type;       int is_Atom1_ff_type;
   std::string Atom2_ff_type;       int is_Atom2_ff_type;
   std::string Atom3_ff_type;       int is_Atom3_ff_type;
   std::string Atom4_ff_type;       int is_Atom4_ff_type;

   int         Atom1_ff_int_type;   int is_Atom1_ff_int_type;
   int         Atom2_ff_int_type;   int is_Atom2_ff_int_type;
   int         Atom3_ff_int_type;   int is_Atom3_ff_int_type;
   int         Atom4_ff_int_type;   int is_Atom4_ff_int_type;

   int         Dihedral_type_index; int is_Dihedral_type_index;
//-------------------- Data -----------------------------

   double Dihedral_vphi;            int is_Dihedral_vphi;
   double Dihedral_vphi1;           int is_Dihedral_vphi1;
   double Dihedral_vphi2;           int is_Dihedral_vphi2;
   double Dihedral_vphi3;           int is_Dihedral_vphi3;
   double Dihedral_phase;           int is_Dihedral_phase;
   int    Dihedral_mult;            int is_Dihedral_mult;

   Dihedral_Record(){

   is_Atom1_ff_type = 0;
   is_Atom2_ff_type = 0;
   is_Atom3_ff_type = 0;
   is_Atom4_ff_type = 0;

   is_Atom1_ff_int_type = 0;
   is_Atom2_ff_int_type = 0;
   is_Atom3_ff_int_type = 0;
   is_Atom4_ff_int_type = 0;

   is_Dihedral_type_index = 0;

   is_Dihedral_vphi  = 0;
   is_Dihedral_vphi1 = 0;
   is_Dihedral_vphi2 = 0;
   is_Dihedral_vphi3 = 0;
   is_Dihedral_phase = 0;
   is_Dihedral_mult  = 0;

   }

   int set(object at);
   int show_info();
   void merge(const Dihedral_Record&);
   Dihedral_Record& operator=(const Dihedral_Record&);
   void save(boost::property_tree::ptree& pt,std::string path);
   void load(boost::property_tree::ptree& pt,std::string path,int& status);


};

void save(boost::property_tree::ptree& pt,std::string path,vector<Dihedral_Record>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Dihedral_Record>& vt,int& status);


class Fragment_Record{

public:

//-------------------- Type markers ----------------------
  std::string Fragment_ff_type;       int is_Fragment_ff_type;
  int         Fragment_ff_int_type;   int is_Fragment_ff_int_type;
//-------------------- Data -----------------------------
  // General fragment properties
  // Gay-Berne properties:
  double Fragment_di;                 int is_Fragment_di;  // fragment depth 
  double Fragment_li;                 int is_Fragment_li;  // fragment length
  double Fragment_e0;                 int is_Fragment_e0;  // epsilon - interaction strength
  double Fragment_rat;                int is_Fragment_rat; // ratio of eps(E) to eps(S)
                                                           // eps(E)- well depth of end-to-end/face-to-face
                                                           // eps(S)- well depth of side-by-side
  double Fragment_dw;                 int is_Fragment_dw;  // correction parameter
  double Fragment_mu;                 int is_Fragment_mu;
  double Fragment_nu;                 int is_Fragment_nu;

  // Constructor
  Fragment_Record(){
    is_Fragment_ff_type = 0;
    is_Fragment_ff_int_type = 0;
    is_Fragment_di = 0;
    is_Fragment_li = 0;
    is_Fragment_e0 = 0;
    is_Fragment_rat = 0;
    is_Fragment_dw = 0;
    is_Fragment_mu = 0;
    is_Fragment_nu = 0;
  }
  int set(object at);
  int show_info();
  void merge(const Fragment_Record&);
  Fragment_Record& operator=(const Fragment_Record&);
  void save(boost::property_tree::ptree& pt,std::string path);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);


};

void save(boost::property_tree::ptree& pt,std::string path,vector<Fragment_Record>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Fragment_Record>& vt,int& status);


class ForceField{

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const ForceField&); // Copies the content which is defined
  void extract_dictionary(boost::python::dict);

public:

   std::string ForceField_Name;       int is_ForceField_Name;

   std::string bond_functional;       int is_bond_functional;
   std::string angle_functional;      int is_angle_functional;
   std::string dihedral_functional;   int is_dihedral_functional;
   std::string oop_functional;        int is_oop_functional;
   std::string vdw_functional;        int is_vdw_functional;
   std::string elec_functional;       int is_elec_functional;
   std::string mb_functional;         int is_mb_functional;
   std::string cg_functional;         int is_cg_functional;
   std::string mb_excl_functional;    int is_mb_excl_functional;

//   std::string stress_opt;            int is_sress_opt;

   std::string system_pbc;            int is_system_pbc;
   int  pbc_degree_x;                 int is_pbc_degree_x;
   int  pbc_degree_y;                 int is_pbc_degree_y;
   int  pbc_degree_z;                 int is_pbc_degree_z;
   int  reciprocal_degree_x;          int is_reciprocal_degree_x;
   int  reciprocal_degree_y;          int is_reciprocal_degree_y;
   int  reciprocal_degree_z;          int is_reciprocal_degree_z;
   double  R_vdw_on;                  int is_R_vdw_on;
   double  R_vdw_off;                 int is_R_vdw_off;
   double  R_elec_on;                 int is_R_elec_on;
   double  R_elec_off;                int is_R_elec_off;
   double  R_vlist;                   int is_R_vlist;    // Verlet list radius
   double  elec_etha;                 int is_elec_etha;

   std::string sigma_comb_rule;       int is_sigma_comb_rule;
   std::string epsilon_comb_rule;     int is_epsilon_comb_rule;
   double vdw_scale12;                int is_vdw_scale12;
   double vdw_scale13;                int is_vdw_scale13;
   double vdw_scale14;                int is_vdw_scale14;
   double elec_scale12;               int is_elec_scale12;
   double elec_scale13;               int is_elec_scale13;
   double elec_scale14;               int is_elec_scale14;

   std::vector<Atom_Record>     Atom_Records;
   std::vector<Bond_Record>     Bond_Records;
   std::vector<Angle_Record>    Angle_Records;
   std::vector<Dihedral_Record> Dihedral_Records;
   std::vector<Dihedral_Record> Improper_Records;
   std::vector<Fragment_Record> Fragment_Records;

  //----------- Basic class operations ---------------------------
  // Defined in ForceField.cpp
  ForceField();                   // constructor
  ForceField(boost::python::dict);
  ForceField(const ForceField&);  // copy-constructor
 ~ForceField();                   // destructor

  ForceField& operator=(const ForceField&); // assignment operator
  void show_info();
//  void set(object)
  void set(boost::python::dict);
  void save(boost::property_tree::ptree& pt,std::string path);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);



   int Atom_Record_Index(int);
   int Atom_Record_Index(std::string);
   int Atom_Record_Index_by_Element(std::string);
   int Atom_Record_Index_by_Element(std::string,vector<int>&);
   int Atom_Record_Index_by_Element(int,vector<int>&);

   int Bond_Record_Index(int,int);
   int Bond_Record_Index(int,int,int);
   int Bond_Record_Index(int,int,int,vector<int>&);
   int Bond_Record_Index(std::string,std::string);
   int Bond_Record_Index(std::string,std::string,int,vector<int>&);
   int Bond_Record_Index_by_Element(int,int,int,vector<int>&);
   int Bond_Record_Index_by_Element(std::string,std::string,int,vector<int>&);
   
   int Angle_Record_Index(int,int,int);
   int Angle_Record_Index(int,int,int,int,int,vector<int>&);
   int Angle_Record_Index(std::string,std::string,std::string);
   int Angle_Record_Index(std::string,std::string,std::string,int,int,vector<int>&);

   int Dihedral_Record_Index(int,int,int,int);
   int Dihedral_Record_Index(int,int,int,int,int,int,vector<int>&);
   int Dihedral_Record_Index(std::string,std::string,std::string,std::string);
   int Dihedral_Record_Index(std::string,std::string,std::string,std::string,int,int,vector<int>&);
   int Dihedral_Record_Index(int,int);

   int Improper_Record_Index(int,int,int,int);
//   int Improper_Record_Index(int,int,int,int,int,vector<int>&);
//   int Improper_Record_Index(std::string,std::string,std::string,std::string,int,int,vector<int>&);
//   int Improper_Record_Index(int,int);

   int Fragment_Record_Index(int);
   int Fragment_Record_Index(std::string);




   int Add_Atom_Record(Atom_Record rec);
   int Add_Bond_Record(Bond_Record rec);
   int Add_Angle_Record(Angle_Record rec);
   int Add_Angle_Record(Angle_Record rec,int);
   int Add_Dihedral_Record(Dihedral_Record rec);
   int Add_Dihedral_Record(Dihedral_Record rec,int);
   int Add_Improper_Record(Dihedral_Record rec);
   int Add_Fragment_Record(Fragment_Record rec);
//   int Add_Improper_Record(Dihedral_Record rec,int);


   //-------- Methods -----------------
   int show_atom_records();
   int show_bond_records();
   int show_angle_records();
   int show_dihedral_records();
   int show_improper_records();
   int show_fragment_records();


   // Make them private
   string tip3p_type(string elt,int geometry,string func_grp,int min_ring,int& coordination); 
   string uff_type(string elt,int geometry,string func_grp,int min_ring,int& coordination);
   string dreiding_type(string elt,int geometry,string func_grp,int min_ring,int& coordination);
   string tripos_type(string elt,int geometry,string func_grp,int min_ring,int& coordination);
   string gaff_type(string elt,int geometry,string func_grp,int min_ring,int& coordination);
   string mmff94_type(string elt,int geometry,string func_grp,int min_ring,int& coordination);

   // Defined in ForceField_methods.cpp
   void bond_r0_rule(string,string,double,double&,int&);
   void bond_K_rule(string,string,double,double&,int&);
   void bond_D_rule(string,string,double,double&,int&);
   void bond_alpha_rule(string,string,double,double&,int&);

   // Defined in ForceField_methods1.cpp
   void vdw_sigma_rule(string ff_type1,string ff_type2,double& sigma,int& is_sigma);
   void vdw_epsilon_rule(string ff_type1,string ff_type2,double& epsilon,int& is_epsilon);

   // Defined in ForceField_methods2.cpp
   void set_functionals(boost::python::dict);

   // Defined in ForceField_methods3.cpp
   void angle_theta_0_rule(std::string ff_type1,std::string ff_type2,std::string ff_type3,
                           double bond_order12,double bond_order23,int coordination,
                           double& theta_0,int& is_theta_0);
   void angle_k_theta_rule(std::string ff_type1,std::string ff_type2,std::string ff_type3,
                           double bond_order12,double bond_order23,int coordination,
                           double& k_theta,int& is_k_theta);

   // Defined in ForceField_methods4.cpp
   void dihedral_rule(std::string ff_type1,std::string ff_type2,
                      std::string ff_type3,std::string ff_type4,
                      double bond_order12,double bond_order23,double bond_order34,
                      double& Vphi,int& is_Vphi,
                      double& phi0,int& is_phi0,
                      int& n,int& is_n);

   // Defined in ForceField_methods5.cpp
   void oop_rule(std::string ff_type2, double& K, int& is_K, 
                 double& C0, int& is_C0,  double& C1, int& is_C1,
                 double& C2, int& is_C2);



   // These are public:
   string get_atom_type(string elt,int geometry,string func_grp,int min_ring,int& coordination){
     string res = "";
     if(ForceField_Name == "UFF"){ res=uff_type(elt,geometry,func_grp,min_ring,coordination);   }
     else if(ForceField_Name == "DREIDING"){ res=dreiding_type(elt,geometry,func_grp,min_ring,coordination);   }
     else if(ForceField_Name == "GAFF"){ res=gaff_type(elt,geometry,func_grp,min_ring,coordination);   }
     else if(ForceField_Name == "MMFF94"){ res=mmff94_type(elt,geometry,func_grp,min_ring,coordination);   }
     else if(ForceField_Name == "TRIPOS"){ res=tripos_type(elt,geometry,func_grp,min_ring,coordination);   }
     else if(ForceField_Name == "TIP3P"){ res=tip3p_type(elt,geometry,func_grp,min_ring,coordination);   }
     return res;
   }
   // ForceField_methods.cpp
   int is_valid_atom_type(std::string);
   int is_valid_fragment_type(std::string);
   int get_bond_parameters(string, string, double, map<string,double>&);// ForceField_methods.cpp
   int get_angle_parameters(string,string,string,double,double,int,map<string,double>&);  // ForceField_methods3.cpp
   int get_dihedral_parameters(string,string,string,string,double,double,double,map<string,double>&); // ForceField_methods4.cpp
   int get_oop_parameters(string,string,string,string,map<string,double>&); // ForceField_methods5.cpp
   // ForceField_methods1.cpp
   int get_vdw_parameters(string,string,string,map<string,double>&); // ForceField_methods1.cpp
   int get_elec_parameters(string,string,string,double,double,int,int,map<string,double>&); // ForceField_methods6.cpp
   int get_vdw_parameters(int sz,vector<string> types,double** epsilon, double** sigma); // ForceField_methods1.cpp
   // ForceField_methods7.cpp
   int set_ff_charges(int, vector<string>, VECTOR**, double**);
   // ForceField_methods8.cpp
   int get_mb_parameters(map<string,double>& prms);
   int get_mb_excl_parameters(map<string,double>& prms);
   // ForceField_methods9.cpp
   int set_ff_epsilon_and_sigma(int, vector<string> ,double**, double**);
   // ForceField_method10.cpp
   int get_cg_parameters(map<string,double>& prms);

};


}// namespace libforcefield
}// namespace libhamiltonian_mm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian



#endif
