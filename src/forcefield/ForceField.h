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

#ifndef FORCE_FIELDS_H
#define FORCE_FIELDS_H



#include "Atom_Record.h"
#include "Bond_Record.h"
#include "Angle_Record.h"
#include "Dihedral_Record.h"
#include "Fragment_Record.h"


#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"
#include "../Units.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;


namespace libforcefield{



//class Atom_Record;
//class Bond_Record;
//class Angle_Record;
//class Dihedral_Record;
//class Fragment_Record;
//class ForceField;



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

  void show_info() const;
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
   void oop_rule(std::string ff_type1, std::string ff_type2, std::string ff_type3, std::string ff_type4,
                 double& K, int& is_K, double& C0, int& is_C0,  double& C1, int& is_C1,
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
}// namespace liblibra



#endif
