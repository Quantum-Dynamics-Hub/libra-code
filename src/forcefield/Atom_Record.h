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

#ifndef ATOM_RECORD_H
#define ATOM_RECORD_H


//#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"


/// liblibra namespace
namespace liblibra{

//using namespace liblinalg;
using namespace libio;

namespace libforcefield{

class Atom_Record;

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
   int show_info() const;
   void merge(const Atom_Record&);
   Atom_Record(const Atom_Record& at){ *this = at; }

   Atom_Record& operator=(const Atom_Record&);
   void save(boost::property_tree::ptree& pt,std::string path);
   void load(boost::property_tree::ptree& pt,std::string path,int& status);

   friend bool operator == (const Atom_Record& e1, const Atom_Record& e2) {    return &e1 == &e2;  }
   friend bool operator != (const Atom_Record& e1, const Atom_Record& e2) {    return !(e1==e2);   } 


};

void save(boost::property_tree::ptree& pt,std::string path,vector<Atom_Record>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Atom_Record>& vt,int& status);

typedef std::vector< Atom_Record > Atom_RecordList; ///< Type containing the vector of Electronic objects


}// namespace libforcefield
}// namespace liblibra



#endif
