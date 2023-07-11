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

#ifndef BOND_ORDER_H
#define BOND_ORDER_H


//#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"


/// liblibra namespace
namespace liblibra{

//using namespace liblinalg;
using namespace libio;

namespace libforcefield{



class Bond_Record;

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
   Bond_Record(const Bond_Record& at){ *this = at; }
   Bond_Record& operator=(const Bond_Record&);
   void save(boost::property_tree::ptree& pt,std::string path);
   void load(boost::property_tree::ptree& pt,std::string path,int& status);

   friend bool operator == (const Bond_Record& e1, const Bond_Record& e2) {    return &e1 == &e2;  }
   friend bool operator != (const Bond_Record& e1, const Bond_Record& e2) {    return !(e1==e2);   } 



};

void save(boost::property_tree::ptree& pt,std::string path,vector<Bond_Record>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Bond_Record>& vt,int& status);

typedef std::vector< Bond_Record > Bond_RecordList; ///< Type containing the vector of Electronic objects



}// namespace libforcefield
}// namespace liblibra



#endif
