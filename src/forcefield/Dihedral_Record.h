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

#ifndef DIHEDRAL_RECORD_H
#define DIHEDRAL_RECORD_H


//#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"


/// liblibra namespace
namespace liblibra{

//using namespace liblinalg;
using namespace libio;

namespace libforcefield{

class Dihedral_Record;

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
   Dihedral_Record(const Dihedral_Record& at){ *this = at; }
   Dihedral_Record& operator=(const Dihedral_Record&);
   void save(boost::property_tree::ptree& pt,std::string path);
   void load(boost::property_tree::ptree& pt,std::string path,int& status);

   friend bool operator == (const Dihedral_Record& e1, const Dihedral_Record& e2) {    return &e1 == &e2;  }
   friend bool operator != (const Dihedral_Record& e1, const Dihedral_Record& e2) {    return !(e1==e2);   } 



};

void save(boost::property_tree::ptree& pt,std::string path,vector<Dihedral_Record>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Dihedral_Record>& vt,int& status);

typedef std::vector< Dihedral_Record > Dihedral_RecordList; ///< Type containing the vector of Electronic objects



}// namespace libforcefield
}// namespace liblibra



#endif
