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

#ifndef ANGLE_RECORD_H
#define ANGLE_RECORD_H


//#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"


/// liblibra namespace
namespace liblibra{

//using namespace liblinalg;
using namespace libio;

namespace libforcefield{

class Angle_Record;

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
   int Angle_coordination;          int is_Angle_coordination;

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
   is_Angle_coordination = 0;

   }

   int set(object at);
   int show_info();
   void merge(const Angle_Record&);
   Angle_Record(const Angle_Record& at){ *this = at; }
   Angle_Record& operator=(const Angle_Record&);
   void save(boost::property_tree::ptree& pt,std::string path);
   void load(boost::property_tree::ptree& pt,std::string path,int& status);

   friend bool operator == (const Angle_Record& e1, const Angle_Record& e2) {    return &e1 == &e2;  }
   friend bool operator != (const Angle_Record& e1, const Angle_Record& e2) {    return !(e1==e2);   } 


};

void save(boost::property_tree::ptree& pt,std::string path,vector<Angle_Record>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Angle_Record>& vt,int& status);

typedef std::vector< Angle_Record > Angle_RecordList; ///< Type containing the vector of Electronic objects




}// namespace libforcefield
}// liblibra


#endif
