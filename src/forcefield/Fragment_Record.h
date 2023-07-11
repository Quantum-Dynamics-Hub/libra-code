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

#ifndef FRAGMENT_RECORD_H
#define FRAGMENT_RECORD_H


//#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"


/// liblibra namespace
namespace liblibra{

//using namespace liblinalg;
using namespace libio;

namespace libforcefield{


class Fragment_Record; 

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
  Fragment_Record(const Fragment_Record& at){ *this = at; }
  Fragment_Record& operator=(const Fragment_Record&);
  void save(boost::property_tree::ptree& pt,std::string path);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);

   friend bool operator == (const Fragment_Record& e1, const Fragment_Record& e2) {    return &e1 == &e2;  }
   friend bool operator != (const Fragment_Record& e1, const Fragment_Record& e2) {    return !(e1==e2);   } 



};

void save(boost::property_tree::ptree& pt,std::string path,vector<Fragment_Record>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Fragment_Record>& vt,int& status);

typedef std::vector< Fragment_Record > Fragment_RecordList; ///< Type containing the vector of Electronic objects



}// namespace libforcefield
}// namespace liblibra



#endif
