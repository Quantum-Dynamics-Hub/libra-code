#ifndef UNIVERSE_H
#define UNIVERSE_H

#include "../../mmath/libmmath.h"
using namespace libmmath;

#include "Element.h"

namespace libchemobjects{
namespace libuniverse{



class Universe{

  //----------- Databases -----------------
  // Periodic Table of Elements
  map<std::string,Element>     PeriodicTable;

  // GMP or other electronegativities and some atomic properties
  //map<std::string,_GMP>     GMP;

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const Universe&); // Copies the content which is defined

public:

  //----------- Basic class operations ---------------------------
  // Defined in Universe.cpp
  Universe();                   // constructor
  Universe(const Universe&);    // copy-constructor
 ~Universe();                   // destructor

  Universe& operator=(const Universe&); // assignment operator
  void show_info();
  void set(object);

  //-------------- Getters, setters -----------------------------
  // Defined in Universe_aux.cpp
  void Add_Element_To_Periodic_Table(Element rec);
  Element Get_Element(std::string);

};


}// namespace libuniverse
}// namespace libchemobjects



#endif // UNIVERSE_H

