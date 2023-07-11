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
/**
  \file Element.h
  \brief The file describes the Element class and related functions    
*/

#ifndef ELEMENT_H
#define ELEMENT_H
#include "../../io/libio.h"

/// liblibra namespace
namespace liblibra{


/// libchemobjects namespace
namespace libchemobjects{

/// libuniverse namespace
namespace libuniverse{


class Element{
/**
  \brief The class for keeping fundamental (and some custom) properties of chemical elements in (some) Universe
*/

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const Element&); // Copies the content which is defined

public:
  std::string Elt_name;           int is_Elt_name;   ///< actually it is an element symbol: H, He, Li, etc.
  int         Elt_id;             int is_Elt_id;     ///< Usually is the same as Elt_number
  double      Elt_mass;           int is_Elt_mass;   ///< Atomic mass of the element
  int         Elt_number;         int is_Elt_number; ///< The number of the element in the periodic table of elements
  int         Elt_nucleus_charge; int is_Elt_nucleus_charge; ///< The nucleus charge of given element
  int         Elt_period;         int is_Elt_period; ///< period in periodic table of elements
  int         Elt_group;          int is_Elt_group;  ///< group in periodic table of elements
  std::string Elt_block;          int is_Elt_block;  ///< element block; S, P, D or F
  int         Elt_red;            int is_Elt_red;    ///< value of red - for visualization/coloring applications
  int         Elt_green;          int is_Elt_green;  ///< value of green - for visualization/coloring applications
  int         Elt_blue;           int is_Elt_blue;   ///< value of blue - for visualization/coloring applications
  double      Elt_bond_length;    int is_Elt_bond_length; ///< closest distance in the element
  double      Elt_radius;         int is_Elt_radius; ///< some measure of the atomic radius (based on electron density)

  //----------- Basic class operations ---------------------------
  // Defined in Element.cpp
  Element();                ///< constructor
  Element(const Element&);  ///< copy-constructor
 ~Element();                ///< destructor
  Element& operator=(const Element&); ///< assignment operator

  void show_info();
  void set(object);
  void save(boost::property_tree::ptree& pt,std::string path);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);

};

void save(boost::property_tree::ptree& pt,std::string path,vector<Element>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Element>& vt,int& status);


}// namespace libuniverse
}// namespace libchemobjects

}// liblibra

#endif // ELEMENT_H

