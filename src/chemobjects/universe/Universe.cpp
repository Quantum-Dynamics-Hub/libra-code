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
  \file Universe.cpp
  \brief The file implements the Universe class member functions    
*/

#include "Universe.h"

/// liblibra namespace
namespace liblibra{

/// libchemobjects namespace
namespace libchemobjects{

/// libuniverse namespace
namespace libuniverse{


void Universe::init_variables(){

}

void Universe::copy_content(const Universe& u){
//  if(u.PeriodicTable.size()>0){ PeriodicTable = u.PeriodicTable; }

}

Universe::Universe(){
/**
  \brief Constructor
  Initialize variables to default values
*/

  // Initialize variables to default values
  init_variables();
}

Universe::Universe(const Universe& u){
/**
  \brief Copy constructor
  Initialize variables to default values
  Copy content of the input object which is defined
*/
  // Initialize variables to default values
  init_variables();
  // Copy content of the input object which is defined
  copy_content(u);
}

Universe& Universe::operator=(const Universe& u){
/**
  \brief Assignment operator
  Initialize variables to default values
  Copy content of the input object which is defined
*/
  // Initialize variables to default values
  init_variables();
  // Copy content of the input object which is defined
  copy_content(u);
  return *this;
}


Universe::~Universe(){ ;; }


}// namespace libuniverse
}// namespace libchemobjects

}// liblibra
