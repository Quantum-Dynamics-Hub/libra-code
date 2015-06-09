#include "Universe.h"

namespace libchemobjects{
namespace libuniverse{


void Universe::init_variables(){

}

void Universe::copy_content(const Universe& u){
//  if(u.PeriodicTable.size()>0){ PeriodicTable = u.PeriodicTable; }

}

Universe::Universe(){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
}

Universe::Universe(const Universe& u){
  /********************
    Copy constructor
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of bar object which is defined
  copy_content(u);
}

Universe& Universe::operator=(const Universe& u){
  /********************
    Assignment operator
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of bar object which is defined
  copy_content(u);
  return *this;
}


Universe::~Universe(){ ;; }


}// namespace libuniverse
}// namespace libchemobjects


