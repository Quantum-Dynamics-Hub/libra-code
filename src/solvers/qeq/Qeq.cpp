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

#include "Qeq.h"

void Qeq::set(object at){

 set_value(is_idempotential,    idempotential,    at,"idempotential");
 set_value(is_solution_method,  solution_method,  at,"solution_method");
 set_value(is_integral_method,  integral_method,  at,"integral_method");
 set_value(is_epsilon,          epsilon,          at,"epsilon");
 set_value(is_threshold,        threshold,        at,"threshold");
 set_value(is_MaxCount,         MaxCount,         at,"MaxCount");

}

void Qeq::show_info(){

  std::cout<<"CQeq properties:"<<std::endl;
  if(is_idempotential)   {std::cout<<"idempotential = "<<idempotential<<std::endl;   }
  if(is_solution_method) {std::cout<<"solution_method = "<<solution_method<<std::endl; }
  if(is_integral_method) {std::cout<<"integral_method = "<<integral_method<<std::endl; }
  if(is_epsilon)         {std::cout<<"epsilon = "<<epsilon<<std::endl;   }
  if(is_threshold)       {std::cout<<"threshold = "<<threshold<<std::endl;   }
  if(is_MaxCount)        {std::cout<<"MaxCount = "<<MaxCount<<std::endl;   }

  std::cout<<std::endl;

}

void Qeq::init_variables(){
   // Default parameters
  idempotential = "parr_pearson";       is_idempotential = 1;
  solution_method = 0;                  is_solution_method = 1;  // standard QEq methods (no dJad/aq derivatives)
  integral_method = 7;                  is_integral_method = 1;  // Louwen-Vogt equation
  epsilon =2.0;                         is_epsilon = 1;
  threshold = 0.001;                    is_threshold = 1;
  MaxCount =1000;                       is_MaxCount = 1;
}
 
void Qeq::copy_content(const Qeq& qeq){
  if(qeq.is_idempotential)  { idempotential   = qeq.idempotential;   is_idempotential = 1;}
  if(qeq.is_solution_method){ solution_method = qeq.solution_method; is_solution_method = 1;}
  if(qeq.is_integral_method){ integral_method = qeq.integral_method; is_integral_method = 1;}
  if(qeq.is_epsilon)        { epsilon         = qeq.epsilon;         is_epsilon = 1;}
  if(qeq.is_threshold)      { threshold       = qeq.threshold;       is_threshold = 1;}
  if(qeq.is_MaxCount)       { MaxCount        = qeq.MaxCount;        is_MaxCount = 1;}
}

Qeq::Qeq(){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();


}

Qeq::Qeq(const Qeq& qeq){
  /********************
    Copy constructor
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of qeq object which is defined
  copy_content(qeq);
}

Qeq& Qeq::operator=(const Qeq& qeq){
  /********************
    Assignment operator
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of qeq object which is defined
  copy_content(qeq);
  return *this;
}

Qeq::~Qeq(){;;}

 
