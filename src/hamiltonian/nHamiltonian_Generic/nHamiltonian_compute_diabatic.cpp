/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian_compute_diabatic.cpp
  \brief The file implements calculations of the diabatic Hamiltonians
    
*/


#include <stdlib.h>

#include "nHamiltonian.h"
#include "../Hamiltonian_Model/libhamiltonian_model.h"
#include "../../io/libio.h"

/// liblibra namespace
namespace liblibra{

/// libhamiltonian namespace 
namespace libhamiltonian{

using namespace libhamiltonian_model;
using namespace libio;

/// libhamiltonian_generic namespace 
namespace libhamiltonian_generic{


namespace bp = boost::python;


bp::object import_py_funct(const std::string& module, const std::string& path, bp::object& globals)
{
    bp::dict locals;
    locals["module_name"] = module;
    locals["path"]        = path;

    bp::exec("import imp\n"
             "new_module = imp.load_module(module_name, open(path), path, ('py', 'U', imp.PY_SOURCE))\n",
             globals,
             locals);
    return locals["new_module"];
}



void nHamiltonian::compute_diabatic(int model, vector<double>& q, vector<double>& params){

  

  if(model==100){   model_1S_1D_poly2(ham_dia, ovlp_dia, d1ham_dia, dc1_dia, q, params);  }
  if(model==100){   model_1S_1D_poly4(ham_dia, ovlp_dia, d1ham_dia, dc1_dia, q, params);  }

}



void nHamiltonian::compute_diabatic(bp::object py_funct, bp::object q, bp::object params){
/**
  The Python function should correspond to the following C++ signature:

  py_funct(CMATRIX& _Hdia, CMATRIX& _Sdia, vector<CMATRIX>& _d1ham_dia, 
           vector<CMATRIX>& _dc1_dia, vector<double>& q, bp::object params);

*/

  // Temporary storage
  CMATRIX _Hdia(ndia,ndia);
  CMATRIX _Sdia(ndia,ndia); 
  vector<CMATRIX> _d1ham_dia(nnucl, CMATRIX(ndia,ndia));
  vector<CMATRIX> _dc1_dia(nnucl, CMATRIX(ndia,ndia));

  // Call the Python function with such arguments
  bp::object obj = py_funct(q, params);  

  // Extract all the computed properties
  int has_attr=0;
  has_attr = (int)hasattr(obj,"ham_dia");        
  if(has_attr){    _Hdia = extract<CMATRIX>(obj.attr("ham_dia"));    }

  has_attr=0;
  has_attr = (int)hasattr(obj,"ovlp_dia");        
  if(has_attr){    _Sdia = extract<CMATRIX>(obj.attr("ovlp_dia"));    }

  has_attr=0;
  has_attr = (int)hasattr(obj,"d1ham_dia");        
  if(has_attr){    _d1ham_dia = extract<CMATRIXList>(obj.attr("d1ham_dia"));    }

  has_attr=0;
  has_attr = (int)hasattr(obj,"dc1_dia");        
  if(has_attr){    _dc1_dia = extract<CMATRIXList>(obj.attr("dc1_dia"));    }


  // Copy the temporary results into the Hamiltonian-bound memory
  *ham_dia = _Hdia;
  *ovlp_dia = _Sdia;
  
  for(int i=0;i<nnucl;i++){
    *d1ham_dia[i] = _d1ham_dia[i];
    *dc1_dia[i]   = _dc1_dia[i];
  }


}




}// namespace libhamiltonian_generic
}// namespace libhamiltonian
}// liblibra

