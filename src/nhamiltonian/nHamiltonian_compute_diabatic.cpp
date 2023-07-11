/*********************************************************************************
* Copyright (C) 2018-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian_compute_diabatic.cpp
  \brief The file implements calculations of the diabatic Hamiltonians
    
*/


#if defined(USING_PCH)
#include "../pch.h"
#else
#include <stdlib.h>
#include <omp.h>
#endif 

#include "nHamiltonian.h"
//#include "../Hamiltonian_Model/libhamiltonian_model.h"
#include "../io/libio.h"

/// liblibra namespace
namespace liblibra{

/// libnhamiltonian namespace 
namespace libnhamiltonian{

///using namespace libhamiltonian_model;
using namespace libio;


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

  compute_diabatic(model, q, params, 0);

}

void nHamiltonian::compute_diabatic(int model, vector<double>& q, vector<double>& params, int lvl){

  if(level==lvl){
/* 
    if(model==100){   model_1S_1D_poly2(ham_dia, ovlp_dia, d1ham_dia, dc1_dia, q, params);  }
    if(model==100){   model_1S_1D_poly4(ham_dia, ovlp_dia, d1ham_dia, dc1_dia, q, params);  }
*/
  }
  else if(lvl>level){
  
    for(int i=0;i<children.size();i++){
      children[i]->compute_diabatic(model, q, params, lvl);
    }

  }

  else{
    cout<<"WARNING in nHamiltonian::compute_diabatic\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }



}


//void nHamiltonian::compute_diabatic(bp::object py_funct, bp::object q, bp::object params){
void nHamiltonian::compute_diabatic(bp::object py_funct, MATRIX& q, bp::object params){
/**
  Performs the diabatic properties calculation at the top-most level of the Hamiltonians 
  hierarchy. See the description of the more general function prototype for more info.
*/ 

  compute_diabatic(py_funct, q, params, 0);

}



//void nHamiltonian::compute_diabatic(bp::object py_funct, bp::object q, bp::object params, int lvl){
void nHamiltonian::compute_diabatic(bp::object py_funct, MATRIX& q, bp::object params, int lvl){
/**
  This function will call the <py_funct> function defined in Python and taking the signature:

  def py_funct(q, params, full_id)
      return obj

  The object <obj> returned by the function should contain variables named:
  "ham_dia" (CMATRIX), "ovld_dia" (CMATRIX), "d1ham_dia" (list of CMATRIX), "dc1_dia" (list of CMATRIX)
  of the specified types

  q - is the object containing coordinates of nuclei. Since almost any C++ type exposed to Python is seen 
  by Python as a Python object, you can use any convenient Libra type to store the coordinates, for instance
  a MATRIX object. 

  params - is an object containing any parameters needed by the <py_funct> Python function to perform
  the specified calculations. Again, there are many options to pass the parameters on the Python side, but 
  it may be convenient to use Python dictionary, since the parameters will be mostly used at the 
  Python level of calculations.

  full_id - is vector<int> object containing the "path" of the calling node in the hierarchy of Hamiltonians
  this information may be needed by the py_funct function in order to do the proper calculations

  lvl - is the level of the Hamiltonians in the hierarchy of Hamiltonians to be executed by this call 
  You can only request to perform computations at the same of higher (children) levels that the level of 
  the Hamiltonian from which you call this function. 

  Contentwise, this function computes the diabatic properties and populates the corresponding
  storage.

*/


  if(level==lvl){

    // Temporary storage


  
    // Call the Python function with such arguments
    bp::object obj = py_funct(q, params, get_full_id() );  
  
    // Extract all the computed properties
    int has_attr=0;
    has_attr = (int)hasattr(obj,"ham_dia");        
    if(has_attr){  

      check_cmatrix(obj, "ham_dia", ndia, ndia);
      *ham_dia = extract<CMATRIX>(obj.attr("ham_dia")); 
    }
  
    has_attr=0;
    has_attr = (int)hasattr(obj,"ovlp_dia");        
    if(has_attr){ 

      check_cmatrix(obj, "ovlp_dia", ndia, ndia);
      *ovlp_dia = extract<CMATRIX>(obj.attr("ovlp_dia"));   
    }


    has_attr=0;
    has_attr = (int)hasattr(obj,"nac_dia");        
    if(has_attr){    

      check_cmatrix(obj, "nac_dia", ndia, ndia);
      *nac_dia = extract<CMATRIX>(obj.attr("nac_dia"));    
    }

    has_attr=0;
    has_attr = (int)hasattr(obj,"hvib_dia");        
    if(has_attr){    

      check_cmatrix(obj, "hvib_dia", ndia, ndia);
      *hvib_dia = extract<CMATRIX>(obj.attr("hvib_dia"));    
    }
    
    has_attr=0;
    has_attr = (int)hasattr(obj,"dc1_dia");        
    if(has_attr){  

      check_cmatrix_list(obj, "dc1_dia", ndia, ndia, nnucl);

      vector<CMATRIX> _dc1_dia(nnucl, CMATRIX(ndia,ndia));
      _dc1_dia = extract<CMATRIXList>(obj.attr("dc1_dia"));   

      for(int i=0;i<nnucl;i++){   *dc1_dia[i]   = _dc1_dia[i];    }
    }


    has_attr=0;
    has_attr = (int)hasattr(obj,"d1ham_dia");        
    if(has_attr){  

      check_cmatrix_list(obj, "d1ham_dia", ndia, ndia, nnucl);

      vector<CMATRIX> _d1ham_dia(nnucl, CMATRIX(ndia,ndia));
      _d1ham_dia = extract<CMATRIXList>(obj.attr("d1ham_dia"));   

      for(int i=0;i<nnucl;i++){    *d1ham_dia[i] = _d1ham_dia[i];  }
    }


    has_attr=0;
    has_attr = (int)hasattr(obj,"d2ham_dia");        
    if(has_attr){  

      check_cmatrix_list(obj, "d2ham_dia", ndia, ndia, nnucl*nnucl);

      vector<CMATRIX> _d2ham_dia(nnucl*nnucl, CMATRIX(ndia,ndia));
      _d2ham_dia = extract<CMATRIXList>(obj.attr("d2ham_dia"));   

      for(int i=0;i<nnucl;i++){    *d2ham_dia[i] = _d2ham_dia[i];  }
    }

    has_attr=0;
    has_attr = (int)hasattr(obj,"time_overlap_dia"); 
    if(has_attr){    

      check_cmatrix(obj, "time_overlap_adi", ndia, ndia);
      *time_overlap_dia = extract<CMATRIX>(obj.attr("time_overlap_dia"));    
    }


  
    
  }// if lvl == level

  else if(lvl>level){

    int nthreads = omp_get_num_threads();
    #pragma omp parallel num_threads( nthreads )
    {
        #pragma omp for    
        for(int i=0;i<children.size();i++){
          children[i]->compute_diabatic(py_funct,q,params,lvl);
        }
    }// pragma

  }

  else{
    cout<<"WARNING in nHamiltonian::compute_diabatic\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }



}




}// namespace libnhamiltonian
}// liblibra

