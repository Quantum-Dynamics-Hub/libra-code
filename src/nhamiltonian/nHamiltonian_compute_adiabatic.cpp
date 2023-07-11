/*********************************************************************************
* Copyright (C) 2017-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian_compute_adiabatic.cpp
  \brief The file implements the calculations of the adiabatic Hamiltonian (e.g. as a transformation)
    
*/


#if defined(USING_PCH)
#include "../pch.h"
#else
#include <stdlib.h>
#include <omp.h>
#endif 

#include "nHamiltonian.h"
#include "../math_meigen/libmeigen.h"


/// liblibra namespace
namespace liblibra{

/// libnhamiltonian namespace 
namespace libnhamiltonian{


using namespace liblinalg;
using namespace libmeigen;



void nHamiltonian::update_ordering(vector<int>& perm_t, int lvl){

  int i;

  if(level==lvl){

    update_permutation(perm_t, ordering_adi);
    //ordering_adi[0] = perm_t;



    // Apply the re-ordering to all the underlying properties:
    if(basis_transform_mem_status){
      basis_transform->permute_cols(perm_t);
    }
    
    if(ham_adi_mem_status){
      ham_adi->permute_cols(perm_t);
      ham_adi->permute_rows(perm_t);
    }
    
    if(nac_adi_mem_status){
      nac_adi->permute_cols(perm_t);
      nac_adi->permute_rows(perm_t);
    }
    
    if(hvib_adi_mem_status){
      hvib_adi->permute_cols(perm_t);
      hvib_adi->permute_rows(perm_t);
    }
    
    
    for(i=0; i<dc1_adi.size(); i++){
      if(dc1_adi_mem_status[i]){
        dc1_adi[i]->permute_cols(perm_t);
        dc1_adi[i]->permute_rows(perm_t);
      }
    }
    
    for(i=0; i<d1ham_adi.size(); i++){
      if(d1ham_adi_mem_status[i]){
        d1ham_adi[i]->permute_cols(perm_t);
        d1ham_adi[i]->permute_rows(perm_t);
      }
    }
    
    for(i=0; i<d2ham_adi.size(); i++){
      if(d2ham_adi_mem_status[i]){
        d2ham_adi[i]->permute_cols(perm_t);
        d2ham_adi[i]->permute_rows(perm_t);
      }
    }


  }// level==lvl

  else if(lvl>level){
  
    for(int i=0;i<children.size();i++){
      children[i]->update_ordering(perm_t, lvl);
    }

  }// lvl >level

  else{
    cout<<"WARNING in void nHamiltonian::update_order(vector<int>& perm_t, int lvl):\n"; 
    cout<<"Can not run the function in the parent Hamiltonian from the\
     child node\n";    
  }
 

}



void nHamiltonian::update_ordering(vector<int>& perm_t){

  update_ordering(perm_t, 0);

}


/**

void nHamiltonian::apply_phase_corrections(CMATRIX* phase_corr, int lvl){

//  phase_corr - CMATRIX(nadi, 1) - the changes of cumulative phases of all eigenvectors. 
  


  int i, j;
  complex<double> fji;

  if(level==lvl){

    if(basis_transform_mem_status){
      if(phase_corr->n_rows!=basis_transform->n_cols){
        cout<<"ERROR in void nHamiltonian::update_phases(vector<double>& phases, int lvl): \
        the number of elements in the <phases> input should be equal to the number of columns of \
        the basis transformation matrix\nExiting...\n";
        exit(0);
      }

      for(i=0; i<basis_transform->n_cols; i++){
        for(j=0; j<basis_transform->n_rows; j++){
          basis_transform->scale(j,i, std::conj(phase_corr->get(i,0)) );
        }
      }

    }// basis_transform


    if(nac_adi_mem_status){
      if(phase_corr->n_rows!=nac_adi->n_cols){
        cout<<"ERROR in void nHamiltonian::update_phases(vector<double>& phases, int lvl): \
        the number of elements in the <phases> input should be equal to the number of columns of \
        the nac_adi matrix\nExiting...\n";
        exit(0);
      }

      for(i=0; i<nac_adi->n_cols; i++){
        for(j=0; j<nac_adi->n_rows; j++){          
          fji = phase_corr->get(j,0) * std::conj(phase_corr->get(i,0));
          nac_adi->scale(j,i, fji);
        }
      }
    }// nac_adi


    if(hvib_adi_mem_status){
      if(phase_corr->n_rows!=hvib_adi->n_cols){
        cout<<"ERROR in void nHamiltonian::update_phases(vector<double>& phases, int lvl): \
        the number of elements in the <phases> input should be equal to the number of columns of \
        the hvib_adi matrix\nExiting...\n";
        exit(0);
      }

      for(i=0; i<hvib_adi->n_cols; i++){
        for(j=0; j<hvib_adi->n_rows; j++){          
          fji = phase_corr->get(j,0) * std::conj(phase_corr->get(i,0));
          hvib_adi->scale(j,i, fji);
        }
      }
    }// nac_adi

  }// level==lvl

  else if(lvl>level){
  
    for(int i=0;i<children.size();i++){
      children[i]->apply_phase_corrections(phase_corr, lvl);
    }

  }// lvl >level

  else{
    cout<<"WARNING in void nHamiltonian::update_phases(vector<double>& phases, int lvl) :\n"; 
    cout<<"Can not run the function in the parent Hamiltonian from the\
     child node\n";    
  }

}

void nHamiltonian::apply_phase_corrections(CMATRIX& phase_corr, int lvl){

  apply_phase_corrections(&phase_corr, lvl);
}


void nHamiltonian::apply_phase_corrections(CMATRIX* phase_corr){

  apply_phase_corrections(phase_corr, 0);
}

void nHamiltonian::apply_phase_corrections(CMATRIX& phase_corr){

  apply_phase_corrections(&phase_corr, 0);
}


*/


/**

CMATRIX nHamiltonian::update_phases(CMATRIX& U_prev, int lvl){

//  This function computes the phase corrections to all current eigenvectors
//  w.r.t. those in U_prev.
//  It also applies the correction to the present eigenvectors. 
//  Finally, this function computes the phases by which the adiabatic amplitudes 
//  should be updated and returns these corrections as a column-vector.
//  
//  basis_transform - the |psi(t')>
//  U_prev = |psi(t)>
//
//  t' > t



  if(lvl==level){

    // Memorize the cumulative phase correction up to this step
    CMATRIX* cum_phase_corr_prev; 
    cum_phase_corr_prev = new CMATRIX(nadi, 1);
    *cum_phase_corr_prev = *cum_phase_corr;
  
    // Compute cumulative phase corrections
    *cum_phase_corr = compute_phase_corrections1(*basis_transform, U_prev, phase_corr_ovlp_tol);

    // Update the wavefunction and wavefunction-dependent properties with 
    // the new phase corrections    
    apply_phase_corrections(cum_phase_corr, lvl);

    // Compute the phase correction to the evolving amplitudes
    CMATRIX ampl_corr(nadi, 1);

    for(int i=0;i<nadi; i++){
      complex<double> scl(1.0, 0.0);

      if( abs(cum_phase_corr_prev->get(i,0) ) > 0.0 ){
        scl = cum_phase_corr->get(i,0)/cum_phase_corr_prev->get(i,0); 
      }
      ampl_corr.set(i,0, scl);
    }
    delete cum_phase_corr_prev;
 
    return ampl_corr;

  }// level==lvl

  else if(lvl>level){
  
    for(int i=0;i<children.size();i++){
      children[i]->apply_phase_corrections(U_prev, lvl);
    }

  }// lvl >level

  else{
    cout<<"WARNING in void nHamiltonian::update_phases(CMATRIX& U_prev, int lvl) :\n"; 
    cout<<"Can not run the function in the parent Hamiltonian from the child node\n";    
  }
}

CMATRIX nHamiltonian::update_phases(CMATRIX& U_prev){

  return update_phases(U_prev, 0);
}
*/



void nHamiltonian::compute_adiabatic(int der_lvl){

  compute_adiabatic(der_lvl, 0);

}

void nHamiltonian::compute_adiabatic(int der_lvl, int lvl){
/**
  Compute the adiabatic Hamiltonian

  |psi_adi> = |psi_dia> * U

  Stationary SE:  H |psi_adi> =  |psi_adi> * E

  project on <psi_adi|, keeping in mind that <psi_adi|psi_adi> = I

  <psi_dia| H |psi_adi> =  <psi_dia|psi_adi> E  

  <psi_dia| H |psi_dia> * U =  <psi_dia|psi_dia> * U * E  

  H_dia * U = ovlp_dia * U * H_adi

  Here, U = basis_transform

  Assume: H does not depend on U

  der_lvl >= 0 - only diabatic - to - adiabatic transform
  der_lvl >= 1 - forces and derivative couplings
  der_lvl >= 2 - Hessian

*/

  int i,j;

  if(level==lvl){
 
  if(der_lvl>=0){

    if(ham_dia_mem_status==0){ cout<<"Error in compute_adiabatic(): the diabatic Hamiltonian matrix is not allocated \
    but it is needed for the calculations\n"; exit(0); }
    if(ovlp_dia_mem_status==0){ cout<<"Error in compute_adiabatic(): the overlap matrix of the diabatic states is not allocated \
    but it is needed for the calculations\n"; exit(0); }

    if(ham_adi_mem_status==0){ cout<<"Error in compute_adiabatic(): the adiabatic Hamiltonian matrix is not allocated \
    but it is used to collect the results of the calculations\n"; exit(0); }
    if(basis_transform_mem_status==0){ cout<<"Error in compute_adiabatic(): the basis_transform (eigenvector) matrix is\
    not allocated but it is used to collect the results of the calculations\n"; exit(0); }

    eigen_algo = 0;  /// Should be this option, otherwise the static calculations are incorrect

    if(nadi==1 && ndia==1){
      *ham_adi = *ham_dia;
      basis_transform->set(0,0, 1.0, 0.0);
    }
    else{   
      if(eigen_algo==0){
          solve_eigen(ham_dia, ovlp_dia, ham_adi, basis_transform, 1);  
      }
      else if(eigen_algo==1){
          solve_eigen_nosort(ham_dia, ham_adi, basis_transform, 1);       ///< references
      }
    }

    if(der_lvl>=1){

      // Now compute the derivative couplings (off-diagonal, multiplied by energy difference) and adiabatic gradients (diagonal)
      CMATRIX* tmp; tmp = new CMATRIX(nadi,nadi);
      CMATRIX* dtilda; dtilda = new CMATRIX(nadi,nadi);
      
      for(int n=0;n<nnucl;n++){

        if(d1ham_dia_mem_status[n]==0){ cout<<"Error in compute_adiabatic(): the derivatives of the diabatic Hamiltonian \
        matrix w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for calculations \n"; exit(0); }

        if(d1ham_adi_mem_status[n]==0){ cout<<"Error in compute_adiabatic(): the derivatives of the adiabatic Hamiltonian \
        matrix w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for collecting results \n"; exit(0); }


        if(dc1_dia_mem_status[n]==0){ cout<<"Error in compute_adiabatic(): the derivatives couplings matrix in the diabatic \
        basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for collecting results \n"; exit(0); }

        if(dc1_adi_mem_status[n]==0){ cout<<"Error in compute_adiabatic(): the derivatives couplings matrix in the adiabatic \
        basis w.r.t. the nuclear DOF "<<n<<" is not allocated but is needed for collecting results \n"; exit(0); }


        // E.g. see the derivations here: https://github.com/alexvakimov/Derivatory/blob/master/theory_NAC.pdf
        // also: http://www.theochem.ruhr-uni-bochum.de/~nikos.doltsinis/nic_10_doltsinis.pdf
        *tmp = (*basis_transform).H() * (*d1ham_dia[n]) * (*basis_transform);

//        *dtilda = (*basis_transform) * (*dc1_dia[n]).H() * (*basis_transform).H() * (*ham_adi);
        *dtilda = (*basis_transform).H() * (*dc1_dia[n]) * (*basis_transform) * (*ham_adi);
        *dtilda = (*dtilda + (*dtilda).H() );

        *tmp -= *dtilda;

        // Adiabatic "forces"
        *d1ham_adi[n] = 0.0;
        for(i=0;i<nadi;i++){ d1ham_adi[n]->set(i,i, tmp->get(i,i)); }

        // Adiabatic derivative couplings
        *dc1_adi[n] = 0.0;

        for(i=0;i<nadi;i++){
          for(j=i+1;j<nadi;j++){

            //if(i==j){  dc1_adi[n]->set(i,j, 0.0, 0.0); }
            //else{
 
              double dE = (ham_adi->get(j,j) - ham_adi->get(i,i) ).real();

              if(fabs(dE)<1e-25){ dE = 1e-25; }
              complex<double> val = tmp->get(i,j)/dE;
              
              dc1_adi[n]->set(i,j, val);
              dc1_adi[n]->set(j,i,-val);


              /*
              if(fabs(dE)<1e-100){ 
                //dE = 1e-10 * (dE>0.0 ? 1.0 : -1.0);                 
                dE = 1e-100;
                dc1_adi[n]->set(i,j, tmp->get(i,j)/dE );
                //dc1_adi[n]->set(i,j, 0.0, 0.0 );
              }else{
          
                dc1_adi[n]->set(i,j, tmp->get(i,j)/dE );
              }
              */

            //}

          }// for j
        }// for i
      }// for n

      delete tmp;
      delete dtilda;

    }// der_lvl>=1
  }// der_lvl>=0

  }// level == lvl

  else if(lvl>level){
  
    for(int i=0;i<children.size();i++){
      children[i]->compute_adiabatic(der_lvl, lvl);
    }

  }// lvl >level

  else{
    cout<<"WARNING in nHamiltonian::compute_adiabatic\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }
     
}



//void nHamiltonian::compute_adiabatic(bp::object py_funct, bp::object q, bp::object params){
void nHamiltonian::compute_adiabatic(bp::object py_funct, MATRIX& q, bp::object params){
/**
  Performs the adiabatic properties calculation at the top-most level of the Hamiltonians 
  hierarchy. See the description of the more general function prototype for more info.
*/ 

  compute_adiabatic(py_funct, q, params, 0);

}


//void nHamiltonian::compute_adiabatic(bp::object py_funct, bp::object q, bp::object params, int lvl){
void nHamiltonian::compute_adiabatic(bp::object py_funct, MATRIX& q, bp::object params, int lvl){
/**
  This function will call the <py_funct> function defined in Python and taking the signature:

  def py_funct(q, params, full_id)
      return obj

  The object <obj> returned by the function should contain variables named:
  "ham_adi" (CMATRIX), "d1ham_adi" (list of CMATRIX), "dc1_adi" (list of CMATRIX)
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

  Contentwise, this function computes the adiabatic properties and populates the corresponding
  storage.

*/

  if(level==lvl){

    // Call the Python function with such arguments
    bp::object obj = py_funct(q, params, get_full_id() );  

 
    // Extract all the computed properties
    int has_attr=0;
    has_attr = (int)hasattr(obj,"ham_adi");        
    if(has_attr){    

      check_cmatrix(obj, "ham_adi", nadi, nadi);
      *ham_adi = extract<CMATRIX>(obj.attr("ham_adi"));    
    }

    has_attr=0;
    has_attr = (int)hasattr(obj,"nac_adi");        
    if(has_attr){    

      check_cmatrix(obj, "nac_adi", nadi, nadi);
      *nac_adi = extract<CMATRIX>(obj.attr("nac_adi"));    
    }

    has_attr=0;
    has_attr = (int)hasattr(obj,"hvib_adi");        
    if(has_attr){    

      check_cmatrix(obj, "hvib_adi", nadi, nadi);
      *hvib_adi = extract<CMATRIX>(obj.attr("hvib_adi"));    
    }

    has_attr=0;
    has_attr = (int)hasattr(obj,"basis_transform");        
    if(has_attr){    

      check_cmatrix(obj, "basis_transform", nadi, nadi);
      *basis_transform = extract<CMATRIX>(obj.attr("basis_transform"));    
    }

    has_attr=0;
    has_attr = (int)hasattr(obj,"time_overlap_adi"); 
    if(has_attr){    

      check_cmatrix(obj, "time_overlap_adi", nadi, nadi);
      *time_overlap_adi = extract<CMATRIX>(obj.attr("time_overlap_adi"));    
    }


    has_attr=0;
    has_attr = (int)hasattr(obj,"dc1_adi");        
    if(has_attr){

      check_cmatrix_list(obj, "dc1_adi", nadi, nadi, nnucl);

      vector<CMATRIX> _dc1_adi(nnucl, CMATRIX(nadi,nadi));
      _dc1_adi = extract<CMATRIXList>(obj.attr("dc1_adi"));    

      for(int i=0;i<nnucl;i++){   *dc1_adi[i]   = _dc1_adi[i];     }
    }

  
    has_attr=0;
    has_attr = (int)hasattr(obj,"d1ham_adi");        
    if(has_attr){

      check_cmatrix_list(obj, "d1ham_adi", nadi, nadi, nnucl);

      vector<CMATRIX> _d1ham_adi(nnucl, CMATRIX(nadi,nadi));
      _d1ham_adi = extract<CMATRIXList>(obj.attr("d1ham_adi"));    

      for(int i=0;i<nnucl;i++){   *d1ham_adi[i] = _d1ham_adi[i];   }
    }


    has_attr=0;
    has_attr = (int)hasattr(obj,"d2ham_adi");        
    if(has_attr){

      check_cmatrix_list(obj, "d2ham_adi", nadi, nadi, nnucl*nnucl);

      vector<CMATRIX> _d2ham_adi(nnucl*nnucl, CMATRIX(nadi,nadi));
      _d2ham_adi = extract<CMATRIXList>(obj.attr("d2ham_adi"));    

      for(int i=0;i<nnucl*nnucl;i++){   *d2ham_adi[i] = _d2ham_adi[i];   }
    }
  
  
  
  
  }// if lvl == level

  else if(lvl>level){

    int nthreads = omp_get_num_threads();
    #pragma omp parallel num_threads( nthreads )
    {
        #pragma omp for  
        for(int i=0;i<children.size();i++){
          children[i]->compute_adiabatic(py_funct,q,params,lvl);
        }
    }// pragma

  }

  else{
    cout<<"WARNING in nHamiltonian::compute_adiabatic\n"; 
    cout<<"Can not run evaluation of function in the parent Hamiltonian from the\
     child node\n";    
  }


}



}// namespace libnhamiltonian
}// liblibra

