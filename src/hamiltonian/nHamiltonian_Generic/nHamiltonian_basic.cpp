/*********************************************************************************
* Copyright (C) 2017-2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian_basic.cpp
  \brief The file implements the generic methods of the nHamiltonian class: getters and setters
    
*/

#include "nHamiltonian.h"
#include <stdlib.h>


/// liblibra namespace
namespace liblibra{

/// libhamiltonian namespace 
namespace libhamiltonian{

/// libhamiltonian_generic namespace 
namespace libhamiltonian_generic{

/*
nHamiltonian::nHamiltonian(){ 
// Default constructor of base (Hamiltonian) class


  level = 0;
  id = 0;

//  parent = NULL;

  ovlp_dia = NULL;             ovlp_dia_mem_status = 0; 

  ham_dia = NULL;              ham_dia_mem_status = 0;
  nac_dia = NULL;              nac_dia_mem_status = 0;
  hvib_dia = NULL;             hvib_dia_mem_status = 0;

  ham_adi = NULL;              ham_adi_mem_status = 0;
  nac_adi = NULL;              nac_adi_mem_status = 0;
  hvib_adi = NULL;             hvib_adi_mem_status = 0;

  basis_transform = NULL;      basis_transform_mem_status = 0;

}
*/

nHamiltonian::nHamiltonian(int ndia_, int nadi_, int nnucl_){ 
/** Constructor of the nHamiltonian class
*/

  cout<<"nHamiltonian constructor at"<<this<<endl;

//  nHamiltonian();

  int n;

  level = 0;
  id = 0;

//  parent = NULL;

  ovlp_dia = NULL;             ovlp_dia_mem_status = 0; 

  ham_dia = NULL;              ham_dia_mem_status = 0;
  nac_dia = NULL;              nac_dia_mem_status = 0;
  hvib_dia = NULL;             hvib_dia_mem_status = 0;

  ham_adi = NULL;              ham_adi_mem_status = 0;
  nac_adi = NULL;              nac_adi_mem_status = 0;
  hvib_adi = NULL;             hvib_adi_mem_status = 0;

  basis_transform = NULL;      basis_transform_mem_status = 0;

  cum_phase_corr = NULL;       cum_phase_corr_mem_status = 0;


  ndia = ndia_;                   
  nadi = nadi_;
  nnucl = nnucl_;

  ordering_adi = new vector<int>(nadi, 0);
  for(n=0;n<nadi;n++){ ordering_adi[0][n] = n; }

  dc1_dia = vector<CMATRIX*>(nnucl, NULL); 
  dc1_dia_mem_status = vector<int>(nnucl, 0);

  d1ham_dia = vector<CMATRIX*>(nnucl, NULL);
  d1ham_dia_mem_status = vector<int>(nnucl, 0);

  d2ham_dia = vector<CMATRIX*>(nnucl*nnucl, NULL);
  d2ham_dia_mem_status = vector<int>(nnucl*nnucl, 0);

  dc1_adi = vector<CMATRIX*>(nnucl, NULL);
  dc1_adi_mem_status = vector<int>(nnucl, 0);

  d1ham_adi = vector<CMATRIX*>(nnucl, NULL);
  d1ham_adi_mem_status = vector<int>(nnucl, 0);

  d2ham_adi = vector<CMATRIX*>(nnucl*nnucl, NULL);
  d2ham_adi_mem_status = vector<int>(nnucl*nnucl, 0);


}



nHamiltonian::~nHamiltonian(){ 
/**
  Deallocate memory only if it was allocated internally
*/

  cout<<"nHamiltonian destructor"<<this<<endl;
  
//  id = 0;
//  level = 0;
  int n;
 
  delete ordering_adi;  ordering_adi = NULL;
 
  if(ovlp_dia_mem_status == 1){ delete ovlp_dia;  ovlp_dia = NULL; ovlp_dia_mem_status = 0;}

  for(n;n<dc1_dia.size();n++){
    if(dc1_dia_mem_status[n] == 1){ delete dc1_dia[n];  dc1_dia[n] = NULL; dc1_dia_mem_status[n] = 0;}
  } 
  dc1_dia.clear();
  dc1_dia_mem_status.clear();

  if(ham_dia_mem_status == 1){ delete ham_dia; ham_dia = NULL; ham_dia_mem_status = 0;}
  if(nac_dia_mem_status == 1){ delete nac_dia; nac_dia = NULL; nac_dia_mem_status = 0;}
  if(hvib_dia_mem_status == 1){ delete hvib_dia; hvib_dia = NULL; hvib_dia_mem_status = 0;}

  for(n;n<d1ham_dia.size();n++){
    if(d1ham_dia_mem_status[n] == 1){ delete d1ham_dia[n];  d1ham_dia[n] = NULL; d1ham_dia_mem_status[n] = 0;}
  } 
  d1ham_dia.clear();
  d1ham_dia_mem_status.clear();

  for(n;n<d2ham_dia.size();n++){
    if(d2ham_dia_mem_status[n] == 1){ delete d2ham_dia[n];  d2ham_dia[n] = NULL; d2ham_dia_mem_status[n] = 0;}
  } 
  d2ham_dia.clear();
  d2ham_dia_mem_status.clear();


  for(n;n<dc1_adi.size();n++){
    if(dc1_adi_mem_status[n] == 1){ delete dc1_adi[n];  dc1_adi[n] = NULL; dc1_adi_mem_status[n] = 0;}
  } 
  dc1_adi.clear();
  dc1_adi_mem_status.clear();


  if(ham_adi_mem_status == 1){ delete ham_adi; ham_adi = NULL; ham_adi_mem_status = 0; }
  if(nac_adi_mem_status == 1){ delete nac_adi; nac_adi = NULL; nac_adi_mem_status = 0;}
  if(hvib_adi_mem_status == 1){ delete hvib_adi; hvib_adi = NULL; hvib_adi_mem_status = 0;}


  for(n;n<d1ham_adi.size();n++){
    if(d1ham_adi_mem_status[n] == 1){ delete d1ham_adi[n];  d1ham_adi[n] = NULL; d1ham_adi_mem_status[n] = 0;}
  } 
  d1ham_adi.clear();
  d1ham_adi_mem_status.clear();

  for(n;n<d2ham_adi.size();n++){
    if(d2ham_adi_mem_status[n] == 1){ delete d2ham_adi[n];  d2ham_adi[n] = NULL; d2ham_adi_mem_status[n] = 0;}
  } 
  d2ham_adi.clear();
  d2ham_adi_mem_status.clear();

  if(basis_transform_mem_status == 1){ delete basis_transform; basis_transform = NULL; basis_transform_mem_status = 0;}

  if(cum_phase_corr_mem_status == 1){ delete cum_phase_corr; cum_phase_corr = NULL; cum_phase_corr_mem_status = 0;}

//  if(next!=NULL){
//    for(n=0;n<next.size();n++){  next[n]->~nHamiltonian(); }
//  }
//  next = NULL;

}


void nHamiltonian::init_all(int lvl){ 
/**
  This function will allocate memory for all the variables - in case, we don't need to 
  set them up in Python. This may lead to a simpler and cleaner Python code and may
  also be a bit more efficient in terms of performance
*/

  if(lvl>=0){
    init_ovlp_dia();
    init_ham_dia();
    init_nac_dia();
    init_hvib_dia();

    init_ham_adi();
    init_nac_adi();
    init_hvib_adi();

    init_basis_transform();
    init_cum_phase_corr();
  }

  if(lvl>=1){
    init_dc1_dia();
    init_d1ham_dia();

    init_dc1_adi();
    init_d1ham_adi();
  }

  if(lvl>=2){
    init_d2ham_dia();
    init_d2ham_adi();
  }


}


void nHamiltonian::set_levels(int lvl_){
/**
  This function set the level of a given Hamiltonian to the specified value
  and propagates this value to all children of this Hamiltonian
*/
  cout<<"in set_levels with lvl = "<<lvl_<<endl;
  level = lvl_;

  for(int i=0; i<children.size(); i++){ 
    cout<<" == child "<<i<<endl;
    children[i]->set_levels(lvl_+1);
  }

}


void nHamiltonian::add_child(nHamiltonian& child){
/**
  Associate an external Hamiltonian <child> with the present Hamiltonian 
  such that the <child> becomes connected with <this> by a child-parent relationship.
*/

  child.id = children.size();
  cout<<"in add_child, new child id = "<<child.id<<"\n";
  children.push_back(&child);

  child.parent = this;
  cout<<"parent = "<<child.parent<<endl;
  child.set_levels(level+1);

}

vector<int> nHamiltonian::get_full_id(){

  vector<int> res(1,0);

  if(level>0){
    res = parent->get_full_id();
    res.push_back(id);
  }

  return res;

}

/*************************************************************************

          SETTERS  :         DIABATIC PARAMETERS

**************************************************************************/

/************************ OVERLAPS *****************************/

void nHamiltonian::init_ovlp_dia(){
/**
  Allocate memory for the overlap matrix in the diabatic basis
*/
  if(ovlp_dia_mem_status==0){
    ovlp_dia = new CMATRIX(ndia, ndia);
    ovlp_dia_mem_status = 1;
  }
  else{ 
    cout<<"WARNING in init_ovlp_dia: memory is already allocated\n";
  }

}

void nHamiltonian::set_ovlp_dia_by_ref(CMATRIX& ovlp_dia_){
/**
  Setup of the overlap matrix in the diabatic basis
*/

//  set_X1_by_ref(ovlp_dia, ovlp_dia_, ovlp_dia_mem_status, ndia, ndia);

  check_mat_dimensions(&ovlp_dia_, ndia, ndia);

  if(ovlp_dia_mem_status==0){ ovlp_dia = &ovlp_dia_; } /// Not allocated - not a problem
  else if(ovlp_dia_mem_status==1){ delete ovlp_dia; ovlp_dia = &ovlp_dia_;   }
  else if(ovlp_dia_mem_status==2){ ovlp_dia = &ovlp_dia_; }

  ovlp_dia_mem_status = 2; // Allocated externally



}

void nHamiltonian::set_ovlp_dia_by_val(CMATRIX& ovlp_dia_){
/**
  Setup of the overlap matrix in the diabatic basis
*/

//  set_X1_by_val(ovlp_dia, ovlp_dia_, ovlp_dia_mem_status, ndia, ndia);

  check_mat_dimensions(&ovlp_dia_,ndia,ndia);

  if(ovlp_dia_mem_status==0){  ovlp_dia = new CMATRIX(ndia,ndia);  } 
  else if(ovlp_dia_mem_status==1){ check_mat_dimensions(ovlp_dia,ndia,ndia);  }
  else if(ovlp_dia_mem_status==2){  ovlp_dia = new CMATRIX(ndia,ndia);    }

  *ovlp_dia = ovlp_dia_;
  ovlp_dia_mem_status = 1; // Allocated internally

}


/************************ DERIVATIVE COUPLINGS *****************************/

void nHamiltonian::init_dc1_dia(){
/**
  Allocate memory for the derivative couplings in the diabatic basis
*/

  for(int n=0;n<nnucl;n++){

    if(dc1_dia_mem_status[n]==0){
      dc1_dia[n] = new CMATRIX(ndia, ndia); 
      dc1_dia_mem_status[n] = 1;
    }
    else{ 
      cout<<"WARNING in init_dc1_dia: memory for element"<< n <<" is already allocated\n";
    }
  }

}


void nHamiltonian::set_dc1_dia_by_ref(vector<CMATRIX>& dc1_dia_){
/**
  Setup of the derivative coupling matrices in the diabatic basis w.r.t all nuclear DOFs
*/

//  set_X2_by_ref(dc1_dia, dc1_dia_, dc1_dia_mem_status, ndia, ndia, nnucl);

  check_vector_dimensions(dc1_dia, dc1_dia_, dc1_dia_mem_status, nnucl);

  for(int n=0;n<nnucl;n++){
    check_mat_dimensions(&dc1_dia_[n], ndia, ndia);

    if(dc1_dia_mem_status[n]==0){ dc1_dia[n] = &dc1_dia_[n]; } 
    else if(dc1_dia_mem_status[n]==1){ delete dc1_dia[n]; dc1_dia[n] = &dc1_dia_[n];   }
    else if(dc1_dia_mem_status[n]==2){ dc1_dia[n] = &dc1_dia_[n]; }

    dc1_dia_mem_status[n] = 2; // Allocated externally
  }

}

void nHamiltonian::set_dc1_dia_by_val(vector<CMATRIX>& dc1_dia_){
/**
  Setup of the derivative coupling matrices in the diabatic basis w.r.t all nuclear DOFs
*/
//  set_X2_by_val(dc1_dia, dc1_dia_, dc1_dia_mem_status, ndia, ndia, nnucl);

  check_vector_dimensions(dc1_dia, dc1_dia_, dc1_dia_mem_status, nnucl);

  for(int n=0;n<nnucl;n++){

    check_mat_dimensions(&dc1_dia_[n], ndia, ndia);

    if(dc1_dia_mem_status[n]==0){  dc1_dia[n] = new CMATRIX(ndia,ndia);  } 
    else if(dc1_dia_mem_status[n]==1){ check_mat_dimensions(dc1_dia[n],ndia,ndia);  }
    else if(dc1_dia_mem_status[n]==2){ dc1_dia[n] = new CMATRIX(ndia,ndia);    }

    *dc1_dia[n] = dc1_dia_[n];
    dc1_dia_mem_status[n] = 1; // Allocated internally
  }

}

/************************ HAMILTONIAN *****************************/
void nHamiltonian::init_ham_dia(){
/**
  Allocate memory for the Hamiltonian matrix in the diabatic basis
*/

  if(ham_dia_mem_status==0){
    ham_dia = new CMATRIX(ndia, ndia);
    ham_dia_mem_status = 1;
  }
  else{ 
    cout<<"WARNING in init_ham_dia: memory is already allocated\n";
  }

}


void nHamiltonian::set_ham_dia_by_ref(CMATRIX& ham_dia_){
/**
  Setup of the Hamiltonian matrix in the diabatic basis
*/

  check_mat_dimensions(&ham_dia_, ndia, ndia);

  if(ham_dia_mem_status==0){ ham_dia = &ham_dia_; } /// Not allocated - not a problem
  else if(ham_dia_mem_status==1){ delete ham_dia; ham_dia = &ham_dia_;   }
  else if(ham_dia_mem_status==2){ ham_dia = &ham_dia_; }

  ham_dia_mem_status = 2; // Allocated externally


}

void nHamiltonian::set_ham_dia_by_val(CMATRIX& ham_dia_){
/**
  Setup of the Hamiltonian matrix in the diabatic basis
*/

//  set_X1_by_val(ham_dia, ham_dia_, ham_dia_mem_status, ndia, ndia);

  check_mat_dimensions(&ham_dia_,ndia,ndia);

  if(ham_dia_mem_status==0){  ham_dia = new CMATRIX(ndia,ndia);  } 
  else if(ham_dia_mem_status==1){ check_mat_dimensions(ham_dia,ndia,ndia);  }
  else if(ham_dia_mem_status==2){  ham_dia = new CMATRIX(ndia,ndia);    }

  *ham_dia = ham_dia_;
  ham_dia_mem_status = 1; // Allocated internally


}


/************************ NAC *****************************/
void nHamiltonian::init_nac_dia(){
/**
  Allocate memory for the time-derivative NAC matrix in the diabatic basis
*/

  if(nac_dia_mem_status==0){
    nac_dia = new CMATRIX(ndia, ndia);
    nac_dia_mem_status = 1;
  }
  else{ 
    cout<<"WARNING in init_nac_dia: memory is already allocated\n";
  }

}

void nHamiltonian::set_nac_dia_by_ref(CMATRIX& nac_dia_){
/**
  Setup of the nonadiabatic (time-derivative) coupling matrix in the diabatic basis
*/

  check_mat_dimensions(&nac_dia_, ndia, ndia);

  if(nac_dia_mem_status==0){ nac_dia = &nac_dia_; } /// Not allocated - not a problem
  else if(nac_dia_mem_status==1){ delete nac_dia; nac_dia = &nac_dia_;   }
  else if(nac_dia_mem_status==2){ nac_dia = &nac_dia_; }

  nac_dia_mem_status = 2; // Allocated externally


}

void nHamiltonian::set_nac_dia_by_val(CMATRIX& nac_dia_){
/**
  Setup of the nonadiabatic (time-derivative) coupling matrix in the diabatic basis
*/

//  set_X1_by_val(ham_dia, ham_dia_, ham_dia_mem_status, ndia, ndia);

  check_mat_dimensions(&nac_dia_,ndia,ndia);

  if(nac_dia_mem_status==0){  nac_dia = new CMATRIX(ndia,ndia);  } 
  else if(nac_dia_mem_status==1){ check_mat_dimensions(nac_dia,ndia,ndia);  }
  else if(nac_dia_mem_status==2){  nac_dia = new CMATRIX(ndia,ndia);    }

  *nac_dia = nac_dia_;
  nac_dia_mem_status = 1; // Allocated internally


}

/************************ Vibronic Hamiltonian *****************************/
void nHamiltonian::init_hvib_dia(){
/**
  Allocate memory for the vibronic Hamiltonian matrix in the diabatic basis
*/

  if(hvib_dia_mem_status==0){
    hvib_dia = new CMATRIX(ndia, ndia);
    hvib_dia_mem_status = 1;
  }
  else{ 
    cout<<"WARNING in init_hvib_dia: memory is already allocated\n";
  }

}

void nHamiltonian::set_hvib_dia_by_ref(CMATRIX& hvib_dia_){
/**
  Setup of the vibronic Hamiltonian matrix in the diabatic basis
*/

  check_mat_dimensions(&hvib_dia_, ndia, ndia);

  if(hvib_dia_mem_status==0){ hvib_dia = &hvib_dia_; } /// Not allocated - not a problem
  else if(hvib_dia_mem_status==1){ delete hvib_dia; hvib_dia = &hvib_dia_;   }
  else if(hvib_dia_mem_status==2){ hvib_dia = &hvib_dia_; }

  hvib_dia_mem_status = 2; // Allocated externally


}

void nHamiltonian::set_hvib_dia_by_val(CMATRIX& hvib_dia_){
/**
  Setup of the vibronic Hamiltonian matrix in the diabatic basis
*/

//  set_X1_by_val(ham_dia, ham_dia_, ham_dia_mem_status, ndia, ndia);

  check_mat_dimensions(&hvib_dia_,ndia,ndia);

  if(hvib_dia_mem_status==0){  hvib_dia = new CMATRIX(ndia,ndia);  } 
  else if(hvib_dia_mem_status==1){ check_mat_dimensions(hvib_dia,ndia,ndia);  }
  else if(hvib_dia_mem_status==2){  hvib_dia = new CMATRIX(ndia,ndia);    }

  *hvib_dia = hvib_dia_;
  hvib_dia_mem_status = 1; // Allocated internally


}




/************************ 1-ST DERIVATIVES OF HAMILTONIAN *****************************/

void nHamiltonian::init_d1ham_dia(){
/**
  Allocate memory for the 1-st derivatives of the Hamiltonian matrix in the diabatic basis w.r.t. all nuclear DOFs
*/

  for(int n=0;n<nnucl;n++){

    if(d1ham_dia_mem_status[n]==0){
      d1ham_dia[n] = new CMATRIX(ndia, ndia); 
      d1ham_dia_mem_status[n] = 1;
    }
    else{ 
      cout<<"WARNING in init_d1ham_dia: memory for element"<< n <<" is already allocated\n";
    }
  }

}


void nHamiltonian::set_d1ham_dia_by_ref(vector<CMATRIX>& d1ham_dia_){
/**
  Setup of the 1-st derivatives of the Hamiltonian matrix in the diabatic basis w.r.t all nuclear DOFs
*/

//  set_X2_by_ref(d1ham_dia, d1ham_dia_, d1ham_dia_mem_status, ndia, ndia, nnucl);

  check_vector_dimensions(d1ham_dia, d1ham_dia_, d1ham_dia_mem_status, nnucl);

  for(int n=0;n<nnucl;n++){
    check_mat_dimensions(&d1ham_dia_[n], ndia, ndia);

    if(d1ham_dia_mem_status[n]==0){ d1ham_dia[n] = &d1ham_dia_[n]; } 
    else if(d1ham_dia_mem_status[n]==1){ delete d1ham_dia[n]; d1ham_dia[n] = &d1ham_dia_[n];   }
    else if(d1ham_dia_mem_status[n]==2){ d1ham_dia[n] = &d1ham_dia_[n]; }

    d1ham_dia_mem_status[n] = 2; // Allocated externally
  }


}

void nHamiltonian::set_d1ham_dia_by_val(vector<CMATRIX>& d1ham_dia_){
/**
  Setup of the 1-st derivatives of the Hamiltonian matrix in the diabatic basis w.r.t all nuclear DOFs
*/

//  set_X2_by_val(d1ham_dia, d1ham_dia_, d1ham_dia_mem_status, ndia, ndia, nnucl);

  check_vector_dimensions(d1ham_dia, d1ham_dia_, d1ham_dia_mem_status, nnucl);

  for(int n=0;n<nnucl;n++){

    check_mat_dimensions(&d1ham_dia_[n], ndia, ndia);

    if(d1ham_dia_mem_status[n]==0){  d1ham_dia[n] = new CMATRIX(ndia,ndia);  } 
    else if(d1ham_dia_mem_status[n]==1){ check_mat_dimensions(d1ham_dia[n],ndia,ndia);  }
    else if(d1ham_dia_mem_status[n]==2){ d1ham_dia[n] = new CMATRIX(ndia,ndia);    }

    *d1ham_dia[n] = d1ham_dia_[n];
    d1ham_dia_mem_status[n] = 1; // Allocated internally
  }


}

/************************ 2-ND DERIVATIVES OF HAMILTONIAN *****************************/

void nHamiltonian::init_d2ham_dia(){
/**
  Allocate memory for the 2-nd derivatives of the Hamiltonian matrix in the diabatic basis w.r.t. all nuclear DOFs
*/

  for(int n=0;n<nnucl*nnucl;n++){

    if(d2ham_dia_mem_status[n]==0){
      d2ham_dia[n] = new CMATRIX(ndia, ndia); 
      d2ham_dia_mem_status[n] = 1;
    }
    else{ 
      cout<<"WARNING in init_d2ham_dia: memory for element"<< n <<" is already allocated\n";
    }
  }

}


void nHamiltonian::set_d2ham_dia_by_ref(vector<CMATRIX>& d2ham_dia_){
/**
  Setup of the 2-nd derivatives of the Hamiltonian matrix in the diabatic basis w.r.t all nuclear DOFs
*/

//  set_X2_by_ref(d2ham_dia, d2ham_dia_, d2ham_dia_mem_status, ndia, ndia, nnucl);


  check_vector_dimensions(d2ham_dia, d2ham_dia_, d2ham_dia_mem_status, nnucl*nnucl);

  for(int n=0;n<nnucl*nnucl;n++){
    check_mat_dimensions(&d2ham_dia_[n], ndia, ndia);

    if(d2ham_dia_mem_status[n]==0){ d2ham_dia[n] = &d2ham_dia_[n]; } 
    else if(d2ham_dia_mem_status[n]==1){ delete d2ham_dia[n]; d2ham_dia[n] = &d2ham_dia_[n];   }
    else if(d2ham_dia_mem_status[n]==2){ d2ham_dia[n] = &d2ham_dia_[n]; }

    d2ham_dia_mem_status[n] = 2; // Allocated externally
  }


}

void nHamiltonian::set_d2ham_dia_by_val(vector<CMATRIX>& d2ham_dia_){
/**
  Setup of the 2-nd derivatives of the Hamiltonian matrix in the diabatic basis w.r.t all nuclear DOFs
*/

//  set_X2_by_val(d2ham_dia, d2ham_dia_, d2ham_dia_mem_status, ndia, ndia, nnucl);

  check_vector_dimensions(d2ham_dia, d2ham_dia_, d2ham_dia_mem_status, nnucl*nnucl);

  for(int n=0;n<nnucl*nnucl;n++){

    check_mat_dimensions(&d2ham_dia_[n], ndia, ndia);

    if(d2ham_dia_mem_status[n]==0){  d2ham_dia[n] = new CMATRIX(ndia,ndia);  } 
    else if(d2ham_dia_mem_status[n]==1){ check_mat_dimensions(d2ham_dia[n],ndia,ndia);  }
    else if(d2ham_dia_mem_status[n]==2){ d2ham_dia[n] = new CMATRIX(ndia,ndia);    }

    *d2ham_dia[n] = d2ham_dia_[n];
    d2ham_dia_mem_status[n] = 1; // Allocated internally
  }



}


/*************************************************************************

          SETTERS  :         ADIABATIC PARAMETERS

**************************************************************************/



/************************ DERIVATIVE COUPLING *****************************/

void nHamiltonian::init_dc1_adi(){
/**
  Allocate memory for the derivative coupling matrices in the adiabatic basis w.r.t. all nuclear DOFs
*/

  for(int n=0;n<nnucl;n++){

    if(dc1_adi_mem_status[n]==0){
      dc1_adi[n] = new CMATRIX(nadi, nadi); 
      dc1_adi_mem_status[n] = 1;
    }
    else{ 
      cout<<"WARNING in init_dc1_adi: memory for element"<< n <<" is already allocated\n";
    }
  }

}


void nHamiltonian::set_dc1_adi_by_ref(vector<CMATRIX>& dc1_adi_){
/**
  Setup of the derivative coupling matrices in the adiabatic basis w.r.t all nuclear DOFs
*/

//  set_X2_by_ref(dc1_adi, dc1_adi_, dc1_adi_mem_status, nadi, nadi, nnucl);

  check_vector_dimensions(dc1_adi, dc1_adi_, dc1_adi_mem_status, nnucl);

  for(int n=0;n<nnucl;n++){
    check_mat_dimensions(&dc1_adi_[n], nadi, nadi);

    if(dc1_adi_mem_status[n]==0){ dc1_adi[n] = &dc1_adi_[n]; } 
    else if(dc1_adi_mem_status[n]==1){ delete dc1_adi[n]; dc1_adi[n] = &dc1_adi_[n];   }
    else if(dc1_adi_mem_status[n]==2){ dc1_adi[n] = &dc1_adi_[n]; }

    dc1_adi_mem_status[n] = 2; // Allocated externally
  }



}

void nHamiltonian::set_dc1_adi_by_val(vector<CMATRIX>& dc1_adi_){
/**
  Setup of the derivative coupling matrices in the adiabatic basis w.r.t all nuclear DOFs
*/

//  set_X2_by_val(dc1_adi, dc1_adi_, dc1_adi_mem_status, nadi, nadi, nnucl);

  check_vector_dimensions(dc1_adi, dc1_adi_, dc1_adi_mem_status, nnucl);

  for(int n=0;n<nnucl;n++){

    check_mat_dimensions(&dc1_adi_[n], nadi, nadi);

    if(dc1_adi_mem_status[n]==0){  dc1_adi[n] = new CMATRIX(nadi,nadi);  } 
    else if(dc1_adi_mem_status[n]==1){ check_mat_dimensions(dc1_adi[n],nadi,nadi);  }
    else if(dc1_adi_mem_status[n]==2){ dc1_adi[n] = new CMATRIX(nadi,nadi);    }

    *dc1_adi[n] = dc1_adi_[n];
    dc1_adi_mem_status[n] = 1; // Allocated internally
  }


}



/************************ HAMILTONIAN *****************************/

void nHamiltonian::init_ham_adi(){
/**
  Allocate memory for the Hamiltonian matrix in the adiabatic basis
*/

  if(ham_adi_mem_status==0){
    ham_adi = new CMATRIX(nadi, nadi);
    ham_adi_mem_status = 1;
  }
  else{ 
    cout<<"WARNING in init_ham_adi: memory is already allocated\n";
  }

}


void nHamiltonian::set_ham_adi_by_ref(CMATRIX& ham_adi_){
/**
  Setup of the Hamiltonian matrix in the adiabatic basis
*/

//  set_X1_by_ref(ham_adi, ham_adi_, ham_adi_mem_status, nadi, nadi);
  check_mat_dimensions(&ham_adi_,nadi,nadi);

  if(ham_adi_mem_status==0){  ham_adi = &ham_adi_;  } 
  else if(ham_adi_mem_status==1){ delete ham_adi; ham_adi = &ham_adi_;  }
  else if(ham_adi_mem_status==2){ ham_adi = &ham_adi_;    }

  ham_adi_mem_status = 2; // Allocated externally

}

void nHamiltonian::set_ham_adi_by_val(CMATRIX& ham_adi_){
/**
  Setup of the Hamiltonian matrix in the adiabatic basis
*/

//  set_X1_by_val(ham_adi, ham_adi_, ham_adi_mem_status, nadi, nadi);
  check_mat_dimensions(&ham_adi_,nadi,nadi);

  if(ham_adi_mem_status==0){  ham_adi = new CMATRIX(nadi,nadi);  } 
  else if(ham_adi_mem_status==1){ check_mat_dimensions(ham_adi,nadi,nadi);  }
  else if(ham_adi_mem_status==2){  ham_adi = new CMATRIX(nadi,nadi);    }

  *ham_adi = ham_adi_;
  ham_adi_mem_status = 1; // Allocated internally

}


/************************ NAC *****************************/

void nHamiltonian::init_nac_adi(){
/**
  Allocate memory for the nonadiabatic (time-derivative) coupling matrix in the adiabatic basis
*/

  if(nac_adi_mem_status==0){
    nac_adi = new CMATRIX(nadi, nadi);
    nac_adi_mem_status = 1;
  }
  else{ 
    cout<<"WARNING in init_nac_adi: memory is already allocated\n";
  }

}


void nHamiltonian::set_nac_adi_by_ref(CMATRIX& nac_adi_){
/**
  Setup of the nonadiabatic (time-derivative) coupling matrix in the adiabatic basis
*/

  check_mat_dimensions(&nac_adi_, nadi, nadi);

  if(nac_adi_mem_status==0){ nac_adi = &nac_adi_; } /// Not allocated - not a problem
  else if(nac_adi_mem_status==1){ delete nac_adi; nac_adi = &nac_adi_;   }
  else if(nac_adi_mem_status==2){ nac_adi = &nac_adi_; }

  nac_adi_mem_status = 2; // Allocated externally


}

void nHamiltonian::set_nac_adi_by_val(CMATRIX& nac_adi_){
/**
  Setup of the nonadiabatic (time-derivative) coupling matrix in the adiabatic basis
*/

  check_mat_dimensions(&nac_adi_,nadi,nadi);

  if(nac_adi_mem_status==0){  nac_adi = new CMATRIX(nadi,nadi);  } 
  else if(nac_adi_mem_status==1){ check_mat_dimensions(nac_adi,nadi,nadi);  }
  else if(nac_adi_mem_status==2){  nac_adi = new CMATRIX(nadi,nadi);    }

  *nac_adi = nac_adi_;
  nac_adi_mem_status = 1; // Allocated internally


}

/************************ Vibronic Hamiltonian *****************************/

void nHamiltonian::init_hvib_adi(){
/**
  Allocate memory for the vibronic Hamiltonian matrix in the adiabatic basis
*/

  if(hvib_adi_mem_status==0){
    hvib_adi = new CMATRIX(nadi, nadi);
    hvib_adi_mem_status = 1;
  }
  else{ 
    cout<<"WARNING in init_hvib_adi: memory is already allocated\n";
  }

}


void nHamiltonian::set_hvib_adi_by_ref(CMATRIX& hvib_adi_){
/**
  Setup of the vibronic Hamiltonian matrix in the adiabatic basis
*/

  check_mat_dimensions(&hvib_adi_, nadi, nadi);

  if(hvib_adi_mem_status==0){ hvib_adi = &hvib_adi_; } /// Not allocated - not a problem
  else if(hvib_adi_mem_status==1){ delete hvib_adi; hvib_adi = &hvib_adi_;   }
  else if(hvib_adi_mem_status==2){ hvib_adi = &hvib_adi_; }

  hvib_adi_mem_status = 2; // Allocated externally


}

void nHamiltonian::set_hvib_adi_by_val(CMATRIX& hvib_adi_){
/**
  Setup of the vibronic Hamiltonian matrix in the adiabatic basis
*/

  check_mat_dimensions(&hvib_adi_,nadi,nadi);

  if(hvib_adi_mem_status==0){  hvib_adi = new CMATRIX(nadi,nadi);  } 
  else if(hvib_adi_mem_status==1){ check_mat_dimensions(hvib_adi,nadi,nadi);  }
  else if(hvib_adi_mem_status==2){  hvib_adi = new CMATRIX(nadi,nadi);    }

  *hvib_adi = hvib_adi_;
  hvib_adi_mem_status = 1; // Allocated internally


}




void nHamiltonian::init_d1ham_adi(){
/**
  Allocate memory for the 1-st derivatives of the Hamiltonian matrix w.r.t. all nuclear DOFs
*/

  for(int n=0;n<nnucl;n++){

    if(d1ham_adi_mem_status[n]==0){
      d1ham_adi[n] = new CMATRIX(nadi, nadi); 
      d1ham_adi_mem_status[n] = 1;
    }
    else{ 
      cout<<"WARNING in init_d1ham_adi: memory for element"<< n <<" is already allocated\n";
    }
  }

}


void nHamiltonian::set_d1ham_adi_by_ref(vector<CMATRIX>& d1ham_adi_){
/**
  Setup of the 1-st derivatives of the Hamiltonian matrix in the adiabatic basis w.r.t all nuclear DOFs
*/

//  set_X2_by_ref(d1ham_adi, d1ham_adi_, d1ham_adi_mem_status, nadi, nadi, nnucl);

  check_vector_dimensions(d1ham_adi, d1ham_adi_, d1ham_adi_mem_status, nnucl);

  for(int n=0;n<nnucl;n++){
    check_mat_dimensions(&d1ham_adi_[n], nadi, nadi);

    if(d1ham_adi_mem_status[n]==0){ d1ham_adi[n] = &d1ham_adi_[n]; } 
    else if(d1ham_adi_mem_status[n]==1){ delete d1ham_adi[n]; d1ham_adi[n] = &d1ham_adi_[n];   }
    else if(d1ham_adi_mem_status[n]==2){ d1ham_adi[n] = &d1ham_adi_[n]; }

    d1ham_adi_mem_status[n] = 2; // Allocated externally
  }



}


void nHamiltonian::set_d1ham_adi_by_val(vector<CMATRIX>& d1ham_adi_){
/**
  Setup of the 1-st derivatives of the Hamiltonian matrix in the adiabatic basis w.r.t all nuclear DOFs
*/

//  set_X2_by_val(d1ham_adi, d1ham_adi_, d1ham_adi_mem_status, nadi, nadi, nnucl);

  check_vector_dimensions(d1ham_adi, d1ham_adi_, d1ham_adi_mem_status, nnucl);

  for(int n=0;n<nnucl;n++){

    check_mat_dimensions(&d1ham_adi_[n], nadi, nadi);

    if(d1ham_adi_mem_status[n]==0){  d1ham_adi[n] = new CMATRIX(nadi,nadi);  } 
    else if(d1ham_adi_mem_status[n]==1){ check_mat_dimensions(d1ham_adi[n],nadi,nadi);  }
    else if(d1ham_adi_mem_status[n]==2){ d1ham_adi[n] = new CMATRIX(nadi,nadi);    }

    *d1ham_adi[n] = d1ham_adi_[n];
    d1ham_adi_mem_status[n] = 1; // Allocated internally
  }



}


void nHamiltonian::init_d2ham_adi(){
/**
  Allocate memory for the 2-nd derivatives of the Hamiltonian matrix w.r.t. all nuclear DOFs
*/

  for(int n=0;n<nnucl*nnucl;n++){

    if(d2ham_adi_mem_status[n]==0){
      d2ham_adi[n] = new CMATRIX(nadi, nadi); 
      d2ham_adi_mem_status[n] = 1;
    }
    else{ 
      cout<<"WARNING in init_d2ham_adi: memory for element"<< n <<" is already allocated\n";
    }
  }

}


void nHamiltonian::set_d2ham_adi_by_ref(vector<CMATRIX>& d2ham_adi_){
/**
  Setup of the 2-nd derivatives of the Hamiltonian matrix in the adiabatic basis w.r.t all nuclear DOFs
*/

//  set_X2_by_ref(d2ham_adi, d2ham_adi_, d2ham_adi_mem_status, nadi, nadi, nnucl);


  check_vector_dimensions(d2ham_adi, d2ham_adi_, d2ham_adi_mem_status, nnucl*nnucl);

  for(int n=0;n<nnucl*nnucl;n++){
    check_mat_dimensions(&d2ham_adi_[n], nadi, nadi);

    if(d2ham_adi_mem_status[n]==0){ d2ham_adi[n] = &d2ham_adi_[n]; } 
    else if(d2ham_adi_mem_status[n]==1){ delete d2ham_adi[n]; d2ham_adi[n] = &d2ham_adi_[n];   }
    else if(d2ham_adi_mem_status[n]==2){ d2ham_adi[n] = &d2ham_adi_[n]; }

    d2ham_adi_mem_status[n] = 2; // Allocated externally
  }



}

void nHamiltonian::set_d2ham_adi_by_val(vector<CMATRIX>& d2ham_adi_){
/**
  Setup of the 2-nd derivatives of the Hamiltonian matrix in the adiabatic basis w.r.t all nuclear DOFs
*/

//  set_X2_by_val(d2ham_adi, d2ham_adi_, d2ham_adi_mem_status, nadi, nadi, nnucl);

  check_vector_dimensions(d2ham_adi, d2ham_adi_, d2ham_adi_mem_status, nnucl*nnucl);

  for(int n=0;n<nnucl*nnucl;n++){

    check_mat_dimensions(&d2ham_adi_[n], nadi, nadi);

    if(d2ham_adi_mem_status[n]==0){  d2ham_adi[n] = new CMATRIX(nadi,nadi);  } 
    else if(d2ham_adi_mem_status[n]==1){ check_mat_dimensions(d2ham_adi[n],nadi,nadi);  }
    else if(d2ham_adi_mem_status[n]==2){ d2ham_adi[n] = new CMATRIX(nadi,nadi);    }

    *d2ham_adi[n] = d2ham_adi_[n];
    d2ham_adi_mem_status[n] = 1; // Allocated internally
  }



}



void nHamiltonian::init_basis_transform(){
/**
  Allocate memory for the basis transformation matrix
*/

  if(basis_transform_mem_status==0){
    basis_transform = new CMATRIX(ndia, nadi);
    basis_transform_mem_status = 1;
  }
  else{ 
    cout<<"WARNING in init_basis_transform: memory is already allocated\n";
  }

}

void nHamiltonian::set_basis_transform_by_ref(CMATRIX& basis_transform_){
/**
  Setup of the basis transform matrix
*/

//  set_X1_by_ref(basis_transform, basis_transform_, basis_transform_mem_status, ndia, nadi);
  check_mat_dimensions(&basis_transform_,ndia,nadi);

  if(basis_transform_mem_status==0){  basis_transform = &basis_transform_;  } 
  else if(basis_transform_mem_status==1){ delete basis_transform; basis_transform = &basis_transform_;  }
  else if(basis_transform_mem_status==2){ basis_transform = &basis_transform_;    }

  basis_transform_mem_status = 2; // Allocated externally

}

void nHamiltonian::set_basis_transform_by_val(CMATRIX& basis_transform_){
/**
  Setup of the basis transform matrix
*/

//  set_X1_by_val(basis_transform, basis_transform_, basis_transform_mem_status, ndia, nadi);

  check_mat_dimensions(&basis_transform_,ndia,nadi);

  if(basis_transform_mem_status==0){  basis_transform = new CMATRIX(ndia,nadi);  } 
  else if(basis_transform_mem_status==1){ check_mat_dimensions(basis_transform,ndia,nadi);  }
  else if(basis_transform_mem_status==2){  basis_transform = new CMATRIX(ndia,nadi);    }

  *basis_transform = basis_transform_;
  basis_transform_mem_status = 1; // Allocated internally

}


void nHamiltonian::set_ordering_adi_by_ref(vector<int>& ordering_adi_){
/**
  Setup of the variable that takes care of the ordering of the adiabatic states 
*/

  if(ordering_adi_.size()!=nadi){
    cout<<"ERROR in void nHamiltonian::set_ordering_adi_by_ref(vector<int>& ordering_adi_): \
          The external variable must be allocated and its size should be = "<<nadi<<". The \
          current size is = "<<ordering_adi_.size()<<"\nExiting...\n";
    exit(0);
  }
  
  ordering_adi = &ordering_adi_;

}

void nHamiltonian::set_ordering_adi_by_val(vector<int>& ordering_adi_){
/**
  Setup of the variable that takes care of the ordering of the adiabatic states 
*/

  if(ordering_adi_.size()!=nadi){
    cout<<"ERROR in void nHamiltonian::set_ordering_adi_by_ref(vector<int>& ordering_adi_): \
          The external variable must be allocated and its size should be = "<<nadi<<". The \
          current size is = "<<ordering_adi_.size()<<"\nExiting...\n";
    exit(0);
  }
  
  *ordering_adi = ordering_adi_;

}


void nHamiltonian::init_cum_phase_corr(){
/**
  Allocate memory for the 
*/

  if(cum_phase_corr_mem_status==0){
    cum_phase_corr = new CMATRIX(nadi, 1);
    cum_phase_corr_mem_status = 1;
  }
  else{ 
    cout<<"WARNING in void nHamiltonian::init_cum_phase_corr(): memory is already allocated\n";
  }

}

void nHamiltonian::set_cum_phase_corr_by_ref(CMATRIX& cum_phase_corr_){
/**
  Setup of the cumulative phase corrections 
*/

//  set_X1_by_ref(basis_transform, basis_transform_, basis_transform_mem_status, ndia, nadi);
  check_mat_dimensions(&cum_phase_corr_,nadi,1);

  if(cum_phase_corr_mem_status==0){  cum_phase_corr = &cum_phase_corr_;  } 
  else if(cum_phase_corr_mem_status==1){ delete cum_phase_corr; cum_phase_corr = &cum_phase_corr_;  }
  else if(cum_phase_corr_mem_status==2){ cum_phase_corr = &cum_phase_corr_;    }

  cum_phase_corr_mem_status = 2; // Allocated externally

}

void nHamiltonian::set_cum_phase_corr_by_val(CMATRIX& cum_phase_corr_){
/**
  Setup of the basis transform matrix
*/

//  set_X1_by_val(basis_transform, basis_transform_, basis_transform_mem_status, ndia, nadi);
  check_mat_dimensions(&cum_phase_corr_,nadi,1);

  if(cum_phase_corr_mem_status==0){  cum_phase_corr = new CMATRIX(nadi,1);  } 
  else if(cum_phase_corr_mem_status==1){ check_mat_dimensions(cum_phase_corr,nadi,1);  }
  else if(cum_phase_corr_mem_status==2){  cum_phase_corr = new CMATRIX(nadi,1);    }

  *cum_phase_corr = cum_phase_corr_;
  cum_phase_corr_mem_status = 1; // Allocated internally

}



/*************************************************************************

          GETTERS  :         DIABATIC PARAMETERS

**************************************************************************/

CMATRIX nHamiltonian::get_ovlp_dia(){ 
/**
  Return the overlap matrix in the diabatic basis
*/
  if(ovlp_dia_mem_status==0){
    cout<<"Error in get_ovlp_dia: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *ovlp_dia; 
}

CMATRIX nHamiltonian::get_ovlp_dia(vector<int>& id_){ 
/**
  Return the overlap matrix in the diabatic basis
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_ovlp_dia();    }
    else{ cout<<"ERROR in get_ovlp_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_ovlp_dia(next);
  }

}



CMATRIX nHamiltonian::get_dc1_dia(int i){ 
/**
  Return the derivative matrix in the diabatic basis
*/
  if(dc1_dia_mem_status[i]==0){
    cout<<"Error in get_dc1_dia: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *dc1_dia[i]; 
}

CMATRIX nHamiltonian::get_dc1_dia(int i, vector<int>& id_){ 
/**
  Return the derivative matrix in the diabatic basis
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_dc1_dia(i);    }
    else{ cout<<"ERROR in get_dc1_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_dc1_dia(i,next);
  }
}


CMATRIX nHamiltonian::get_ham_dia(){ 
/**
  Return the Hamiltonian matrix in the diabatic basis
*/
  if(ham_dia_mem_status==0){
    cout<<"Error in get_ham_dia: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *ham_dia; 
}

CMATRIX nHamiltonian::get_ham_dia(vector<int>& id_){ 
/**
  Return the Hamiltonian matrix in the diabatic basis
*/

  if(id_.size()==1){
    if(id_[0]==id){   return get_ham_dia();     }
    else{ cout<<"ERROR in get_ham_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_ham_dia(next);
  }

}





CMATRIX nHamiltonian::get_nac_dia(){ 
/**
  Return the nonadiabatic (time-derivative) coupling matrix in the diabatic basis
*/
  if(nac_dia_mem_status==0){
    cout<<"Error in get_nac_dia: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *nac_dia; 
}

CMATRIX nHamiltonian::get_nac_dia(vector<int>& id_){ 
/**
  Return the nonadiabatic (time-derivative) coupling matrix in the diabatic basis
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_nac_dia();    }
    else{ cout<<"ERROR in get_nac_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_nac_dia(next);
  }

}



CMATRIX nHamiltonian::get_hvib_dia(){ 
/**
  Return the vibronic Hamiltonian matrix in the diabatic basis
*/
  if(hvib_dia_mem_status==0){
    cout<<"Error in get_hvib_dia: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *hvib_dia; 
}

CMATRIX nHamiltonian::get_hvib_dia(vector<int>& id_){ 
/**
  Return the vibronic Hamiltonian matrix in the diabatic basis
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_hvib_dia();    }
    else{ cout<<"ERROR in get_hvib_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_hvib_dia(next);
  }

}



CMATRIX nHamiltonian::get_d1ham_dia(int i){ 
/**
  Return the 1-st derivative of the Hamiltonian matrix in the diabatic basis w.r.t. nuclear DOFs
*/
  if(d1ham_dia_mem_status[i]==0){
    cout<<"Error in get_d1ham_dia: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *d1ham_dia[i]; 
}

CMATRIX nHamiltonian::get_d1ham_dia(int i, vector<int>& id_){ 
/**
  Return the 1-st derivative of the Hamiltonian matrix in the diabatic basis w.r.t. nuclear DOFs
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_d1ham_dia(i);    }
    else{ cout<<"ERROR in get_d1ham_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_d1ham_dia(i, next);
  }

}


CMATRIX nHamiltonian::get_d2ham_dia(int i){ 
/**
  Return the 2-nd derivative of the Hamiltonian matrix in the diabatic basis w.r.t. nuclear DOFs
*/
  if(d2ham_dia_mem_status[i]==0){
    cout<<"Error in get_d2ham_dia: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *d2ham_dia[i]; 
}

CMATRIX nHamiltonian::get_d2ham_dia(int i, vector<int>& id_){ 
/**
  Return the 2-nd derivative of the Hamiltonian matrix in the diabatic basis w.r.t. nuclear DOFs
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_d2ham_dia(i);    }
    else{ cout<<"ERROR in get_d2ham_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_d2ham_dia(i, next);
  }

}


CMATRIX nHamiltonian::get_d2ham_dia(int i,int j){ 
/**
  Return the 2-nd derivative of the Hamiltonian matrix in the diabatic basis w.r.t. nuclear DOFs
*/
  if(d2ham_dia_mem_status[i*nnucl+j]==0){
    cout<<"Error in get_d2ham_dia: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *d2ham_dia[i*nnucl+j]; 
}

CMATRIX nHamiltonian::get_d2ham_dia(int i, int j, vector<int>& id_){ 
/**
  Return the 2-nd derivative of the Hamiltonian matrix in the diabatic basis w.r.t. nuclear DOFs
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_d2ham_dia(i,j);    }
    else{ cout<<"ERROR in get_d2ham_dia: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_d2ham_dia(i, j, next);
  }

}


/*************************************************************************

          GETTERS  :         ADIABATIC PARAMETERS

**************************************************************************/

CMATRIX nHamiltonian::get_dc1_adi(int i){ 
/**
  Return the derivative matrix in the adiabatic basis
*/
  if(dc1_adi_mem_status[i]==0){
    cout<<"Error in get_dc1_adi: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *dc1_adi[i]; 
}

CMATRIX nHamiltonian::get_dc1_adi(int i, vector<int>& id_){ 
/**
  Return the derivative matrix in the adiabatic basis
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_dc1_adi(i);    }
    else{ cout<<"ERROR in get_dc1_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_dc1_adi(i,next);
  }
}



CMATRIX nHamiltonian::get_ham_adi(){ 
/**
  Return the Hamiltonian matrix in the adiabatic basis
*/
  if(ham_adi_mem_status==0){
    cout<<"Error in get_ham_adi: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *ham_adi; 
}

CMATRIX nHamiltonian::get_ham_adi(vector<int>& id_){ 
/**
  Return the Hamiltonian matrix in the adiabatic basis
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_ham_adi();    }
    else{ cout<<"ERROR in get_ham_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_ham_adi(next);
  }
}


CMATRIX nHamiltonian::get_nac_adi(){ 
/**
  Return the nonadiabatic (time-derivative) coupling matrix in the adiabatic basis
*/
  if(nac_adi_mem_status==0){
    cout<<"Error in get_nac_adi: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *nac_adi; 
}

CMATRIX nHamiltonian::get_nac_adi(vector<int>& id_){ 
/**
  Return the nonadiabatic (time-derivative) coupling matrix in the adiabatic basis
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_nac_adi();    }
    else{ cout<<"ERROR in get_nac_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_nac_adi(next);
  }
}


CMATRIX nHamiltonian::get_hvib_adi(){ 
/**
  Return the vibronic Hamiltonian matrix in the adiabatic basis
*/
  if(hvib_adi_mem_status==0){
    cout<<"Error in get_hvib_adi: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *hvib_adi; 
}

CMATRIX nHamiltonian::get_hvib_adi(vector<int>& id_){ 
/**
  Return the vibronic Hamiltonian matrix in the adiabatic basis
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_hvib_adi();    }
    else{ cout<<"ERROR in get_hvib_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_hvib_adi(next);
  }
}


CMATRIX nHamiltonian::get_d1ham_adi(int i){ 
/**
  Return the 1-st derivative of the Hamiltonian matrix in the adiabatic basis w.r.t. nuclear DOFs
*/
  if(d1ham_adi_mem_status[i]==0){
    cout<<"Error in get_d1ham_adi: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *d1ham_adi[i]; 
}

CMATRIX nHamiltonian::get_d1ham_adi(int i, vector<int>& id_){ 
/**
  Return the 1-st derivative of the Hamiltonian matrix in the adiabatic basis w.r.t. nuclear DOFs
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_d1ham_adi(i);    }
    else{ cout<<"ERROR in get_d1ham_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_d1ham_adi(i,next);
  }
}


CMATRIX nHamiltonian::get_d2ham_adi(int i){ 
/**
  Return the 2-nd derivative of the Hamiltonian matrix in the adiabatic basis w.r.t. nuclear DOFs
*/
  if(d2ham_adi_mem_status[i]==0){
    cout<<"Error in get_d2ham_adi: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *d2ham_adi[i]; 
}

CMATRIX nHamiltonian::get_d2ham_adi(int i, vector<int>& id_){ 
/**
  Return the 2-nd derivative of the Hamiltonian matrix in the adiabatic basis w.r.t. nuclear DOFs
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_d2ham_adi(i);    }
    else{ cout<<"ERROR in get_d2ham_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_d2ham_adi(i,next);
  }
}


CMATRIX nHamiltonian::get_d2ham_adi(int i,int j){ 
/**
  Return the 2-nd derivative of the Hamiltonian matrix in the adiabatic basis w.r.t. nuclear DOFs
*/
  if(d2ham_adi_mem_status[i*nnucl+j]==0){
    cout<<"Error in get_d2ham_adi: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *d2ham_adi[i*nnucl+j]; 
}

CMATRIX nHamiltonian::get_d2ham_adi(int i, int j, vector<int>& id_){ 
/**
  Return the 2-nd derivative of the Hamiltonian matrix in the adiabatic basis w.r.t. nuclear DOFs
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_d2ham_adi(i,j);    }
    else{ cout<<"ERROR in get_d2ham_adi: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_d2ham_adi(i,j,next);
  }
}




CMATRIX nHamiltonian::get_basis_transform(){ 
/**
  Return the diabatic-to-adiabatic transformation matrix
*/
  if(basis_transform_mem_status==0){
    cout<<"Error in get_basis_transform: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *basis_transform; 
}


CMATRIX nHamiltonian::get_basis_transform(vector<int>& id_){ 
/**
  Return the diabatic-to-adiabatic transformation matrix
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_basis_transform();    }
    else{ cout<<"ERROR in get_basis_transform: No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_basis_transform(next);
  }
}


vector<int> nHamiltonian::get_ordering_adi(){ 
/**
  Return the permutation describing the ordering of the adiabatic states
*/
  return *ordering_adi; 
}

vector<int> nHamiltonian::get_ordering_adi(vector<int>& id_){ 
/**
  Return the permutation describing the ordering of the adiabatic states
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_ordering_adi();    }
    else{ cout<<"ERROR in vector<int> nHamiltonian::get_ordering_adi(vector<int>& id_): \
                No Hamiltonian matching the requested id\nExiting...\n";
          exit(0);   }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_ordering_adi(next);
  }

}




CMATRIX nHamiltonian::get_cum_phase_corr(){ 
/**
  Return the cumulative phase corrections for all eigenvectors
*/
  if(cum_phase_corr_mem_status==0){
    cout<<"ERROR in CMATRIX nHamiltonian::get_cum_phase_corr(): The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *basis_transform; 
}


CMATRIX nHamiltonian::get_cum_phase_corr(vector<int>& id_){ 
/**
  Return the cumulative phase corrections for all eigenvectors
*/
  if(id_.size()==1){
    if(id_[0]==id){   return get_cum_phase_corr();    }
    else{ cout<<"ERROR in CMATRIX nHamiltonian::get_cum_phase_corr(vector<int>& id_):\
                 No Hamiltonian matching the requested id\n"; exit(0); }
  }
  else{
    vector<int> next(id_.begin()+1,id_.end());
    return children[id_[1]]->get_cum_phase_corr(next);
  }
}








}// namespace libhamiltonian_generic
}// namespace libhamiltonian
}// liblibra

