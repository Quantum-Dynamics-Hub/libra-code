/*********************************************************************************
* Copyright (C) 2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian.cpp
  \brief The file implements the generic Hamiltonian class
    
*/

#include "nHamiltonian.h"
#include <stdlib.h>


/// liblibra namespace
namespace liblibra{

/// libhamiltonian namespace 
namespace libhamiltonian{

/// libhamiltonian_generic namespace 
namespace libhamiltonian_generic{



nHamiltonian::nHamiltonian(int ndia_, int nadi_, int nnucl_){ 
/** Default constructor of base (Hamiltonian) class
*/
  int n;

  ndia = ndia_;                   
  nadi = nadi_;
  nnucl = nnucl_;

  ovlp_dia = NULL;             ovlp_dia_mem_status = 0; 

  ampl_dia = NULL;             ampl_dia_mem_status = 0; 

  dc1_dia = vector<CMATRIX*>(nnucl, NULL); 
  dc1_dia_mem_status = vector<int>(nnucl, 0);

  ham_dia = NULL;              ham_dia_mem_status = 0;

  d1ham_dia = vector<CMATRIX*>(nnucl, NULL);
  d1ham_dia_mem_status = vector<int>(nnucl, 0);

  d2ham_dia = vector<CMATRIX*>(nnucl*nnucl, NULL);
  d2ham_dia_mem_status = vector<int>(nnucl*nnucl, 0);

  den_mat_dia = NULL;          den_mat_dia_mem_status = 0;




  dc1_adi = vector<CMATRIX*>(nnucl, NULL);
  dc1_adi_mem_status = vector<int>(nnucl, 0);

  ampl_adi = NULL;             ampl_adi_mem_status = 0; 

  ham_adi = NULL;              ham_adi_mem_status = 0;

  d1ham_adi = vector<CMATRIX*>(nnucl, NULL);
  d1ham_adi_mem_status = vector<int>(nnucl, 0);

  d2ham_adi = vector<CMATRIX*>(nnucl*nnucl, NULL);
  d2ham_adi_mem_status = vector<int>(nnucl*nnucl, 0);

  den_mat_adi = NULL;          den_mat_adi_mem_status = 0;

  basis_transform = NULL;      basis_transform_mem_status = 0;
}

nHamiltonian::~nHamiltonian(){ 
/**
  Deallocate memory only if it was allocated internally
*/

  int n;
 
  if(ovlp_dia_mem_status == 1){ delete ovlp_dia;  ovlp_dia = NULL; ovlp_dia_mem_status = 0;}

  if(ampl_dia_mem_status == 1){ delete ampl_dia; ampl_dia = NULL; ampl_dia_mem_status = 0;}

  for(n;n<dc1_dia.size();n++){
    if(dc1_dia_mem_status[n] == 1){ delete dc1_dia[n];  dc1_dia[n] = NULL; dc1_dia_mem_status[n] = 0;}
  } 
  dc1_dia.clear();
  dc1_dia_mem_status.clear();

  if(ham_dia_mem_status == 1){ delete ham_dia; ham_dia = NULL; ham_dia_mem_status = 0;}

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

  if(den_mat_dia_mem_status == 1){ delete den_mat_dia; den_mat_dia = NULL; den_mat_dia_mem_status = 0;}



  if(ampl_adi_mem_status == 1){ delete ampl_adi; ampl_adi = NULL; ampl_adi_mem_status = 0;}

  for(n;n<dc1_adi.size();n++){
    if(dc1_adi_mem_status[n] == 1){ delete dc1_adi[n];  dc1_adi[n] = NULL; dc1_adi_mem_status[n] = 0;}
  } 
  dc1_adi.clear();
  dc1_adi_mem_status.clear();


  if(ham_adi_mem_status == 1){ delete ham_adi; ham_adi = NULL; ham_adi_mem_status = 0; }

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

  if(den_mat_adi_mem_status == 1){ delete den_mat_adi; den_mat_adi = NULL; den_mat_adi_mem_status = 0;}

  if(basis_transform_mem_status == 1){ delete basis_transform; basis_transform = NULL; basis_transform_mem_status = 0;}


}



/*************************************************************************

          SETTERS  :         DIABATIC PARAMETERS

**************************************************************************/

/************************ AMPLITUDES *****************************/

void nHamiltonian::set_ampl_dia_by_ref(CMATRIX& ampl_dia_){
/**
  Setup of the amplitudes of the diabatic states
*/

//  set_X1_by_ref(ampl_dia, ampl_dia_, ampl_dia_mem_status, ndia, 1);

  check_mat_dimensions(&ampl_dia_, ndia, 1);

  if(ampl_dia_mem_status==0){ ampl_dia = &ampl_dia_; } /// Not allocated - not a problem
  else if(ampl_dia_mem_status==1){ delete ampl_dia; ampl_dia = &ampl_dia_;   }
  else if(ampl_dia_mem_status==2){ ampl_dia = &ampl_dia_; }

  ampl_dia_mem_status = 2; // Allocated externally


}

void nHamiltonian::set_ampl_dia_by_val(CMATRIX& ampl_dia_){
/**
  Setup of the amplitudes of the diabatic states
*/

//  set_X1_by_val(ampl_dia, ampl_dia_, ampl_dia_mem_status, ndia, 1);

  check_mat_dimensions(&ampl_dia_,ndia,1);

  if(ampl_dia_mem_status==0){  ampl_dia = new CMATRIX(ndia,1);  } 
  else if(ampl_dia_mem_status==1){ check_mat_dimensions(ampl_dia,ndia,1);  }
  else if(ampl_dia_mem_status==2){  ampl_dia = new CMATRIX(ndia,1);    }

  *ampl_dia = ampl_dia_;
  ampl_dia_mem_status = 1; // Allocated internally

}


/************************ OVERLAPS *****************************/

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

/************************ 1-ST DERIVATIVES OF HAMILTONIAN *****************************/

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

/************************ DENSITY MATRIX *****************************/

void nHamiltonian::set_den_mat_dia_by_ref(CMATRIX& den_mat_dia_){
/**
  Setup of the density matrix in the diabatic basis
*/

//  set_X1_by_ref(den_mat_dia, den_mat_dia_, den_mat_dia_mem_status, ndia, ndia);

  check_mat_dimensions(&den_mat_dia_, ndia, ndia);

  if(den_mat_dia_mem_status==0){ den_mat_dia = &den_mat_dia_; } /// Not allocated - not a problem
  else if(den_mat_dia_mem_status==1){ delete den_mat_dia; den_mat_dia = &den_mat_dia_;   }
  else if(den_mat_dia_mem_status==2){ den_mat_dia = &den_mat_dia_; }

  den_mat_dia_mem_status = 2; // Allocated externally



}

void nHamiltonian::set_den_mat_dia_by_val(CMATRIX& den_mat_dia_){
/**
  Setup of the density matrix in the diabatic basis
*/

//  set_X1_by_val(den_mat_dia, den_mat_dia_, den_mat_dia_mem_status, ndia, ndia);
  check_mat_dimensions(&den_mat_dia_,ndia,ndia);

  if(den_mat_dia_mem_status==0){  den_mat_dia = new CMATRIX(ndia,ndia);  } 
  else if(den_mat_dia_mem_status==1){ check_mat_dimensions(den_mat_dia,ndia,ndia);  }
  else if(den_mat_dia_mem_status==2){  den_mat_dia = new CMATRIX(ndia,ndia);    }

  *den_mat_dia = den_mat_dia_;
  den_mat_dia_mem_status = 1; // Allocated internally

}



/*************************************************************************

          SETTERS  :         ADIABATIC PARAMETERS

**************************************************************************/

/************************ AMPLITUDES *****************************/

void nHamiltonian::set_ampl_adi_by_ref(CMATRIX& ampl_adi_){
/**
  Setup of the amplitudes of the adiabatic states
*/

//  set_X1_by_ref(ampl_adi, ampl_adi_, ampl_adi_mem_status, nadi, 1);

  check_mat_dimensions(&ampl_adi_, nadi, 1);

  if(ampl_adi_mem_status==0){ ampl_adi = &ampl_adi_; } /// Not allocated - not a problem
  else if(ampl_adi_mem_status==1){ delete ampl_adi; ampl_adi = &ampl_adi_;   }
  else if(ampl_adi_mem_status==2){ ampl_adi = &ampl_adi_; }

  ampl_adi_mem_status = 2; // Allocated externally


}

void nHamiltonian::set_ampl_adi_by_val(CMATRIX& ampl_adi_){
/**
  Setup of the amplitudes of the adiabatic states
*/

//  set_X1_by_val(ampl_adi, ampl_adi_, ampl_adi_mem_status, nadi, 1);

  check_mat_dimensions(&ampl_adi_,nadi,1);

  if(ampl_adi_mem_status==0){  ampl_adi = new CMATRIX(nadi,1);  } 
  else if(ampl_adi_mem_status==1){ check_mat_dimensions(ampl_adi,nadi,1);  }
  else if(ampl_adi_mem_status==2){  ampl_adi = new CMATRIX(nadi,1);    }

  *ampl_adi = ampl_adi_;
  ampl_adi_mem_status = 1; // Allocated internally

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



void nHamiltonian::set_den_mat_adi_by_ref(CMATRIX& den_mat_adi_){
/**
  Setup of the density matrix in the adiabatic basis
*/

  //set_X1_by_ref(den_mat_adi, den_mat_adi_, den_mat_adi_mem_status, nadi, nadi);
  check_mat_dimensions(&den_mat_adi_,nadi,nadi);

  if(den_mat_adi_mem_status==0){  den_mat_adi = &den_mat_adi_;  } 
  else if(den_mat_adi_mem_status==1){ delete den_mat_adi; den_mat_adi = &den_mat_adi_;  }
  else if(den_mat_adi_mem_status==2){ den_mat_adi = &den_mat_adi_;    }

  den_mat_adi_mem_status = 2; // Allocated externally

}

void nHamiltonian::set_den_mat_adi_by_val(CMATRIX& den_mat_adi_){
/**
  Setup of the density matrix in the adiabatic basis
*/

//  set_X1_by_val(den_mat_adi, den_mat_adi_, den_mat_adi_mem_status, nadi, nadi);
  check_mat_dimensions(&den_mat_adi_,nadi,nadi);

  if(den_mat_adi_mem_status==0){  den_mat_adi = new CMATRIX(nadi,nadi);  } 
  else if(den_mat_adi_mem_status==1){ check_mat_dimensions(den_mat_adi,nadi,nadi);  }
  else if(den_mat_adi_mem_status==2){  den_mat_adi = new CMATRIX(nadi,nadi);    }

  *den_mat_adi = den_mat_adi_;
  den_mat_adi_mem_status = 1; // Allocated internally

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



/*************************************************************************

          GETTERS  :         DIABATIC PARAMETERS

**************************************************************************/

CMATRIX nHamiltonian::get_ampl_dia(){ 
/**
  Return the amplitudes of the diabatic basis functions in the overall wavefunction
*/
  if(ampl_dia_mem_status==0){
    cout<<"Error in get_ampl_dia: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *ampl_dia; 
}


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

CMATRIX nHamiltonian::get_den_mat_dia(){ 
/**
  Return the density matrix in the diabatic basis 
*/
  if(den_mat_dia_mem_status==0){
    cout<<"Error in get_den_mat_dia: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *den_mat_dia; 
}


/*************************************************************************

          GETTERS  :         ADIABATIC PARAMETERS

**************************************************************************/

CMATRIX nHamiltonian::get_ampl_adi(){ 
/**
  Return the amplitudes of the adiabatic basis functions in the overall wavefunction
*/
  if(ampl_adi_mem_status==0){
    cout<<"Error in get_ampl_adi: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *ampl_adi; 
}


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

CMATRIX nHamiltonian::get_den_mat_adi(){ 
/**
  Return the density matrix in the adiabatic basis 
*/
  if(den_mat_adi_mem_status==0){
    cout<<"Error in get_den_mat_adi: The matrix is not allocated anywhere\nExiting...\n";
    exit(0);
  }

  return *den_mat_adi; 
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






}// namespace libhamiltonian_generic
}// namespace libhamiltonian
}// liblibra

