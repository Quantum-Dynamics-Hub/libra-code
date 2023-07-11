/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Hamiltonian_MM.cpp
  \brief The file implements functions for molecular-mechanical Hamiltonian calculations as well
  as the classes for organizing such computations in an object-oriented way: basic functionality of the class
*/

#include "Hamiltonian_MM.h"

/// liblibra namespace
namespace liblibra{


/// libhamiltonian namespace
namespace libatomistic{

/// libhamiltonian_mm namespace
namespace libhamiltonian_mm{

//========================================================
// Basic functionality of the class Hamiltonian_MM
//========================================================


void Hamiltonian_MM::show_info(){
/**
  Prints out the Hamiltonian_MM object properties
*/

  std::cout<<"Hamiltonian_MM properties:"<<std::endl;
  std::cout<<"is_active = "<<is_active<<std::endl; 
  if(is_int_type){std::cout<<"int_type = "<<int_type<<std::endl;   }
  if(is_functional){std::cout<<"functional = "<<functional<<std::endl;   }
  if(is_respa_type){std::cout<<"respa_type = "<<respa_type<<std::endl; }
  if(Box!=NULL){ std::cout<<"Box = "<<*Box<<endl; }
  std::cout<<"kx = "<<kx<<" ky = "<<ky<<" kz = "<<kz<<std::endl;
//  std::cout<<"kx2 = "<<kx2<<" ky2 = "<<ky2<<" kz2 = "<<kz2<<std::endl;

  if(is_int_type && is_functional){
    if(int_type==0){
      cout<<"id1 = "<<data_bond->id1<<" id2 = "<<data_bond->id2<<"\n";
      cout<<"&r1 = "<<data_bond->r1<<" &r2 = "<<data_bond->r2<<"\n";
      cout<<"&f1 = "<<data_bond->f1<<" &f2 = "<<data_bond->f2<<"\n";
      cout<<"K = "<<data_bond->K<<" D = "<<data_bond->D<<" r0 = "<<data_bond->r0<<" alpha = "<<data_bond->alpha<<"\n";
    }
    else if(int_type==1){
      cout<<"id1 = "<<data_angle->id1<<" id2 = "<<data_angle->id2<<" id3 = "<<data_angle->id3<<"\n";
      cout<<"&r1 = "<<data_angle->r1<<" &r2 = "<<data_angle->r2<<" &r3 = "<<data_angle->r3<<"\n";
      cout<<"&f1 = "<<data_angle->f1<<" &f2 = "<<data_angle->f2<<" &f3 = "<<data_angle->f3<<"\n";
      cout<<"coordination = "<<data_angle->coordination<<" k_theta = "<<data_angle->k_theta
          <<" theta_0 = "<<data_angle->theta_0<<" cos_theta_0 = "<<data_angle->cos_theta_0
          <<" C0 = "<<data_angle->C0<<" C1 = "<<data_angle->C1<<" C2 = "<<data_angle->C2<<"\n";
    }
    else if(int_type==2){
      cout<<"id1 = "<<data_dihedral->id1<<" id2 = "<<data_dihedral->id2<<" id3 = "<<data_dihedral->id3<<" id4 = "<<data_dihedral->id4<<"\n";
      cout<<"&r1 = "<<data_dihedral->r1<<" &r2 = "<<data_dihedral->r2<<" &r3 = "<<data_dihedral->r3<<" &r4 = "<<data_dihedral->r4<<"\n";
      cout<<"&f1 = "<<data_dihedral->f1<<" &f2 = "<<data_dihedral->f2<<" &f3 = "<<data_dihedral->f3<<" &f4 = "<<data_dihedral->f4<<"\n";
      cout<<"Vphi = "<<data_dihedral->Vphi<<" phi0 = "<<data_dihedral->phi0
          <<" Vphi1 = "<<data_dihedral->Vphi1<<" Vphi2 = "<<data_dihedral->Vphi2<<" Vphi3 = "<<data_dihedral->Vphi3
          <<" opt = "<<data_dihedral->opt<<" n = "<<data_dihedral->n<<"\n";
    }
    else if(int_type==4){
      cout<<"id1 = "<<data_vdw->id1<<" id2 = "<<data_vdw->id2<<"\n";
      cout<<"&r1 = "<<data_vdw->r1<<" &r2 = "<<data_vdw->r2<<"\n";
      cout<<"&f1 = "<<data_vdw->f1<<" &f2 = "<<data_vdw->f2<<"\n";
      cout<<"sigma = "<<data_vdw->sigma<<" epsilon = "<<data_vdw->epsilon<<" D = "<<data_vdw->D
          <<" r0 = "<<data_vdw->r0<<" alpha = "<<data_vdw->alpha<<" scale = "<<data_vdw->scale<<"\n";
      if(data_vdw->is_cutoff){    cout<<" R_on = "<<data_vdw->R_on<<" R_off = "<<data_vdw->R_off<<"\n";    }
    }
  }
  std::cout<<std::endl;
}


void Hamiltonian_MM::init_variables(){
  is_active = 1;
  is_int_type = 0;         
  is_functional = 0; 
  respa_type = 0;      is_respa_type = 1;

  Box = NULL;
  kx = ky = kz = 0;
//  kx2 = ky2 = kz2 = 0;
  energy = 0.0;        is_energy = 1;
  hessian = 0.0;       is_hessian = 1;
  stress_at = 0.0;     is_stress_at = 1;
  stress_fr = 0.0;     is_stress_fr = 1;
  stress_ml = 0.0;     is_stress_ml = 1;

  data_bond = NULL;
  data_angle = NULL;
  data_dihedral = NULL;
  data_oop = NULL;
  data_vdw = NULL;
  data_elec = NULL;
  data_gay_berne = NULL;
  data_mb = NULL;

}

void Hamiltonian_MM::copy_content(const Hamiltonian_MM& in){
  if(in.is_int_type) { int_type = in.int_type; is_int_type = 1; }
  if(in.is_functional) { functional = in.functional; is_functional = 1; }
  if(in.is_respa_type) { respa_type = in.respa_type; is_respa_type = 1; }
  if(in.Box!=NULL){ Box = in.Box; }
  kx = in.kx;  ky = in.ky;  kz = in.kz;
//  kx2 = in.kx2;  ky2 = in.ky2;  kz2 = in.kz2;
  if(in.is_energy){ energy = in.energy; is_energy = 1; }
  if(in.is_hessian){ hessian = in.hessian; is_hessian = 1; }
  if(in.is_stress_at){ stress_at = in.stress_at; is_stress_at = 1; }
  if(in.is_stress_fr){ stress_fr = in.stress_fr; is_stress_fr = 1; }
  if(in.is_stress_ml){ stress_ml = in.stress_ml; is_stress_ml = 1; }

  // Copy contents
  if(in.data_bond!=NULL){ data_bond = new bond_interaction; *data_bond = *in.data_bond; }
  if(in.data_angle!=NULL){ data_angle = new angle_interaction; *data_angle = *in.data_angle; }
  if(in.data_dihedral!=NULL){ data_dihedral = new dihedral_interaction; *data_dihedral = *in.data_dihedral; }
  if(in.data_oop!=NULL){ data_oop = new oop_interaction; *data_oop = *in.data_oop; }
  if(in.data_vdw!=NULL){ data_vdw = new vdw_interaction; *data_vdw = *in.data_vdw; }
  if(in.data_elec!=NULL){ data_elec = new elec_interaction; *data_elec = *in.data_elec; }
  if(in.data_gay_berne!=NULL){ data_gay_berne = new gay_berne_interaction; *data_gay_berne = *in.data_gay_berne; }
  if(in.data_mb!=NULL){ data_mb = new mb_interaction; *data_mb = *in.data_mb; }

  is_active = in.is_active;

}

Hamiltonian_MM::Hamiltonian_MM(){
/**
  Constructor
  Initializes internal variables to default values
*/

  // Initialize variables to default values
  init_variables();
}

Hamiltonian_MM::Hamiltonian_MM(const Hamiltonian_MM& in){
/**
  Copy constructor
  Initializes internal variables to default values
  Copies the defined content of the source object
*/

  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(in);
}

Hamiltonian_MM& Hamiltonian_MM::operator=(const Hamiltonian_MM& in){
/**
  Assignment operator
  Initializes internal variables to default values
  Copies the defined content of the source object
*/

  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(in);
  return *this;
}

Hamiltonian_MM::~Hamiltonian_MM(){
/**
  Destructor: frees the memory allocated for the internal variables
*/

//  if(Box!=NULL){ delete Box; Box = NULL; }
  if(data_bond!=NULL) { delete data_bond; data_bond = NULL; }
  if(data_angle!=NULL){ delete data_angle;data_angle = NULL;}
  if(data_dihedral!=NULL){ delete data_dihedral;data_dihedral = NULL;}
  if(data_oop!=NULL){delete data_oop; data_oop = NULL;}
  if(data_vdw!=NULL){delete data_vdw; data_vdw = NULL;}
  if(data_elec!=NULL){delete data_elec; data_elec = NULL;}
  if(data_gay_berne!=NULL){delete data_gay_berne;data_gay_berne = NULL;}
  if(data_mb!=NULL){ delete data_mb; data_mb = NULL; }

}

int operator == (const Hamiltonian_MM& i1, const Hamiltonian_MM& i2){
  int res = 1;
  if(i1.int_type!=i2.int_type){ res = 0; }
  if(i1.functional!=i2.functional){ res = 0; }
  if(i1.respa_type!=i2.respa_type){ res = 0; }
  if(i1.Box!=i2.Box){ res = 0; }

//  if(i1.kx!=i2.kx){ res = 0; }
//  if(i1.ky!=i2.ky){ res = 0; }
//  if(i1.kz!=i2.kz){ res = 0; }

//  if(i1.kx2!=i2.kx2){ res = 0; }
//  if(i1.ky2!=i2.ky2){ res = 0; }
//  if(i1.kz2!=i2.kz2){ res = 0; }


  if(res){
    if(i1.int_type==0){
      if(i1.data_bond==NULL && i2.data_bond==NULL){;;}
      else if(i1.data_bond!=NULL && i2.data_bond!=NULL){
        if(i1.data_bond->id1!=i2.data_bond->id1){ res = 0; }
        if(i1.data_bond->id2!=i2.data_bond->id2){ res = 0; }
        if(i1.data_bond->r1!=i2.data_bond->r1){ res = 0; }
        if(i1.data_bond->r2!=i2.data_bond->r2){ res = 0; }
        if(i1.data_bond->f1!=i2.data_bond->f1){ res = 0; }
        if(i1.data_bond->f2!=i2.data_bond->f2){ res = 0; }
        if(i1.data_bond->K!=i2.data_bond->K){ res = 0;  }
        if(i1.data_bond->D!=i2.data_bond->D){ res = 0;  }
        if(i1.data_bond->r0!=i2.data_bond->r0){ res = 0; }
        if(i1.data_bond->alpha!=i2.data_bond->alpha){ res = 0; }
      }
      else{ res = 0; }
    }// bond

    else if(i1.int_type==1){
      if(i1.data_angle==NULL && i2.data_angle==NULL){;;}
      else if(i1.data_angle!=NULL && i2.data_angle!=NULL){
        if(i1.data_angle->id1!=i2.data_angle->id1){ res = 0; }
        if(i1.data_angle->id2!=i2.data_angle->id2){ res = 0; }
        if(i1.data_angle->id3!=i2.data_angle->id3){ res = 0; }
        if(i1.data_angle->r1!=i2.data_angle->r1){ res = 0; }
        if(i1.data_angle->r2!=i2.data_angle->r2){ res = 0; }
        if(i1.data_angle->r3!=i2.data_angle->r3){ res = 0; }
        if(i1.data_angle->f1!=i2.data_angle->f1){ res = 0; }
        if(i1.data_angle->f2!=i2.data_angle->f2){ res = 0; }
        if(i1.data_angle->f3!=i2.data_angle->f3){ res = 0; }
        if(i1.data_angle->k_theta!=i2.data_angle->k_theta){ res = 0;  }
        if(i1.data_angle->theta_0!=i2.data_angle->theta_0){ res = 0;  }
        if(i1.data_angle->cos_theta_0!=i2.data_angle->cos_theta_0){ res = 0; }
        if(i1.data_angle->C0!=i2.data_angle->C0){ res = 0; }
        if(i1.data_angle->C1!=i2.data_angle->C1){ res = 0; }
        if(i1.data_angle->C2!=i2.data_angle->C2){ res = 0; }
        if(i1.data_angle->coordination!=i2.data_angle->coordination){ res = 0; }
      }
      else{ res = 0; }
    }// angle

    if(i1.int_type==2){
      if(i1.data_dihedral==NULL && i2.data_dihedral==NULL){;;}
      else if(i1.data_dihedral!=NULL && i2.data_dihedral!=NULL){
        if(i1.data_dihedral->id1!=i2.data_dihedral->id1){ res = 0; }
        if(i1.data_dihedral->id2!=i2.data_dihedral->id2){ res = 0; }
        if(i1.data_dihedral->id3!=i2.data_dihedral->id3){ res = 0; }
        if(i1.data_dihedral->id4!=i2.data_dihedral->id4){ res = 0; }
        if(i1.data_dihedral->r1!=i2.data_dihedral->r1){ res = 0; }
        if(i1.data_dihedral->r2!=i2.data_dihedral->r2){ res = 0; }
        if(i1.data_dihedral->r3!=i2.data_dihedral->r3){ res = 0; }
        if(i1.data_dihedral->r4!=i2.data_dihedral->r4){ res = 0; }
        if(i1.data_dihedral->f1!=i2.data_dihedral->f1){ res = 0; }
        if(i1.data_dihedral->f2!=i2.data_dihedral->f2){ res = 0; }
        if(i1.data_dihedral->f3!=i2.data_dihedral->f3){ res = 0; }
        if(i1.data_dihedral->f4!=i2.data_dihedral->f4){ res = 0; }
        if(i1.data_dihedral->Vphi!=i2.data_dihedral->Vphi){ res = 0;  }
        if(i1.data_dihedral->Vphi1!=i2.data_dihedral->Vphi1){ res = 0;  }
        if(i1.data_dihedral->Vphi2!=i2.data_dihedral->Vphi2){ res = 0;  }
        if(i1.data_dihedral->Vphi3!=i2.data_dihedral->Vphi3){ res = 0;  }
        if(i1.data_dihedral->phi0!=i2.data_dihedral->phi0){ res = 0;  }
        if(i1.data_dihedral->n!=i2.data_dihedral->n){ res = 0; }
        if(i1.data_dihedral->opt!=i2.data_dihedral->opt){ res = 0; }
      }
      else{ res = 0; }
    }// dihedral

    if(i1.int_type==3){
      if(i1.data_oop==NULL && i2.data_oop==NULL){;;}
      else if(i1.data_oop!=NULL && i2.data_dihedral!=NULL){
        if(i1.data_oop->id1!=i2.data_oop->id1){ res = 0; }
        if(i1.data_oop->id2!=i2.data_oop->id2){ res = 0; }
        if(i1.data_oop->id3!=i2.data_oop->id3){ res = 0; }
        if(i1.data_oop->id4!=i2.data_oop->id4){ res = 0; }
        if(i1.data_oop->r1!=i2.data_oop->r1){ res = 0; }
        if(i1.data_oop->r2!=i2.data_oop->r2){ res = 0; }
        if(i1.data_oop->r3!=i2.data_oop->r3){ res = 0; }
        if(i1.data_oop->r4!=i2.data_oop->r4){ res = 0; }
        if(i1.data_oop->f1!=i2.data_oop->f1){ res = 0; }
        if(i1.data_oop->f2!=i2.data_oop->f2){ res = 0; }
        if(i1.data_oop->f3!=i2.data_oop->f3){ res = 0; }
        if(i1.data_oop->f4!=i2.data_oop->f4){ res = 0; }
        if(i1.data_oop->K!=i2.data_oop->K){ res = 0;  }
        if(i1.data_oop->C0!=i2.data_oop->C0){ res = 0;  }
        if(i1.data_oop->C1!=i2.data_oop->C1){ res = 0;  }
        if(i1.data_oop->C2!=i2.data_oop->C2){ res = 0;  }
        if(i1.data_oop->xi_0!=i2.data_oop->xi_0){ res = 0;  }
        if(i1.data_oop->opt!=i2.data_oop->opt){ res = 0; }
      }
      else{ res = 0; }
    }// oop

    else if(i1.int_type==4){
      if(i1.data_vdw==NULL && i2.data_vdw==NULL){;;}
      else if(i1.data_vdw!=NULL && i2.data_vdw!=NULL){
        if(i1.data_vdw->id1!=i2.data_vdw->id1){ res = 0; }
        if(i1.data_vdw->id2!=i2.data_vdw->id2){ res = 0; }
        if(i1.data_vdw->r1!=i2.data_vdw->r1){ res = 0; }
        if(i1.data_vdw->r2!=i2.data_vdw->r2){ res = 0; }
        if(i1.data_vdw->f1!=i2.data_vdw->f1){ res = 0; }
        if(i1.data_vdw->f2!=i2.data_vdw->f2){ res = 0; }
        if(i1.data_vdw->sigma!=i2.data_vdw->sigma){ res = 0;  }
        if(i1.data_vdw->epsilon!=i2.data_vdw->epsilon){ res = 0;  }
        if(i1.data_vdw->D!=i2.data_vdw->D){ res = 0;  }
        if(i1.data_vdw->r0!=i2.data_vdw->r0){ res = 0; }
        if(i1.data_vdw->alpha!=i2.data_vdw->alpha){ res = 0; }
        if(i1.data_vdw->scale!=i2.data_vdw->scale){ res = 0; }
      }
      else{ res = 0; }
    }// vdw

    else if(i1.int_type==5){
      if(i1.data_elec==NULL && i2.data_elec==NULL){;;}
      else if(i1.data_elec!=NULL && i2.data_elec!=NULL){
        if(i1.data_elec->id1!=i2.data_elec->id1){ res = 0; }
        if(i1.data_elec->id2!=i2.data_elec->id2){ res = 0; }
        if(i1.data_elec->r1!=i2.data_elec->r1){ res = 0; }
        if(i1.data_elec->r2!=i2.data_elec->r2){ res = 0; }
        if(i1.data_elec->f1!=i2.data_elec->f1){ res = 0; }
        if(i1.data_elec->f2!=i2.data_elec->f2){ res = 0; }
/*        if(i1.data_elec->q1!=i2.data_elec->q1){ res = 0;  }
        if(i1.data_elec->q2!=i2.data_elec->q2){ res = 0;  }
        if(i1.data_elec->xi1!=i2.data_elec->xi1){ res = 0;  }
        if(i1.data_elec->xi2!=i2.data_elec->xi2){ res = 0; }
        if(i1.data_elec->J!=i2.data_elec->J){ res = 0; }
        if(i1.data_elec->eps!=i2.data_elec->eps){ res = 0; }
        if(i1.data_elec->delta!=i2.data_elec->delta){ res = 0; }
*/
      }
      else{ res = 0; }
    }// elec

    else if(i1.int_type==7){
      if(i1.data_gay_berne==NULL && i2.data_gay_berne==NULL){;;}
      else if(i1.data_gay_berne!=NULL && i2.data_gay_berne!=NULL){
        if(i1.data_gay_berne->id1!=i2.data_gay_berne->id1){ res = 0; }
        if(i1.data_gay_berne->id2!=i2.data_gay_berne->id2){ res = 0; }
        if(i1.data_gay_berne->r1!=i2.data_gay_berne->r1){ res = 0; }
        if(i1.data_gay_berne->r2!=i2.data_gay_berne->r2){ res = 0; }
        if(i1.data_gay_berne->u1!=i2.data_gay_berne->u1){ res = 0; }
        if(i1.data_gay_berne->u2!=i2.data_gay_berne->u2){ res = 0; }
        if(i1.data_gay_berne->f1!=i2.data_gay_berne->f1){ res = 0; }
        if(i1.data_gay_berne->f2!=i2.data_gay_berne->f2){ res = 0; }
        if(i1.data_gay_berne->t1!=i2.data_gay_berne->t1){ res = 0; }
        if(i1.data_gay_berne->t2!=i2.data_gay_berne->t2){ res = 0; }
        if(i1.data_gay_berne->di!=i2.data_gay_berne->di){ res = 0;  }
        if(i1.data_gay_berne->dj!=i2.data_gay_berne->dj){ res = 0;  }
        if(i1.data_gay_berne->li!=i2.data_gay_berne->li){ res = 0;  }
        if(i1.data_gay_berne->lj!=i2.data_gay_berne->lj){ res = 0;  }
        if(i1.data_gay_berne->e0!=i2.data_gay_berne->e0){ res = 0;  }
        if(i1.data_gay_berne->rat!=i2.data_gay_berne->rat){ res = 0;  }
        if(i1.data_gay_berne->dw!=i2.data_gay_berne->dw){ res = 0;  }
        if(i1.data_gay_berne->mu!=i2.data_gay_berne->mu){ res = 0;  }
        if(i1.data_gay_berne->nu!=i2.data_gay_berne->nu){ res = 0;  }
      }
      else{ res = 0; }
    }// gay_berne

    else if(i1.int_type==6){
      if(i1.data_mb==NULL && i2.data_mb==NULL){ ;; }
      else if(i1.data_mb!=NULL && i2.data_mb!=NULL){
        if(i1.data_mb->sz!=i2.data_mb->sz){ res = 0; }
        else{
          for(int i=0;i<i1.data_mb->sz;i++){
            if((i1.data_mb->id+i)!=(i2.data_mb->id+i)){ res = 0; }
          }
        }// else
      }
    }// mb

  }// if res


  return res;
}



}// namespace libhamiltonian_mm
}// namespace libatomistic
}// liblibra
