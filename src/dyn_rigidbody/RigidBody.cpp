/*********************************************************************************
* Copyright (C) 2015-2021 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file RigidBody.cpp
  \brief The file implements the most basic methods (constructors, destructors, load/save, initialization, copying, etc)
  of the RigidBody class 
    
*/

#include "RigidBody.h"

// ======================= Internal methods ============================

/// liblibra namespace
namespace liblibra{

using namespace libio;
using namespace liblinalg;


/// librigidbody namespace
namespace librigidbody{

void RigidBody::init_permutations(){
/**
  \brief Internal function to generate the list of permutation matrices

  Permutation matrices are listed according to:
  http://en.wikipedia.org/wiki/File:Symmetric_group_3;_Cayley_table;_matrices.svg 
*/

  // U[0] * diag(I1,I2,I3) * U[0] = diag(I1,I2,I3)
  // det = 1
  U[0].xx = 1.0; U[0].xy = 0.0; U[0].xz = 0.0;
  U[0].yx = 0.0; U[0].yy = 1.0; U[0].yz = 0.0;
  U[0].zx = 0.0; U[0].zy = 0.0; U[0].zz = 1.0;

  // U[1] * diag(I1,I2,I3) * U[1] = diag(I2,I1,I3)
  // det = -1 => change sign
  U[1].xx = 0.0; U[1].xy = 1.0; U[1].xz = 0.0;
  U[1].yx = 1.0; U[1].yy = 0.0; U[1].yz = 0.0;
  U[1].zx = 0.0; U[1].zy = 0.0; U[1].zz =-1.0;

  // U[2] * diag(I1,I2,I3) * U[2] = diag(I1,I3,I2)
  // det = -1 => change sign
  U[2].xx =-1.0; U[2].xy = 0.0; U[2].xz = 0.0;
  U[2].yx = 0.0; U[2].yy = 0.0; U[2].yz = 1.0;
  U[2].zx = 0.0; U[2].zy = 1.0; U[2].zz = 0.0;

  // U[3] * diag(I1,I2,I3) * U[4] = diag(I3,I1,I2)
  // det = 1
  U[3].xx = 0.0; U[3].xy = 0.0; U[3].xz = 1.0;
  U[3].yx = 1.0; U[3].yy = 0.0; U[3].yz = 0.0;
  U[3].zx = 0.0; U[3].zy = 1.0; U[3].zz = 0.0;

  // U[4] * diag(I1,I2,I3) * U[3] = diag(I2,I3,I1)
  // det = 1
  U[4].xx = 0.0; U[4].xy = 1.0; U[4].xz = 0.0;
  U[4].yx = 0.0; U[4].yy = 0.0; U[4].yz = 1.0;
  U[4].zx = 1.0; U[4].zy = 0.0; U[4].zz = 0.0;

  // U[5] * diag(I1,I2,I3) * U[5] = diag(I3,I2,I1)
  // det = -1 => change sign
  U[5].xx = 0.0; U[5].xy = 0.0; U[5].xz = 1.0;
  U[5].yx = 0.0; U[5].yy =-1.0; U[5].yz = 0.0;
  U[5].zx = 1.0; U[5].zy = 0.0; U[5].zz = 0.0;

  if(P1==NULL){
      P1 = new MATRIX(4,4);
      // P1: (q0,q1,q2,q3) -> (-q1,q0,q3,-q2)
      P1->M[0]  = 0.0; P1->M[1]  =-1.0; P1->M[2]  = 0.0; P1->M[3]  = 0.0;
      P1->M[4]  = 1.0; P1->M[5]  = 0.0; P1->M[6]  = 0.0; P1->M[7]  = 0.0;
      P1->M[8]  = 0.0; P1->M[9]  = 0.0; P1->M[10] = 0.0; P1->M[11] = 1.0;
      P1->M[12] = 0.0; P1->M[13] = 0.0; P1->M[14] =-1.0; P1->M[15] = 0.0;
  }

  if(P2==NULL){
      P2 = new MATRIX(4,4);
      // P2: (q0,q1,q2,q3) -> (-q2,-q3,q0,q1)
      P2->M[0]  = 0.0; P2->M[1]  = 0.0; P2->M[2]  =-1.0; P2->M[3]  = 0.0;
      P2->M[4]  = 0.0; P2->M[5]  = 0.0; P2->M[6]  = 0.0; P2->M[7]  =-1.0;
      P2->M[8]  = 1.0; P2->M[9]  = 0.0; P2->M[10] = 0.0; P2->M[11] = 0.0;
      P2->M[12] = 0.0; P2->M[13] = 1.0; P2->M[14] = 0.0; P2->M[15] = 0.0;
  }

  if(P3==NULL){
      P3 = new MATRIX(4,4);
      // P3: (q0,q1,q2,q3) -> (-q3,q2,-q1,q0)
      P3->M[0]  = 0.0; P3->M[1]  = 0.0; P3->M[2]  = 0.0; P3->M[3]  =-1.0;
      P3->M[4]  = 0.0; P3->M[5]  = 0.0; P3->M[6]  = 1.0; P3->M[7]  = 0.0;
      P3->M[8]  = 0.0; P3->M[9]  =-1.0; P3->M[10] = 0.0; P3->M[11] = 0.0;
      P3->M[12] = 1.0; P3->M[13] = 0.0; P3->M[14] = 0.0; P3->M[15] = 0.0;
  }

}

void RigidBody::init_variables(int is_new){
///  is_new = 1 for Cctor and ctor, =0 for assignment operator

//----------- Parameters ---------------
  MACHPREC = 1e-15;
  IntN = 10;
  tol = 1e-15;
  permutindx = 0;
  invpermutindx = 0;
  IEPS = 1e-3;
  BIG = 1e+10;
  MAX_NO = 1e-15;
  minDet = 1e-10;
  NT = 0;
  SERIES_EXPANSION = 10;

  if(is_new){
    cr = NULL;
    ci = NULL;
    Coeffs = NULL;
  
    P1 = NULL; P2 = NULL; P3 = NULL;
  }

  orderflag = 0;
//------------ Variables -----------------
  rb_centers_size = 0;

  // Total mass
  is_rb_mass = 0;
  is_rb_iM = 0;

  // Center of mass
  rb_cm = 0.0;                is_rb_cm = 1;

  // Linear momentum and velocity of the center of mass
  rb_p = 0.0;                 is_rb_p = 1;
  rb_v = 0.0;                 is_rb_v = 1;

  // Force and torque
  rb_force = 0.0;             is_rb_force = 1;
  rb_torque_e = 0.0;          is_rb_torque_e = 1;

  // Inertia moments
  is_rb_I_I = 0;
  is_rb_I_e = 0; 

  // Inverse inertia moments
  is_rb_invI_I = 0;
  is_rb_invI_e = 0;

  // Rotational constants
  is_rb_A = 0;
  is_rb_B = 0;
  is_rb_C = 0;

  //  Attitude matrix and its transpose
  rb_A_I_to_e.identity();   is_rb_A_I_to_e = 1;
  rb_A_I_to_e_T.identity(); is_rb_A_I_to_e_T = 1;

  // Orientation quaternion and its conjugate momentum
  rb_L.init(0.0,0.0,0.0,1.0);           is_rb_L = 1;
  rb_p_r.init(0.0,0.0,0.0,0.0);         is_rb_p_r = 1;

  // Angular momentum and angular velocity (in body frame)
  rb_l_e.init(0.0,0.0,0.0);            is_rb_l_e = 1;
  rb_w_e.init(0.0,0.0,0.0);            is_rb_w_e = 1;

  is_fixed_translation = 0;
  is_fixed_rotation = 0;

  set_orientation_option = 1;          is_set_orientation_option = 0;


}

void RigidBody::copy_content(const RigidBody& rb){
/// Internal auxiliary function that actually implements the copy constructor

  if(rb.rb_centers_size>0){ rb_centers = rb.rb_centers; rb_centers_size = rb.rb_centers_size;}

  if(rb.is_rb_mass) { rb_mass = rb.rb_mass; is_rb_mass = 1; }
  if(rb.is_rb_iM){ rb_iM = rb.rb_iM; is_rb_iM = 1; }
  
  if(rb.is_rb_cm)   { rb_cm = rb.rb_cm;  is_rb_cm = 1; }

  if(rb.is_rb_p)    { rb_p = rb.rb_p; is_rb_p = 1; }
  if(rb.is_rb_v)    { rb_v = rb.rb_v; is_rb_v = 1; }

  if(rb.is_rb_force)   { rb_force = rb.rb_force; is_rb_force = 1; }
  if(rb.is_rb_torque_e){ rb_torque_e = rb.rb_torque_e; is_rb_torque_e = 1; }

  if(rb.is_rb_I_I)  { rb_I_I = rb.rb_I_I;  is_rb_I_I = 1; }
  if(rb.is_rb_I_e)  { rb_I_e = rb.rb_I_e;  is_rb_I_e = 1; }

  if(rb.is_rb_invI_I)  { rb_invI_I = rb.rb_invI_I;  is_rb_invI_I = 1; }
  if(rb.is_rb_invI_e)  { rb_invI_e = rb.rb_invI_e;  is_rb_invI_e = 1; }

  if(rb.is_rb_A)  { rb_A = rb.rb_A;  is_rb_A = 1; }
  if(rb.is_rb_B)  { rb_B = rb.rb_B;  is_rb_B = 1; }
  if(rb.is_rb_C)  { rb_C = rb.rb_C;  is_rb_C = 1; }

  if(rb.is_rb_A_I_to_e)   { rb_A_I_to_e   = rb.rb_A_I_to_e;    is_rb_A_I_to_e   = 1; }
  if(rb.is_rb_A_I_to_e_T) { rb_A_I_to_e_T = rb.rb_A_I_to_e_T;  is_rb_A_I_to_e_T = 1; }

  if(rb.is_rb_L)    { rb_L   = rb.rb_L;   is_rb_L   = 1; }
  if(rb.is_rb_p_r)  { rb_p_r = rb.rb_p_r; is_rb_p_r = 1; }

  if(rb.is_rb_l_e)  { rb_l_e = rb.rb_l_e; is_rb_l_e = 1; }
  if(rb.is_rb_w_e)  { rb_w_e = rb.rb_w_e; is_rb_w_e = 1; }

  is_fixed_translation = rb.is_fixed_translation;
  is_fixed_rotation = rb.is_fixed_rotation;

  if(is_set_orientation_option) {  set_orientation_option = rb.set_orientation_option; }

//-------------- Copy auxiliary variables ---------------------

  top_l = rb.top_l;  
  top_w = rb.top_w;
  top_wm = rb.top_wm; 
  Iint = rb.Iint;   
  invI = rb.invI;
  I1 = rb.I1; I2 = rb.I2; I3 = rb.I3;
  m = rb.m;  
  eps = rb.eps;
  K = rb.K; Kcompl = rb.Kcompl;  
  q = rb.q;   
  A1 = rb.A1; A2 = rb.A2;
  NT = rb.NT;
  SERIES_EXPANSION = rb.SERIES_EXPANSION;
  if(cr==NULL){  cr = new double[NT]; }
  if(ci==NULL){  ci = new double[NT]; }
  for(int i=0;i<NT;i++){ cr[i] = rb.cr[i]; ci[i] = rb.ci[i]; }
  wp = rb.wp;   
  int arr_sz = (SERIES_EXPANSION*(1+SERIES_EXPANSION)/2);

  orderflag = rb.orderflag;

}

void RigidBody::set(object){

}



//============================= Interface methods =================================


RigidBody::RigidBody(){
/** 
  \brief Constructor

  Initializes parameters and permutation matrices
*/

  // Initialize precision parameters
  init_variables(1);
  // Initialize permutation matrix
  init_permutations();  


}

RigidBody::RigidBody(const RigidBody& rb){
/**
  \brief Copy constructor

  Initializes parameters and permutation matrices (default) 
  and then copies the data from the input object

  \param[in] rb Input parameter

*/
  // Initialize precision parameters
  init_variables(1);
  // Initialize permutation matrix
  init_permutations();
  // Copy content of rb object which is defined
  copy_content(rb);
}

RigidBody& RigidBody::operator=(const RigidBody& rb){
/**
  \brief Assignment operator
*/
  // Initialize precision parameters
  init_variables(0);
  // Initialize permutation matrix
  init_permutations();
  // Copy content of rb object which is defined
  copy_content(rb);
  return *this;
}


RigidBody::~RigidBody(){ 
/**
  \brief Destructor
*/
  if(cr!=NULL) { delete [] cr; cr = NULL; }
  if(ci!=NULL) { delete [] ci; ci = NULL; }
  NT = 0;
  if(P1!=NULL) { delete P1; P1 = NULL; }
  if(P2!=NULL) { delete P2; P2 = NULL; }
  if(P3!=NULL) { delete P3; P3 = NULL; }
  if(Coeffs!=NULL){ delete [] Coeffs; Coeffs = NULL; }
  if(rb_centers.size()>0){ rb_centers.clear(); rb_centers_size = 0; }
}


void RigidBody::show_info(){
/**
  \brief Printing internal variables and parameters of RigidBody object
*/

  std::cout<<"Information about object of type RigidBody\n";
  if(rb_centers_size>0){
    for(int i=0;i<rb_centers_size;i++){
      std::cout<<"  rb_centers["<<i<<"] = "<<rb_centers[i]<<std::endl;
    }
  }
  if(is_rb_mass){ std::cout<<"rb_mass = "<<rb_mass<<std::endl; }
  if(is_rb_iM){ std::cout<<"rb_iM = "<<rb_iM<<endl; }
  if(is_rb_cm){ std::cout<<"rb_cm = "<<rb_cm<<std::endl; }
  if(is_rb_p){ std::cout<<"rb_p = "<<rb_p<<std::endl; }
  if(is_rb_v){ std::cout<<"rb_v = "<<rb_v<<std::endl; }
  if(is_rb_force){std::cout<<"rb_force = "<<rb_force<<std::endl; }
  if(is_rb_I_I){ std::cout<<"rb_I_I = "<<rb_I_I<<std::endl; }
  if(is_rb_I_e){ std::cout<<"rb_I_e = "<<rb_I_e<<std::endl; }
  if(is_rb_invI_I){ std::cout<<"rb_invI_I = "<<rb_invI_I<<std::endl; }
  if(is_rb_invI_e){ std::cout<<"rb_invI_e = "<<rb_invI_e<<std::endl; }
  if(is_rb_A){ std::cout<<"rb_A = "<<rb_A<<std::endl; }
  if(is_rb_B){ std::cout<<"rb_B = "<<rb_B<<std::endl; }
  if(is_rb_C){ std::cout<<"rb_C = "<<rb_C<<std::endl; }
  if(is_rb_A_I_to_e){ std::cout<<"rb_A_I_to_e = "<<rb_A_I_to_e<<std::endl; }
  if(is_rb_A_I_to_e_T){ std::cout<<"rb_A_I_to_e_T = "<<rb_A_I_to_e_T<<std::endl; }
  if(is_rb_L){ std::cout<<"rb_L = "<<rb_L<<std::endl; }
  if(is_rb_p_r){ std::cout<<"rb_p_r = "<<rb_p_r<<std::endl; }
  if(is_rb_l_e){ std::cout<<"rb_l_e = "<<rb_l_e<<std::endl; }
  if(is_rb_w_e){ std::cout<<"rb_w_e = "<<rb_w_e<<std::endl; }
  if(is_rb_torque_e){ std::cout<<"rb_torque_e = "<<rb_torque_e<<std::endl; }
  std::cout<<std::endl;

}


void RigidBody::save(boost::property_tree::ptree& pt,std::string path){
/**
  \brief Save the state of the RigidBody object as a property tree

  Each defined data member is added as a node to the property tree. The nodes are added to 
  the level of the tree controlled by the path variable.
 
  \param[in,out] pt The property tree to which the properties of the RigidBody are added
  \param[in] path The parameter controlling the level of the tree to which the RigidBody members will be added.
*/
  int i;

  if(rb_centers_size>0){  liblinalg::save(pt,path+".rb_centers",rb_centers);    }
  if(is_rb_mass){  libio::save(pt,path+".rb_mass",rb_mass);    }
  if(is_rb_iM){  libio::save(pt,path+".rb_iM",rb_iM);    }
  if(is_rb_cm){  liblinalg::save(pt,path+".rb_cm",rb_cm);    }
  if(is_rb_p){  liblinalg::save(pt,path+".rb_p",rb_p);    }
  if(is_rb_v){  liblinalg::save(pt,path+".rb_v",rb_v);    }
  if(is_rb_force){  liblinalg::save(pt,path+".rb_force",rb_force);    }
  if(is_rb_I_I){  liblinalg::save(pt,path+".rb_I_I",rb_I_I);    }
  if(is_rb_I_e){  liblinalg::save(pt,path+".rb_I_e",rb_I_e);    }
  if(is_rb_invI_I){  liblinalg::save(pt,path+".rb_invI_I",rb_invI_I);    }
  if(is_rb_invI_e){  liblinalg::save(pt,path+".rb_invI_e",rb_invI_e);    }
  if(is_rb_A){  libio::save(pt,path+".rb_A",rb_A);    }
  if(is_rb_B){  libio::save(pt,path+".rb_B",rb_B);    }
  if(is_rb_C){  libio::save(pt,path+".rb_C",rb_C);    }
  if(is_rb_A_I_to_e){  liblinalg::save(pt,path+".rb_A_I_to_e",rb_A_I_to_e);    }
  if(is_rb_A_I_to_e_T){  liblinalg::save(pt,path+".rb_A_I_to_e_T",rb_A_I_to_e_T);    }
  if(is_rb_L){  liblinalg::save(pt,path+".rb_L",rb_L);    }
  if(is_rb_p_r){  liblinalg::save(pt,path+".rb_p_r",rb_p_r);    }
  if(is_rb_l_e){  liblinalg::save(pt,path+".rb_l_e",rb_l_e);    }
  if(is_rb_w_e){  liblinalg::save(pt,path+".rb_w_e",rb_w_e);    }
  if(is_rb_torque_e){  liblinalg::save(pt,path+".rb_torque_e",rb_torque_e);    }
  libio::save(pt,path+".is_fixed_translation",is_fixed_translation);
  libio::save(pt,path+".is_fixed_rotation",is_fixed_rotation);

  // Private variables
  libio::save(pt,path+".MACHPREC",MACHPREC);
  libio::save(pt,path+".IntN",IntN);
  libio::save(pt,path+".tol",tol);
  libio::save(pt,path+".IEPS",IEPS);
  libio::save(pt,path+".BIG",BIG);
  libio::save(pt,path+".MAX_NO",MAX_NO);
  libio::save(pt,path+".minDet",minDet);

  for(i=0;i<8;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    libio::save(pt,path+".U["+rt+"]",U[i]);
  }

  libio::save(pt,path+".permutindx",permutindx);
  libio::save(pt,path+".invpermutindx",invpermutindx);
  libio::save(pt,path+".orderflag",orderflag);
  libio::save(pt,path+".SERIES_EXPANSION",SERIES_EXPANSION);
  liblinalg::save(pt,path+".top_l",top_l);
  liblinalg::save(pt,path+".top_w",top_w);
  liblinalg::save(pt,path+".top_wm",top_wm);
  liblinalg::save(pt,path+".Iint",Iint);
  liblinalg::save(pt,path+".invI",invI);
  libio::save(pt,path+".I1",I1);
  libio::save(pt,path+".I2",I2);
  libio::save(pt,path+".I3",I3);
  libio::save(pt,path+".m",m);
  libio::save(pt,path+".eps",eps);
  libio::save(pt,path+".K",K);
  libio::save(pt,path+".Kcompl",Kcompl);
  libio::save(pt,path+".q",q);
  libio::save(pt,path+".A1",A1);
  libio::save(pt,path+".A2",A2);
  libio::save(pt,path+".NT",NT);

  if(cr!=NULL){
    for(i=0;i<NT;i++){
      stringstream ss(stringstream::in | stringstream::out);
      std::string rt; ss<<i; ss>>rt;
      libio::save(pt,path+".cr["+rt+"]",cr[i]);
    }
  }
  if(ci!=NULL){
    for(i=0;i<NT;i++){
      stringstream ss(stringstream::in | stringstream::out);
      std::string rt; ss<<i; ss>>rt;
      libio::save(pt,path+".ci["+rt+"]",ci[i]);
    }
  }
  libio::save(pt,path+".wp",wp);
  if(P1!=NULL){   liblinalg::save(pt,path+".*P1",*P1); }
  if(P2!=NULL){   liblinalg::save(pt,path+".*P2",*P2); }
  if(P3!=NULL){   liblinalg::save(pt,path+".*P3",*P3); }
  if(Coeffs!=NULL){
    for(i=0;i<SERIES_EXPANSION;i++){
      stringstream ss(stringstream::in | stringstream::out);
      std::string rt; ss<<i; ss>>rt;
      libio::save(pt,path+".Coeffs["+rt+"]",Coeffs[i]);
    }
  }

}

void save(boost::property_tree::ptree& pt,std::string path,vector<RigidBody>& vt){
/**
  \brief Save the state of the vector of RigidBody objects as a property tree

  Each RigidBody object is added as a separate branch. 
 
  \param[in,out] pt The property tree to which the list of the RigidBody objects will be added
  \param[in] path The parameter controlling the level of the tree to which the list of RigidBody will be added.
  \param[in] vt The list of RigidBody objects to be printed out into property tree
*/

  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".RigidBody"+rt);
  }
}


void RigidBody::load(boost::property_tree::ptree& pt,std::string path,int& status){
/**
  \brief Load the state of the RigidBody object from a property tree

  Each data member found in the property tree is extracted as the member of the RigidBody object. The
  status of each found data member is set to 1.
 
  \param[in] pt The property tree from which the properties of the RigidBody will be extracted
  \param[in] path The parameter controlling from which level of the tree we try to extract the Barostat object
  \param[out] status Is the global status of the success of the operation. It is 1 is at least one RigidBody member is found at
              given level of the property tree.
*/


  int st;
  status = 0;

  liblinalg::load(pt,path+".rb_centers",rb_centers,st); if(st==1) { status=1; rb_centers_size = rb_centers.size(); }
  libio::load(pt,path+".rb_mass",rb_mass,is_rb_mass); if(is_rb_mass==1) { status=1;}
  libio::load(pt,path+".rb_iM",rb_iM,is_rb_iM); if(is_rb_iM==1) { status=1;}
  liblinalg::load(pt,path+".rb_cm",rb_cm,is_rb_cm); if(is_rb_cm==1) { status=1;}
  liblinalg::load(pt,path+".rb_p",rb_p,is_rb_p); if(is_rb_p==1) { status=1;}
  liblinalg::load(pt,path+".rb_v",rb_v,is_rb_v); if(is_rb_v==1) { status=1;}
  liblinalg::load(pt,path+".rb_force",rb_force,is_rb_force); if(is_rb_force==1) { status=1;}
  liblinalg::load(pt,path+".rb_I_I",rb_I_I,is_rb_I_I); if(is_rb_I_I==1) { status=1;}
  liblinalg::load(pt,path+".rb_I_e",rb_I_e,is_rb_I_e); if(is_rb_I_e==1) { status=1;}
  liblinalg::load(pt,path+".rb_invI_I",rb_invI_I,is_rb_invI_I); if(is_rb_invI_I==1) { status=1;}
  liblinalg::load(pt,path+".rb_invI_e",rb_invI_e,is_rb_invI_e); if(is_rb_invI_e==1) { status=1;}
  libio::load(pt,path+".rb_A",rb_A,is_rb_A); if(is_rb_A==1) { status=1;}
  libio::load(pt,path+".rb_B",rb_B,is_rb_B); if(is_rb_B==1) { status=1;}
  libio::load(pt,path+".rb_C",rb_C,is_rb_C); if(is_rb_C==1) { status=1;}
  liblinalg::load(pt,path+".rb_A_I_to_e",rb_A_I_to_e,is_rb_A_I_to_e); if(is_rb_A_I_to_e==1) { status=1;}
  liblinalg::load(pt,path+".rb_A_I_to_e_T",rb_A_I_to_e_T,is_rb_A_I_to_e_T); if(is_rb_A_I_to_e_T==1) { status=1;}
  liblinalg::load(pt,path+".rb_L",rb_L,is_rb_L); if(is_rb_L==1) { status=1;}
  liblinalg::load(pt,path+".rb_p_r",rb_p_r,is_rb_p_r); if(is_rb_p_r==1) { status=1;}
  liblinalg::load(pt,path+".rb_l_e",rb_l_e,is_rb_l_e); if(is_rb_l_e==1) { status=1;}
  liblinalg::load(pt,path+".rb_w_e",rb_w_e,is_rb_w_e); if(is_rb_w_e==1) { status=1;}
  liblinalg::load(pt,path+".rb_torque_e",rb_torque_e,is_rb_torque_e); if(is_rb_torque_e==1) { status=1;}
  libio::load(pt,path+".is_fixed_translation",is_fixed_translation,st); if(st==1) { status=1;}
  libio::load(pt,path+".is_fixed_rotation",is_fixed_rotation,st); if(st==1) { status=1;}

  libio::load(pt,path+".MACHPREC",MACHPREC,st); if(st==1) { status=1;}
  libio::load(pt,path+".IntN",IntN,st); if(st==1) { status=1;}
  libio::load(pt,path+".tol",tol,st); if(st==1) { status=1;}
  libio::load(pt,path+".IEPS",IEPS,st); if(st==1) { status=1;}
  libio::load(pt,path+".BIG",BIG,st); if(st==1) { status=1;}
  libio::load(pt,path+".MAX_NO",MAX_NO,st); if(st==1) { status=1;}
  libio::load(pt,path+".minDet",minDet,st); if(st==1) { status=1;}
  libio::load(pt,path+".permutindx",permutindx,st); if(st==1) { status=1;}
  libio::load(pt,path+".invpermutindx",invpermutindx,st); if(st==1) { status=1;}
  libio::load(pt,path+".orderflag",orderflag,st); if(st==1) { status=1;}
  libio::load(pt,path+".SERIES_EXPANSION",SERIES_EXPANSION,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".top_l",top_l,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".top_w",top_w,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".top_wm",top_wm,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".Iint",Iint,st); if(st==1) { status=1;}
  liblinalg::load(pt,path+".invI",invI,st); if(st==1) { status=1;}
  libio::load(pt,path+".I1",I1,st); if(st==1) { status=1;}
  libio::load(pt,path+".I2",I2,st); if(st==1) { status=1;}
  libio::load(pt,path+".I3",I3,st); if(st==1) { status=1;}
  libio::load(pt,path+".m",m,st); if(st==1) { status=1;}
  libio::load(pt,path+".eps",eps,st); if(st==1) { status=1;}
  libio::load(pt,path+".K",K,st); if(st==1) { status=1;}
  libio::load(pt,path+".Kcompl",Kcompl,st); if(st==1) { status=1;}
  libio::load(pt,path+".q",q,st); if(st==1) { status=1;}
  libio::load(pt,path+".A1",A1,st); if(st==1) { status=1;}
  libio::load(pt,path+".A2",A2,st); if(st==1) { status=1;}
  libio::load(pt,path+".NT",NT,st); if(st==1) { status=1;}
  libio::load(pt,path+".wp",wp,st); if(st==1) { status=1;}
//  ::load(pt,path+".",,st); if(st==1) { status=1;}




}


void load(boost::property_tree::ptree& pt,std::string path,vector<RigidBody>& vt,int& status){
/**
  \brief Load the vector of RigidBody objects from a property tree

  Each RigidBody object is extracted from a separate branch. 
 
  \param[in] pt The property tree from which the vector of RigidBody objects will be extracted
  \param[in] path The parameter controlling from which level of the property tree we will try to extract the vector of RigidBody objects
  \param[out] vt The vector of created RigidBody objects
  \param[out] status Is the global status of the success of the operation. It is 1 is at least one RigidBody object is extracted
*/

  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      RigidBody x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}



}// namespace librigidbody
}// liblibra
