#include "Interaction.h"
#include "Cell.h"

void Interaction::set_pbc(MATRIX3x3* box,int kx_,int ky_,int kz_){
  Box = box;
  Box_old = *Box;
  kx = kx_; ky = ky_; kz = kz_;
}

int Interaction::is_origin(){
  return ((kx==0) && (ky==0) && (kz==0));
}

void Interaction::set_interaction_type_and_functional(std::string t,std::string f){

  if(t=="bond"){ int_type = 0; is_int_type = 1;
    if(f=="Harmonic")    { functional = 0; is_functional = 1; }
    else if(f=="Quartic"){ functional = 1; is_functional = 1; }
    else if(f=="Morse")  { functional = 2; is_functional = 1; }
    else{ std::cout<<"Warning: Bond functional "<<f<<" is not implemented\n";  }
  }
  else if(t=="angle"){ int_type = 1; is_int_type = 1; 
    if(f=="Harmonic")                 { functional = 0; is_functional = 1; }
    else if(f=="Fourier")             { functional = 1; is_functional = 1; }
    else if(f=="Fourier_General")     { functional = 2; is_functional = 1; }
    else if(f=="Fourier_Special")     { functional = 3; is_functional = 1; }
    else if(f=="Harmonic_Cos")        { functional = 4; is_functional = 1; }
    else if(f=="Harmonic_Cos_General"){ functional = 5; is_functional = 1; }
    else if(f=="Cubic")               { functional = 6; is_functional = 1; }
    else{ std::cout<<"Warning: Angle functional "<<f<<" is not implemented\n"; }
  }
  else if(t=="dihedral"){ int_type = 2; is_int_type = 1; 
    if((f=="General0")|| (f=="General1") || (f=="General2") || (f=="General3")){ functional = 0; is_functional = 1; }
    else if((f=="Fourier0") || (f=="Fourier1"))                                { functional = 1; is_functional = 1; }
    else{ std::cout<<"Warning: Dihedral/Torsion functional "<<f<<" is not implemented\n";}
  }
  else if(t=="oop"){ int_type = 3; is_int_type = 1; 
    if(f=="Fourier")      { functional = 0; is_functional = 1; }
    else if(f=="Wilson")  { functional = 1; is_functional = 1; }
    else if(f=="Harmonic"){ functional = 2; is_functional = 1; }
    else{ std::cout<<"Warning: Oop functional "<<f<<" is not implemented\n"; }
  }
  else if(t=="vdw"){ int_type = 4; is_int_type = 1; 
    if(f=="LJ")               { functional = 0; is_functional = 1; }
    else if(f=="Buffered14_7"){ functional = 1; is_functional = 1; }
    else if(f=="Morse")       { functional = 2; is_functional = 1; }
    else{ std::cout<<"Warning: Vdw functional "<<f<<" is not implemented\n"; }
  }
  else if(t=="elec"){ int_type = 5; is_int_type = 1; 
    if(f=="Coulomb"){ functional = 0; is_functional = 1; }
    else{ std::cout<<"Warning: Elec functional "<<f<<" is not implemented\n"; }
  }
  else if(t=="mb"){ int_type = 6; is_int_type = 1;
    if(f=="Ewald_3D"){ functional = 0; is_functional = 1; }
    else if(f=="vdw_LJ"){ functional = 1; is_functional = 1; } 
    else if(f=="vdw_LJ1"){ functional = 2; is_functional = 1; }
    else{ std::cout<<"Warning: Many-body potential "<<f<<" is not implemented\n"; }
  }
  else if(t=="cg"){ int_type = 7; is_int_type = 1; 
    if(f=="Gay-Berne"){ functional = 0; is_functional = 1; }
    else{ std::cout<<"Warning: Coarse-grain (fragment-fragment) potential "<<f<<" is not implemented\n"; }
  }
  else if(t=="mb_excl"){ int_type = 8; is_int_type = 1;
    if(f=="vdw_LJ1"){ functional = 0; is_functional = 1; }
    else{ std::cout<<"Warning: Many-body exclusion potential "<<f<<" is not implemented\n"; }
  }

  else{
    std::cout<<"Warning: Interaction type "<<t<<" is not known\n";
  }


}

void Interaction::set_2a_interaction(std::string t,std::string f,
                                     int id1,int id2,
                                     VECTOR& r1,VECTOR& r2,
                                     VECTOR& g1,VECTOR& g2,
                                     VECTOR& m1,VECTOR& m2,
                                     VECTOR& f1,VECTOR& f2, 
                                     double& dr1_2, double& dr2_2,  double& dT_2,
                                     map<std::string,double> params){

  set_interaction_type_and_functional(t,f);

  if(int_type==0){ // bond
    if(data_bond==NULL){  data_bond = new bond_interaction; }
    // Set up parameters
    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="K"){ data_bond->K = it->second; }
      else if(it->first=="D"){ data_bond->D = it->second; }
      else if(it->first=="r0"){ data_bond->r0 = it->second; }
      else if(it->first=="alpha"){ data_bond->alpha = it->second; }
    }
    data_bond->id1 = id1;  data_bond->id2 = id2;
    data_bond->r1 = &r1;   data_bond->r2 = &r2;
    data_bond->g1 = &g1;   data_bond->g2 = &g2;
    data_bond->m1 = &m1;   data_bond->m2 = &m2;
    data_bond->f1 = &f1;   data_bond->f2 = &f2; 
  }// bond
  else if(int_type==4){ // vdw
    if(data_vdw==NULL){  data_vdw = new vdw_interaction; }
    // Set up parameters
    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="sigma"){ data_vdw->sigma = it->second; }
      else if(it->first=="epsilon"){ data_vdw->epsilon = it->second; }
      else if(it->first=="D"){ data_vdw->D = it->second; }
      else if(it->first=="r0"){ data_vdw->r0 = it->second; }
      else if(it->first=="alpha"){ data_vdw->alpha = it->second; }
      else if(it->first=="scale"){ data_vdw->scale = it->second; }
      else if(it->first=="R_on"){ data_vdw->R_on = it->second; }
      else if(it->first=="R_off"){ data_vdw->R_off = it->second; }
      else if(it->first=="R_on2"){ data_vdw->R_on2 = it->second; }
      else if(it->first=="R_off2"){ data_vdw->R_off2 = it->second; }
      else if(it->first=="is_cutoff"){ data_vdw->is_cutoff = it->second; }
      else if(it->first=="time"){ data_vdw->time = it->second; }

    }
    data_vdw->id1 = id1;  data_vdw->id2 = id2;
    data_vdw->r1 = &r1;   data_vdw->r2 = &r2;
    data_vdw->g1 = &g1;   data_vdw->g2 = &g2;
    data_vdw->m1 = &m1;   data_vdw->m2 = &m2;
    data_vdw->f1 = &f1;   data_vdw->f2 = &f2;
    data_vdw->r1_old = r1; data_vdw->r2_old = r2;
    data_vdw->displr1_2 = &dr1_2; data_vdw->displr2_2 = &dr2_2;
    data_vdw->displT_2 = &dT_2;

  }// vdw
  else if(int_type==5){ // elec
    if(data_elec==NULL){  data_elec = new elec_interaction; }
    // Set up parameters
    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="q1"){ data_elec->q1 = it->second; }
      else if(it->first=="q2"){ data_elec->q2 = it->second; }
      else if(it->first=="xi1"){ data_elec->xi1 = it->second; }
      else if(it->first=="xi2"){ data_elec->xi2 = it->second; }
      else if(it->first=="eps"){ data_elec->eps = it->second; }
      else if(it->first=="J"){ data_elec->J = it->second; }
      else if(it->first=="delta"){ data_elec->delta = it->second; }
      else if(it->first=="scale"){ data_elec->scale = it->second; }
      else if(it->first=="R_on"){ data_elec->R_on = it->second; }
      else if(it->first=="R_off"){ data_elec->R_off = it->second; }
      else if(it->first=="R_on2"){ data_elec->R_on2 = it->second; }
      else if(it->first=="R_off2"){ data_elec->R_off2 = it->second; }
      else if(it->first=="is_cutoff"){ data_elec->is_cutoff = it->second; }

    }
    data_elec->id1 = id1;  data_elec->id2 = id2;
    data_elec->r1 = &r1;   data_elec->r2 = &r2;
    data_elec->g1 = &g1;   data_elec->g2 = &g2;
    data_elec->m1 = &m1;   data_elec->m2 = &m2;
    data_elec->f1 = &f1;   data_elec->f2 = &f2;
  }// elec


}
 
void Interaction::set_3a_interaction(std::string t,std::string f,
                                     int id1,int id2,int id3,
                                     VECTOR& r1,VECTOR& r2,VECTOR& r3,
                                     VECTOR& g1,VECTOR& g2,VECTOR& g3,
                                     VECTOR& m1,VECTOR& m2,VECTOR& m3,
                                     VECTOR& f1,VECTOR& f2,VECTOR& f3,
                                      map<std::string,double> params){
  set_interaction_type_and_functional(t,f);

  if(int_type==1){ // angle
    if(data_angle==NULL){  data_angle = new angle_interaction; }
    // Set up parameters
    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="k_theta"){ data_angle->k_theta = it->second; }
      else if(it->first=="theta_0"){ data_angle->theta_0 = it->second; }
      else if(it->first=="cos_theta_0"){ data_angle->cos_theta_0 = it->second; }
      else if(it->first=="C0"){ data_angle->C0 = it->second; }
      else if(it->first=="C1"){ data_angle->C1 = it->second; }
      else if(it->first=="C2"){ data_angle->C2 = it->second; }
      else if(it->first=="coordination"){ data_angle->coordination = (int)it->second; }
    }
    data_angle->id1 = id1;  data_angle->id2 = id2; data_angle->id3 = id3;
    data_angle->r1 = &r1;   data_angle->r2 = &r2;  data_angle->r3 = &r3;
    data_angle->g1 = &g1;   data_angle->g2 = &g2;  data_angle->g3 = &g3;
    data_angle->m1 = &m1;   data_angle->m2 = &m2;  data_angle->m3 = &m3;
    data_angle->f1 = &f1;   data_angle->f2 = &f2;  data_angle->f3 = &f3;
  }// angle


}

void Interaction::set_4a_interaction(std::string t,std::string f,
                                     int id1,int id2,int id3,int id4,
                                     VECTOR& r1,VECTOR& r2,VECTOR& r3,VECTOR& r4,
                                     VECTOR& g1,VECTOR& g2,VECTOR& g3,VECTOR& g4,
                                     VECTOR& m1,VECTOR& m2,VECTOR& m3,VECTOR& m4,
                                     VECTOR& f1,VECTOR& f2,VECTOR& f3,VECTOR& f4,
                                     map<std::string,double> params){

  set_interaction_type_and_functional(t,f);

  if(int_type==2){
    if(data_dihedral==NULL){  data_dihedral = new dihedral_interaction; }
    // Set up parameters
    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="Vphi"){ data_dihedral->Vphi = it->second; }
      else if(it->first=="phi0"){ data_dihedral->phi0 = it->second; }
      else if(it->first=="Vphi1"){ data_dihedral->Vphi1 = it->second; }
      else if(it->first=="Vphi2"){ data_dihedral->Vphi2 = it->second; }
      else if(it->first=="Vphi3"){ data_dihedral->Vphi3 = it->second; }
      else if(it->first=="opt"){ data_dihedral->opt = (int)it->second; }
      else if(it->first=="n"){ data_dihedral->n = (int)it->second; }
    }
    if(f=="General0")     { data_dihedral->opt = 0; }
    else if(f=="General1"){ data_dihedral->opt = 1; }
    else if(f=="General2"){ data_dihedral->opt = 2; }
    else if(f=="General3"){ data_dihedral->opt = 3; }
    else if(f=="Fourier0"){ data_dihedral->opt = 0; }
    else if(f=="Fourier1"){ data_dihedral->opt = 1; }

    data_dihedral->id1 = id1;  data_dihedral->id2 = id2; data_dihedral->id3 = id3; data_dihedral->id4 = id4;
    data_dihedral->r1 = &r1;   data_dihedral->r2 = &r2;  data_dihedral->r3 = &r3;  data_dihedral->r4 = &r4;
    data_dihedral->g1 = &g1;   data_dihedral->g2 = &g2;  data_dihedral->g3 = &g3;  data_dihedral->g4 = &g4;
    data_dihedral->m1 = &m1;   data_dihedral->m2 = &m2;  data_dihedral->m3 = &m3;  data_dihedral->m4 = &m4;
    data_dihedral->f1 = &f1;   data_dihedral->f2 = &f2;  data_dihedral->f3 = &f3;  data_dihedral->f4 = &f4;
  }// dihedral

  else if(int_type==3){
    if(data_oop==NULL){  data_oop = new oop_interaction; }
    // Set up parameters
    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="K"){ data_oop->K = it->second; }
      else if(it->first=="xi_0"){ data_oop->xi_0 = it->second; }
      else if(it->first=="C0"){ data_oop->C0 = it->second; }
      else if(it->first=="C1"){ data_oop->C1 = it->second; }
      else if(it->first=="C2"){ data_oop->C2 = it->second; }
      else if(it->first=="opt"){ data_oop->opt = (int)it->second; }
    }
    data_oop->id1 = id1;  data_oop->id2 = id2; data_oop->id3 = id3; data_oop->id4 = id4;
    data_oop->r1 = &r1;   data_oop->r2 = &r2;  data_oop->r3 = &r3;  data_oop->r4 = &r4;
    data_oop->g1 = &g1;   data_oop->g2 = &g2;  data_oop->g3 = &g3;  data_oop->g4 = &g4;
    data_oop->m1 = &m1;   data_oop->m2 = &m2;  data_oop->m3 = &m3;  data_oop->m4 = &m4;
    data_oop->f1 = &f1;   data_oop->f2 = &f2;  data_oop->f3 = &f3;  data_oop->f4 = &f4;
  }// oop

}

void Interaction::set_mb_interaction(std::string t,std::string p,int sz, int* id, VECTOR** r, VECTOR** g, VECTOR** m, VECTOR** f,double** q,
                                     double** epsilon, double** sigma,
                                     int nexcl, int* excl1, int* excl2, double* scale,double** displr_2, double* displT_2,
                                     vector< vector<excl_scale> >& excl_scales,
                                     map<std::string,double> params){

  set_interaction_type_and_functional(t,p);

  if(data_mb==NULL){  data_mb = new mb_interaction; }
  // Set up parameters
  data_mb->sz = sz;
  data_mb->nexcl = nexcl;
  data_mb->excl1 = excl1;
  data_mb->excl2 = excl2;
  data_mb->scale = scale;
  data_mb->id = id;
  data_mb->r  = r;
  data_mb->g  = g;
  data_mb->m  = m;
  data_mb->f  = f;
  data_mb->q  = q;
  data_mb->epsilon = epsilon;
  data_mb->sigma = sigma;
  data_mb->displr_2 = displr_2; 
  data_mb->displT_2 = displT_2;
  data_mb->excl_scales = excl_scales;

  // Set up general parameters
  for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
    if(it->first=="R_on"){ data_mb->R_on = it->second; }
    else if(it->first=="R_off"){ data_mb->R_off = it->second; }
    else if(it->first=="R_on2"){ data_mb->R_on2 = it->second; }
    else if(it->first=="R_off2"){ data_mb->R_off2 = it->second; }
    else if(it->first=="is_cutoff"){ data_mb->is_cutoff = it->second; }
    else if(it->first=="elec_etha"){ data_mb->elec_etha = it->second; }
    else if(it->first=="time"){ data_mb->time = it->second; }
  }
  data_mb->time = 0;

}

void Interaction::set_2f_interaction(std::string t,std::string f,
                                     int id1,int id2,
                                     VECTOR& r1,VECTOR& r2,VECTOR& u1,VECTOR& u2,
                                     VECTOR& f1,VECTOR& f2,VECTOR& t1,VECTOR& t2,
                                     map<std::string,double> params){

  set_interaction_type_and_functional(t,f);

  if(int_type==7){
    if(data_gay_berne==NULL){  data_gay_berne = new gay_berne_interaction; }
    // Set up parameters
    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="di"){ data_gay_berne->di = it->second; }
      else if(it->first=="dj"){ data_gay_berne->dj = it->second; }
      else if(it->first=="li"){ data_gay_berne->li = it->second; }
      else if(it->first=="lj"){ data_gay_berne->lj = it->second; }
      else if(it->first=="e0"){ data_gay_berne->e0 = it->second; }
      else if(it->first=="rat"){ data_gay_berne->rat = it->second; }
      else if(it->first=="dw"){ data_gay_berne->dw = it->second; }
      else if(it->first=="mu"){ data_gay_berne->mu = it->second; }
      else if(it->first=="nu"){ data_gay_berne->nu = it->second; }
    }
    data_gay_berne->id1 = id1;  data_gay_berne->id2 = id2; 
    data_gay_berne->r1 = &r1;   data_gay_berne->r2 = &r2;  data_gay_berne->u1 = &u1;  data_gay_berne->u2 = &u2;
    data_gay_berne->f1 = &f1;   data_gay_berne->f2 = &f2;  data_gay_berne->t1 = &t1;  data_gay_berne->t2 = &t2;
  }// gay_berne


}

double Interaction::calculate(int& update_displ2){
  return calculate(int_type,update_displ2);
}

double Interaction::calculate(int call_type,int& update_displ2){
  energy = 0.0;
  stress_at = 0.0;
  stress_fr = 0.0;
  stress_ml = 0.0;
  MATRIX3x3 tp;
  update_displ2 = 0;

  if(call_type==int_type){

  if(int_type==0){
    if(is_active){
    VECTOR& r1 = *(data_bond->r1); VECTOR& r2 = *(data_bond->r2);
    VECTOR f1,f2; f1 = f2 = 0.0;
    if(functional==0){ energy = Bond_Harmonic(r1,r2,f1,f2,data_bond->K,data_bond->r0);}
    else if(functional==1){ energy = Bond_Quartic(r1,r2,f1,f2,data_bond->K,data_bond->r0); }
    else if(functional==2){ energy = Bond_Morse(r1,r2,f1,f2,data_bond->D,data_bond->r0,data_bond->alpha); }

    *(data_bond->f1) += f1;
    *(data_bond->f2) += f2;

    tp.tensor_product((*(data_bond->r1) - *(data_bond->r2)) , f1);   stress_at += tp;
    tp.tensor_product((*(data_bond->g1) - *(data_bond->g2)) , f1);   stress_fr += tp;
    tp.tensor_product((*(data_bond->m1) - *(data_bond->m2)) , f1);   stress_ml += tp;

    }
  }// bond
 
  else if(int_type==1){
    if(is_active){

    VECTOR& r1 = *(data_angle->r1); VECTOR& r2 = *(data_angle->r2); VECTOR& r3 = *(data_angle->r3);
    VECTOR f1,f2,f3; f1 = f2 = f3 = 0.0;
//    cout<<"~~~~~~~~~~~~angle interaction "<<data_angle->id1<<" "<<data_angle->id2<<" "<<data_angle->id3<<"~~~~~~~~~~~~~"<<endl;
//    cout<<"&r1 = "<<(data_angle->r1)<<" &r2 = "<<(data_angle->r2)<<" &r3 = "<<(data_angle->r3)<<endl;
//    cout<<" r1 = "<<r1<<" r2 = "<<r2<<" r3 = "<<r3<<endl;
    if(functional==0){ energy = Angle_Harmonic(r1,r2,r3,f1,f2,f3,data_angle->k_theta,data_angle->theta_0); }
    else if(functional==1){ energy = Angle_Fourier(r1,r2,r3,f1,f2,f3,data_angle->k_theta,data_angle->C0, data_angle->C1,data_angle->C2,data_angle->coordination); }
    else if(functional==2){ energy = Angle_Fourier_General(r1,r2,r3,f1,f2,f3,data_angle->k_theta,data_angle->C0, data_angle->C1,data_angle->C2); }
    else if(functional==3){ energy = Angle_Fourier_Special(r1,r2,r3,f1,f2,f3,data_angle->k_theta,data_angle->coordination); }
    else if(functional==4){ energy = Angle_Harmonic_Cos(r1,r2,r3,f1,f2,f3,data_angle->k_theta,data_angle->cos_theta_0,data_angle->coordination); }
    else if(functional==5){ energy = Angle_Harmonic_Cos_General(r1,r2,r3,f1,f2,f3,data_angle->k_theta,data_angle->cos_theta_0); }
    else if(functional==6){ energy = Angle_Cubic(r1,r2,r3,f1,f2,f3,data_angle->k_theta,data_angle->theta_0); }

    *(data_angle->f1) += f1;
    *(data_angle->f2) += f2;
    *(data_angle->f3) += f3;

    tp.tensor_product(*(data_angle->r1) , f1);     stress_at += tp;
    tp.tensor_product(*(data_angle->r2) , f2);     stress_at += tp;
    tp.tensor_product(*(data_angle->r3) , f3);     stress_at += tp;

    tp.tensor_product(*(data_angle->g1) , f1);     stress_fr += tp;
    tp.tensor_product(*(data_angle->g2) , f2);     stress_fr += tp;
    tp.tensor_product(*(data_angle->g3) , f3);     stress_fr += tp;

//    tp.tensor_product(*(data_angle->m1) , f1);     stress_ml += tp;
//    tp.tensor_product(*(data_angle->m2) , f2);     stress_ml += tp;
//    tp.tensor_product(*(data_angle->m3) , f3);     stress_ml += tp;

    }// if is_active
  }// angle

  else if(int_type==2){
    if(is_active){
    VECTOR& r1 = *(data_dihedral->r1); VECTOR& r2 = *(data_dihedral->r2); VECTOR& r3 = *(data_dihedral->r3); VECTOR& r4 = *(data_dihedral->r4);
    VECTOR f1,f2,f3,f4; f1 = f2 = f3 = f4 = 0.0;
    if(functional==0){ energy = Dihedral_General(r1,r2,r3,r4,f1,f2,f3,f4,data_dihedral->Vphi,data_dihedral->phi0,data_dihedral->n,data_dihedral->opt); }
    else if(functional==1){ energy = Dihedral_Fourier(r1,r2,r3,r4,f1,f2,f3,f4,data_dihedral->Vphi1,data_dihedral->Vphi2,data_dihedral->Vphi3,data_dihedral->opt); }

    *(data_dihedral->f1) += f1;
    *(data_dihedral->f2) += f2;
    *(data_dihedral->f3) += f3;
    *(data_dihedral->f4) += f4;

    tp.tensor_product(*(data_dihedral->r1) , f1);     stress_at += tp;
    tp.tensor_product(*(data_dihedral->r2) , f2);     stress_at += tp;
    tp.tensor_product(*(data_dihedral->r3) , f3);     stress_at += tp;
    tp.tensor_product(*(data_dihedral->r4) , f4);     stress_at += tp;

    tp.tensor_product(*(data_dihedral->g1) , f1);     stress_fr += tp;
    tp.tensor_product(*(data_dihedral->g2) , f2);     stress_fr += tp;
    tp.tensor_product(*(data_dihedral->g3) , f3);     stress_fr += tp;
    tp.tensor_product(*(data_dihedral->g4) , f4);     stress_fr += tp;

//    tp.tensor_product(*(data_dihedral->m1) , f1);     stress_ml += tp;
//    tp.tensor_product(*(data_dihedral->m2) , f2);     stress_ml += tp;
//    tp.tensor_product(*(data_dihedral->m3) , f3);     stress_ml += tp;
//    tp.tensor_product(*(data_dihedral->m4) , f4);     stress_ml += tp;

    }
  }// dihedral

  else if(int_type==3){
    if(is_active){
    VECTOR& r1 = *(data_oop->r1); VECTOR& r2 = *(data_oop->r2); VECTOR& r3 = *(data_oop->r3); VECTOR& r4 = *(data_oop->r4);
    VECTOR f1,f2,f3,f4; f1 = f2 = f3 = f4 = 0.0;
    if(functional==0){ energy = OOP_Fourier(r1,r2,r3,r4,f1,f2,f3,f4,data_oop->K,data_oop->C0,data_oop->C1,data_oop->C2,data_oop->opt); }
    else if(functional==1){ energy = OOP_Wilson(r1,r2,r3,r4,f1,f2,f3,f4,data_oop->K,data_oop->xi_0); }
    else if(functional==2){ energy = OOP_Harmonic(r1,r2,r3,r4,f1,f2,f3,f4,data_oop->K); }

    *(data_oop->f1) += f1;
    *(data_oop->f2) += f2;
    *(data_oop->f3) += f3;
    *(data_oop->f4) += f4;

    tp.tensor_product(*(data_oop->r1) , f1);     stress_at += tp;
    tp.tensor_product(*(data_oop->r2) , f2);     stress_at += tp;
    tp.tensor_product(*(data_oop->r3) , f3);     stress_at += tp;
    tp.tensor_product(*(data_oop->r4) , f4);     stress_at += tp;

    tp.tensor_product(*(data_oop->g1) , f1);     stress_fr += tp;
    tp.tensor_product(*(data_oop->g2) , f2);     stress_fr += tp;
    tp.tensor_product(*(data_oop->g3) , f3);     stress_fr += tp;
    tp.tensor_product(*(data_oop->g4) , f4);     stress_fr += tp;

//    tp.tensor_product(*(data_oop->m1) , f1);     stress_ml += tp;
//    tp.tensor_product(*(data_oop->m2) , f2);     stress_ml += tp;
//    tp.tensor_product(*(data_oop->m3) , f3);     stress_ml += tp;
//    tp.tensor_product(*(data_oop->m4) , f4);     stress_ml += tp;

    }
  }// oop
/*
  else if(int_type=="vdw"){
    VECTOR tv1,tv2,tv3; tv1 = tv2 = tv3 = 0.0;
    VECTOR r1,r2,f1,f2;
    double scale = 0.5;
    energy = 0.0;
    if(Box!=NULL){ Box->get_vectors(tv1,tv2,tv3); scale = 1.0; }
    //-------------------------------------------------------------------------
    r1 = *(data_vdw->r1);
    r2 = *(data_vdw->r2) + (kx*tv1 + ky*tv2 + kz*tv3);
    if(potential_type=="LJ"){ energy += Vdw_LJ(r1,r2,f1,f2,data_vdw->sigma,data_vdw->epsilon);}
    else if(potential_type=="Buffered14_7"){ energy += Vdw_Buffered14_7(r1,r2,f1,f2,data_vdw->sigma,data_vdw->epsilon); }
    else if(potential_type=="Morse"){ energy += Vdw_Morse(r1,r2,f1,f2,data_vdw->D,data_vdw->r0,data_vdw->alpha); }

    *(data_vdw->f1) += scale*f1;
    *(data_vdw->f2) += scale*f2;
    //--------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    r1 = *(data_vdw->r1) + (kx*tv1 + ky*tv2 + kz*tv3);
    r2 = *(data_vdw->r2);
    f1 = f2 = 0.0;
    if(potential_type=="LJ"){ energy += Vdw_LJ(r1,r2,f1,f2,data_vdw->sigma,data_vdw->epsilon);}
    else if(potential_type=="Buffered14_7"){ energy += Vdw_Buffered14_7(r1,r2,f1,f2,data_vdw->sigma,data_vdw->epsilon); }
    else if(potential_type=="Morse"){ energy += Vdw_Morse(r1,r2,f1,f2,data_vdw->D,data_vdw->r0,data_vdw->alpha); }

    *(data_vdw->f1) += scale*f1;
    *(data_vdw->f2) += scale*f2;
    //--------------------------------------------------------------------------
    energy *= scale;

  }// vdw
*/
  else if(int_type==4){
    VECTOR tv1,tv2,tv3; tv1 = tv2 = tv3 = 0.0;
    VECTOR f1,f2,f12; f1 = f2 = f12 = 0.0;
    double SW = 1.0; 
    VECTOR dSW; dSW = 0.0;
    double en = 0.0;
    double d12;
    if(Box==NULL){

    if(is_active){  // This also excludes self-pairs!
    // Central cell contribution
    //-------------------------------------------------------------------------
    VECTOR& r1 = *(data_vdw->r1);   
    VECTOR& r2 = *(data_vdw->r2);
    VECTOR rij = r1 - r2;


    SW = 1.0; dSW = 0.0;
    if(data_vdw->is_cutoff){ 
      double d12 = rij.length2();
      if(d12<=data_vdw->R_off2){
        if(d12>=data_vdw->R_on2){
          SWITCH(r1,r2,data_vdw->R_on,data_vdw->R_off,SW,dSW);
        }
      }else{ SW = 0.0; dSW = 0.0; }
    }

    if(SW>0.0){
      if(functional==0){  en = Vdw_LJ(r1,r2,f1,f2,data_vdw->sigma,data_vdw->scale*data_vdw->epsilon);     }      
      else if(functional==1){ en = Vdw_Buffered14_7(r1,r2,f1,f2,data_vdw->sigma,data_vdw->scale*data_vdw->epsilon); }
      else if(functional==2){ en = Vdw_Morse(r1,r2,f1,f2,data_vdw->scale*data_vdw->D,data_vdw->r0,data_vdw->alpha); }

      energy += SW*en;
      f12 = (SW*f1 - en*dSW);
      *(data_vdw->f1) += f12;
      *(data_vdw->f2) -= f12;

      tp.tensor_product((*(data_vdw->r1) - *(data_vdw->r2)) , f12);   stress_at += tp;
      tp.tensor_product((*(data_vdw->g1) - *(data_vdw->g2)) , f12);   stress_fr += tp;
//      tp.tensor_product((*(data_vdw->m1) - *(data_vdw->m2)) , f12);   stress_ml += tp;

    }// SW>0.0
    //--------------------------------------------------------------------------
    }// is active
    }// Box==NULL // Gamma-point only

    else if(Box!=NULL){ // PBC case
      VECTOR T;
      VECTOR r1,r2;
      VECTOR g1,g2,g3; // reciprocal vectors
      VECTOR min_tr;
      VECTOR rij = *(data_vdw->r1) - *(data_vdw->r2);
      int xshift,yshift,zshift; // on how many cell dimensions to translate
                                // atom 2 before considering its images

      Box->get_vectors(tv1,tv2,tv3); 
      Box->inverse().T().get_vectors(g1,g2,g3);
      int nx,ny,nz; nx = ny = nz = 1;
      double scale = 1.0;//0.5;
      double fscale = 1.0;
      double tscale = 1.0;
      if(data_vdw->r1==data_vdw->r2){ scale = 0.5; } // Scaling of self-pairs

      triple central_translation;
      vector<triple> images;
      int n_images;

// No Verlet list
//        Cell cl(tv1,tv2,tv3,data_vdw->R_off); // shell = 2.0
//        cl.calculate(rij,images,central_translation);


// With Verlet list
/*
      data_vdw->time = data_vdw->time%5;
      if(data_vdw->time==0){
        Cell cl(tv1,tv2,tv3,data_vdw->R_off+2.0); // shell = 2.0
        cl.calculate(rij,images,central_translation);
        data_vdw->images = images;
        data_vdw->central_translation = central_translation;
      }else{      
        images = data_vdw->images;
        central_translation = data_vdw->central_translation;
      }
      data_vdw->time++;
*/

// Automatic Verlet list
    double R_skin = 2.0;
    double R_skin2 = 4.0;
    // Criterion
    // Calculate current central translation:
/*
    xshift = floor(rij*g1+0.5);
    yshift = floor(rij*g2+0.5);
    zshift = floor(rij*g3+0.5);
    if(!((xshift==data_vdw->central_translation.n1)&&
         (yshift==data_vdw->central_translation.n2)&&
         (zshift==data_vdw->central_translation.n3)
        )
    ){
      // Update list
      Cell cl(tv1,tv2,tv3,data_vdw->R_off+R_skin); // shell = 2.0
      cl.calculate(rij,images,central_translation);
      data_vdw->images = images;
      data_vdw->central_translation = central_translation;
      // Current variables become old:
      Box_old = *Box;
      data_vdw->r1_old = *data_vdw->r1;
      data_vdw->r2_old = *data_vdw->r2;
    }// central translation is different - recalculate anyways
    else{
*/
  
    // In this case let us check the displacements
//      double d2 = (rij.length2() - (data_vdw->r1_old - data_vdw->r2_old).length2());
//      VECTOR tv1_old,tv2_old,tv3_old;
//      Box_old.get_vectors(tv1_old,tv2_old,tv3_old);
/*
      double a2,b2,c2;
      a2 = tv1.length2();   
      b2 = tv2.length2();
      c2 = tv3.length2();

      double mind = a2;
      if(b2<=mind) { mind = b2; }
      if(c2<=mind) { mind = c2; }
      int Nmax = ceil(((data_vdw->R_off+R_skin) * (data_vdw->R_off+R_skin))/(mind))+1;
      d2 += (a2 - tv1_old.length2() + b2 - tv2_old.length2() + c2 - tv3_old.length2())*Nmax;
*/
/*
// Stupid and inefficient criterion,but hopefully correct:
     int is_update = 0;
     if(data_vdw->images.size()==0){ is_update = 1; }
     else{
         VECTOR T = ((tv1-tv1_old) + (tv2-tv2_old) + (tv3-tv3_old));
         double displ = (rij - (data_vdw->r1_old - data_vdw->r2_old) - T ).length();
         if(displ>R_skin){ is_update = 1; }
     }
*/
// Hopefully more efficient criterion and hopefully correct:
/*
     int is_update = 0;
     if(data_vdw->images.size()==0){ is_update = 1; }
     else{
         VECTOR T = ((tv1-tv1_old) + (tv2-tv2_old) + (tv3-tv3_old));
         double displ  = (*(data_vdw->r1) - data_vdw->r1_old).length2();
                displ += (*(data_vdw->r2) - data_vdw->r2_old).length2();
                displ += T.length2();
         if(displ>R_skin*R_skin){ is_update = 1; }
     }
*/
// Even more efficient version

     int is_update = 0;
     if(data_vdw->images.size()==0){ is_update = 1; }
     else{
// This was too much - overestimated
//         double displ2  = (*data_vdw->displr1_2) * (*data_vdw->displr1_2);
//                displ2 += (*data_vdw->displr2_2) * (*data_vdw->displr2_2);
//                displ2 += (*data_vdw->displT_2) * (*data_vdw->displT_2); 
// This is correct:
         // Calculate new central translation:
         int nxshift = floor(rij*g1+0.5);
         int nyshift = floor(rij*g2+0.5);
         int nzshift = floor(rij*g3+0.5);
         if((nxshift!=data_vdw->central_translation.n1)||
            (nyshift!=data_vdw->central_translation.n2)||
            (nzshift!=data_vdw->central_translation.n3)
           ){ is_update = 1; }
         else{
           double displ2  = (*data_vdw->displr1_2) + (*data_vdw->displr2_2) + (*data_vdw->displT_2);
           if(displ2>R_skin2){ is_update = 1;  }
         }// else
     }// else

//   is_update = 1; // this is no Verlet list option - for checking

      if(is_update){
        update_displ2 = 1;
        // Update list
        Cell cl(tv1,tv2,tv3,data_vdw->R_off+R_skin); // shell = 2.0
        cl.calculate(rij,images,central_translation);        
        data_vdw->images = images;
        data_vdw->central_translation = central_translation;
        // Current variables become old:
        Box_old = *Box;
        data_vdw->r1_old = *data_vdw->r1;
        data_vdw->r2_old = *data_vdw->r2;
      }
      else{
        // Use existing list
        images = data_vdw->images;
        central_translation = data_vdw->central_translation;
    }
//    }// central translation is the same


     
      
      n_images = images.size();
      xshift = central_translation.n1;
      yshift = central_translation.n2;
      zshift = central_translation.n3;

//      xshift = floor(rij*g1+0.5);
//      yshift = floor(rij*g2+0.5);
//      zshift = floor(rij*g3+0.5);
      for(int im=0;im<n_images;im++){

//      for(int Kx=(-nx+xshift);Kx<=(nx+xshift);Kx++){
//        for(int Ky=(-ny+yshift);Ky<=(ny+yshift);Ky++){
//          for(int Kz=(-nz+zshift);Kz<=(nz+zshift);Kz++){
            int Kx = images[im].n1;
            int Ky = images[im].n2;
            int Kz = images[im].n3;         

            T = (Kx*tv1 + Ky*tv2 + Kz*tv3);

            if((Kx==xshift) && (Ky==yshift) && (Kz==zshift)){ 
            // This is gamma point

              if(is_active){  // This also excludes self-pairs!
                // Central cell contribution
                //-------------------------------------------------------------------------
                r1 = *(data_vdw->r1);
                r2 = *(data_vdw->r2) + T;

                SW = 1.0; dSW = 0.0;
                VECTOR rij = r1 - r2;

                if(data_vdw->is_cutoff){
                  double d12 = rij.length2();
                  if(d12<=data_vdw->R_off2){
                    if(d12>=data_vdw->R_on2){
                      SWITCH(r1,r2,data_vdw->R_on,data_vdw->R_off,SW,dSW);
                    }
                  }else{ SW = 0.0; dSW = 0.0; }
                }

                if(SW>0.0){
                  if(functional==0){  en = Vdw_LJ(r1,r2,f1,f2,data_vdw->sigma,data_vdw->scale*data_vdw->epsilon);     }
                  else if(functional==1){ en = Vdw_Buffered14_7(r1,r2,f1,f2,data_vdw->sigma,data_vdw->scale*data_vdw->epsilon); }
                  else if(functional==2){ en = Vdw_Morse(r1,r2,f1,f2,data_vdw->scale*data_vdw->D,data_vdw->r0,data_vdw->alpha); }

                  energy += SW*en;
                  f12 = (SW*f1 - en*dSW);
                  *(data_vdw->f1) += f12;
                  *(data_vdw->f2) -= f12;
 
                  tp.tensor_product((*(data_vdw->r1) - *(data_vdw->r2) - T) , f12);   stress_at += tp;
                  tp.tensor_product((*(data_vdw->g1) - *(data_vdw->g2) - T) , f12);   stress_fr += tp;
                  //tp.tensor_product((*(data_vdw->m1) - *(data_vdw->m2)) , f12);   stress_ml += tp;
                }// SW>0.0
                //--------------------------------------------------------------------------
              }// is active            
            }// if gamma-point

            else{ // all other k-points
              //-------------------------------------------------------------------------
              r1 = *(data_vdw->r1);
              r2 = *(data_vdw->r2) + T;
              SW = 1.0; dSW = 0.0;
              VECTOR rij = r1 - r2;

              if(data_vdw->is_cutoff){
                double d12 = rij.length2();
                if(d12<=data_vdw->R_off2){
                  if(d12>=data_vdw->R_on2){
                    SWITCH(r1,r2,data_vdw->R_on,data_vdw->R_off,SW,dSW);
                  }
                }else{ SW = 0.0; dSW = 0.0; }
              }

             if(SW>0.0){
               f1 = f2 = 0.0;
               if(functional==0){ en = Vdw_LJ(r1,r2,f1,f2,data_vdw->sigma,scale*data_vdw->epsilon);}
               else if(functional==1){ en = Vdw_Buffered14_7(r1,r2,f1,f2,data_vdw->sigma,scale*data_vdw->epsilon); }
               else if(functional==2){ en = Vdw_Morse(r1,r2,f1,f2,scale*data_vdw->D,data_vdw->r0,data_vdw->alpha); }
               energy += SW*en;
               f12 = fscale*(SW*f1 - en*dSW);
               *(data_vdw->f1) += f12;
               *(data_vdw->f2) -= f12; // replica does not contribute to force

               tp.tensor_product((*(data_vdw->r1) - *(data_vdw->r2) - T) , f12);   stress_at += tscale*tp;
               tp.tensor_product((*(data_vdw->g1) - *(data_vdw->g2) - T) , f12);   stress_fr += tscale*tp;
               //tp.tensor_product((*(data_vdw->m1) - *(data_vdw->m2) - T) , f12);   stress_ml += tscale*tp;
             }
             //--------------------------------------------------------------------------
            }// else not Kx==0 && Ky==0 && Kz==0

//          }// Kz
//        }// Ky
//      }// Kx
      }// for im
    }// PBC

  }// vdw

  else if((int_type==5) && (functional==0)){
    VECTOR tv1,tv2,tv3; tv1 = tv2 = tv3 = 0.0;
    VECTOR f1,f2,f12; f1 = f2 = f12 = 0.0;
    double SW = 1.0;
    VECTOR dSW; dSW = 0.0;
    double en = 0.0;
    double d12;
    if(Box==NULL){

    if(is_active){  // This also excludes self-pairs!
    // Central cell contribution
    //-------------------------------------------------------------------------
    VECTOR& r1 = *(data_elec->r1);
    VECTOR& r2 = *(data_elec->r2);
    VECTOR rij = r1 - r2;


    SW = 1.0; dSW = 0.0;
    if(data_elec->is_cutoff){
      double d12 = rij.length2();
      if(d12<=data_elec->R_off2){
        if(d12>=data_elec->R_on2){
          SWITCH(r1,r2,data_elec->R_on,data_elec->R_off,SW,dSW);
        }
      }else{ SW = 0.0; dSW = 0.0; }
    }

    if(SW>0.0){
      if(functional==0){ en = Elec_Coulomb(r1,r2,f1,f2,data_elec->q1,data_elec->q2,data_elec->eps,data_elec->delta);}

      energy += SW*en;
      f12 = (SW*f1 - en*dSW);
      *(data_elec->f1) += f12;
      *(data_elec->f2) -= f12;

      tp.tensor_product((*(data_elec->r1) - *(data_elec->r2)) , f12);   stress_at += tp;
      tp.tensor_product((*(data_elec->g1) - *(data_elec->g2)) , f12);   stress_fr += tp;
//      tp.tensor_product((*(data_vdw->m1) - *(data_vdw->m2)) , f12);   stress_ml += tp;

    }// SW>0.0
    //--------------------------------------------------------------------------
    }// is active
    }// Box==NULL // Gamma-point only

    else if(Box!=NULL){ // PBC case
      VECTOR T;
      VECTOR r1,r2;
      VECTOR g1,g2,g3; // reciprocal vectors
      VECTOR min_tr;
      VECTOR rij = *(data_elec->r1) - *(data_elec->r2);
      int xshift,yshift,zshift; // on how many cell dimensions to translate
                                // atom 2 before considering its images

      Box->get_vectors(tv1,tv2,tv3);
      Box->inverse().T().get_vectors(g1,g2,g3);
      int nx,ny,nz; nx = ny = nz = 1;
      double scale = 1.0;//0.5;
      double fscale = 1.0;
      double tscale = 1.0;
      if(data_elec->r1==data_elec->r2){ scale = 0.5; } // Scaling of self-pairs

      xshift = floor(rij*g1+0.5);
      yshift = floor(rij*g2+0.5);
      zshift = floor(rij*g3+0.5);

      for(int Kx=(-nx+xshift);Kx<=(nx+xshift);Kx++){
        for(int Ky=(-ny+yshift);Ky<=(ny+yshift);Ky++){
          for(int Kz=(-nz+zshift);Kz<=(nz+zshift);Kz++){

            T = (Kx*tv1 + Ky*tv2 + Kz*tv3);

            if((Kx==xshift) && (Ky==yshift) && (Kz==zshift)){
            // This is gamma point

              if(is_active){  // This also excludes self-pairs!
                // Central cell contribution
                //-------------------------------------------------------------------------
                r1 = *(data_elec->r1);
                r2 = *(data_elec->r2) + T;

                SW = 1.0; dSW = 0.0;
                VECTOR rij = r1 - r2;

                if(data_elec->is_cutoff){
                  double d12 = rij.length2();
                  if(d12<=data_elec->R_off2){
                    if(d12>=data_elec->R_on2){
                      SWITCH(r1,r2,data_elec->R_on,data_elec->R_off,SW,dSW);
                    }
                  }else{ SW = 0.0; dSW = 0.0; }
                }

                if(SW>0.0){
                  if(functional==0){ en = Elec_Coulomb(r1,r2,f1,f2,data_elec->q1,data_elec->q2,data_elec->eps,data_elec->delta);}

                  energy += SW*en;
                  f12 = (SW*f1 - en*dSW);
                  *(data_elec->f1) += f12;
                  *(data_elec->f2) -= f12;

                  tp.tensor_product((*(data_elec->r1) - *(data_elec->r2) - T) , f12);   stress_at += tp;
                  tp.tensor_product((*(data_elec->g1) - *(data_elec->g2) - T) , f12);   stress_fr += tp;
                  //tp.tensor_product((*(data_vdw->m1) - *(data_vdw->m2)) , f12);   stress_ml += tp;
                }// SW>0.0
                //--------------------------------------------------------------------------
              }// is active
            }// if gamma-point
            else{ // all other k-points
              //-------------------------------------------------------------------------
              r1 = *(data_elec->r1);
              r2 = *(data_elec->r2) + T;
              SW = 1.0; dSW = 0.0;
              VECTOR rij = r1 - r2;

              if(data_elec->is_cutoff){
                double d12 = rij.length2();
                if(d12<=data_elec->R_off2){
                  if(d12>=data_elec->R_on2){
                    SWITCH(r1,r2,data_elec->R_on,data_elec->R_off,SW,dSW);
                  }
                }else{ SW = 0.0; dSW = 0.0; }
              }

             if(SW>0.0){
               f1 = f2 = 0.0;
               if(functional==0){ en = Elec_Coulomb(r1,r2,f1,f2,data_elec->q1,data_elec->q2,data_elec->eps,data_elec->delta);}

               energy += SW*en;
               f12 = fscale*(SW*f1 - en*dSW);
               *(data_elec->f1) += f12;
               *(data_elec->f2) -= f12; // replica does not contribute to force

               tp.tensor_product((*(data_elec->r1) - *(data_elec->r2) - T) , f12);   stress_at += tscale*tp;
               tp.tensor_product((*(data_elec->g1) - *(data_elec->g2) - T) , f12);   stress_fr += tscale*tp;
               //tp.tensor_product((*(data_vdw->m1) - *(data_vdw->m2) - T) , f12);   stress_ml += tscale*tp;
             }
             //--------------------------------------------------------------------------
            }// else not Kx==0 && Ky==0 && Kz==0
          }// Kz
        }// Ky
      }// Kx
    }// PBC

  }// elec


/*
  else if((int_type==5) && (functional==0)){
    VECTOR tv1,tv2,tv3; tv1 = tv2 = tv3 = 0.0;
    VECTOR r1,r2,f1,f2,f12; f1 = f2 = 0.0;
    double SW = 1.0;
    VECTOR dSW; dSW = 0.0;
    double en = 0.0;
    if(Box==NULL){ // Only gamma-point
      if(is_active){  // This also excludes self-pairs!
      // Central cell contribution
      //-------------------------------------------------------------------------
      r1 = *(data_elec->r1);   r2 = *(data_elec->r2);
      if(data_elec->is_cutoff){ SWITCH(r1,r2,data_elec->R_on,data_elec->R_off,SW,dSW); }

      if(SW>0.0){
        if(functional==0){ en = Elec_Coulomb(r1,r2,f1,f2,data_elec->q1,data_elec->q2,data_elec->eps,data_elec->delta);}

      en = data_elec->scale*en;
      energy += SW*en;
      *(data_elec->f1) += (SW*f1 - en*dSW);
      *(data_elec->f2) += (SW*f2 + en*dSW);
    }
    //--------------------------------------------------------------------------
    }// is_active
    }// Box==NULL

    else if(Box!=NULL){ // PBC case
      VECTOR T;
      Box->get_vectors(tv1,tv2,tv3);
      int nx,ny,nz; nx = ny = nz = 1;
      double scale = 1.0;
      if(data_elec->r1==data_elec->r2){ scale = 0.5; } // Scaling of self-pairs
//      else{ scale = 1.0; }

      for(int Kx=-nx;Kx<=nx;Kx++){
        for(int Ky=-ny;Ky<=ny;Ky++){
          for(int Kz=-nz;Kz<=nz;Kz++){
            if(Kx==0 && Ky==0 && Kz==0){ }
            else{

      T = (Kx*tv1 + Ky*tv2 + Kz*tv3);
      //-------------------------------------------------------------------------
      r1 = *(data_elec->r1);
      r2 = *(data_elec->r2) + T;
      SW = 1.0; dSW = 0.0;

      if(data_elec->is_cutoff){ SWITCH(r1,r2,data_elec->R_on,data_elec->R_off,SW,dSW); }

      if(SW>0.0){
        f1 = f2 = 0.0;
        if(functional==0){ en = Elec_Coulomb(r1,r2,f1,f2,data_elec->q1,data_elec->q2,data_elec->eps,data_elec->delta);}

        energy += scale*SW*en;
        f12 = scale*(SW*f1 - en*dSW);
        *(data_elec->f1) += f12;
        *(data_elec->f2) -= f12; // = +=scale*(SW*f2 + en*dSW); // replica does not contribute to force
      }
      //--------------------------------------------------------------------------
      }// else
     }// Kz
     }// Ky
     }// Kx
   }// Box!=NULL

  }// elec
*/

  else if(int_type==6){
    double en;
    int sz = data_mb->sz;
    VECTOR* r;  r = new VECTOR[sz];
    VECTOR* g;  g = new VECTOR[sz];
    VECTOR* m;  m = new VECTOR[sz];
    VECTOR* f;  f = new VECTOR[sz];
    double* q;  q = new double[sz];
    double* epsilon; epsilon = new double[sz];
    double* sigma;   sigma = new double[sz];
    double* dr2;dr2 = new double[sz];
    for(int i=0;i<sz;i++){
      r[i] = *(data_mb->r[i]);
      g[i] = *(data_mb->g[i]);
      m[i] = *(data_mb->m[i]);
      q[i] = *(data_mb->q[i]);
      epsilon[i] = *(data_mb->epsilon[i]);
      sigma[i] = *(data_mb->sigma[i]);
      f[i] = 0.0;
      dr2[i] = *(data_mb->displr_2[i]);      
    }

    int rec_deg, pbc_deg;
    pbc_deg = 1;
    rec_deg = 1;
    int is_cutoff = 1;
    double R_on, R_off; R_on = 8.0; R_off = 10.0;
    if(data_mb->is_cutoff){ R_on = data_mb->R_on; R_off = data_mb->R_off; }
    int is_update = 0;
    MATRIX3x3 at_st, fr_st,ml_st;
//    cout<<"dr2 = "<<dr2<<endl;
//    cout<<"(data_mb->displT_2) = "<<(data_mb->displT_2)<<endl;
//    cout<<"*(data_mb->displT_2) = "<<*(data_mb->displT_2)<<endl;
//    exit(0);
//    cout<<"in interaction of type = "<<int_type<<endl;
//    cout<<"functional = "<<functional<<endl;
//    vector< vector<quartet> > at_neib;

    if(functional==0){
    en = Elec_Ewald3D(r,g,m,f,at_st,fr_st,ml_st,sz,q,data_mb->nexcl,data_mb->excl1,data_mb->excl2,data_mb->scale,Box,rec_deg,pbc_deg,data_mb->elec_etha,is_cutoff,R_on,R_off,data_mb->time,data_mb->images,data_mb->central_translation,dr2,*(data_mb->displT_2),is_update); 

    }
    else if(functional==1){
    en = Vdw_LJ(r,g,m,f,at_st,fr_st,ml_st,sz,epsilon,sigma,data_mb->nexcl,data_mb->excl1,data_mb->excl2,data_mb->scale,Box,rec_deg,pbc_deg,data_mb->elec_etha,is_cutoff,R_on,R_off,data_mb->time,data_mb->images,data_mb->central_translation,dr2,*(data_mb->displT_2),is_update);
    }

    else if(functional==2){
//    en = Vdw_LJ1(r,g,m,f,at_st,fr_st,ml_st,sz,epsilon,sigma,data_mb->nexcl,data_mb->excl1,data_mb->excl2,data_mb->scale,Box,rec_deg,pbc_deg,data_mb->elec_etha,is_cutoff,R_on,R_off,data_mb->time,data_mb->at_neib,data_mb->central_translation,dr2,*(data_mb->displT_2),is_update);

//    cout<<"data_mb->excl_scales.size = "<<data_mb->excl_scales.size()<<endl;
      en = Vdw_LJ2_no_excl(r,g,m,f,at_st,fr_st,ml_st,sz,epsilon,sigma,data_mb->nexcl,data_mb->excl1,data_mb->excl2,data_mb->scale,Box,rec_deg,pbc_deg,data_mb->elec_etha,is_cutoff,R_on,R_off,data_mb->time,data_mb->excl_scales);
      is_update = 1;

    }


    energy += en;
    stress_at += at_st;
    stress_fr += fr_st;
    stress_ml += ml_st;

    if(is_update){ update_displ2 = 1; }

    for(i=0;i<sz;i++){   *(data_mb->f[i]) += f[i];   }

    delete [] r;
    delete [] g;
    delete [] m;
    delete [] f;
    delete [] q;
    delete [] epsilon;
    delete [] sigma;
    delete [] dr2;
 
  }// mb && Ewald_3D

  else if((int_type==7) && (functional==0)){
    if(is_active){
    VECTOR& r1 = *(data_gay_berne->r1); VECTOR& r2 = *(data_gay_berne->r2);
    VECTOR& u1 = *(data_gay_berne->u1); VECTOR& u2 = *(data_gay_berne->u2);
    VECTOR f1,f2,t1,t2;
    energy = Gay_Berne(r1,r2,u1,u2,f1,f2,t1,t2,data_gay_berne->di,data_gay_berne->dj,data_gay_berne->li,data_gay_berne->lj,
             data_gay_berne->e0,data_gay_berne->rat,data_gay_berne->dw,data_gay_berne->mu,data_gay_berne->nu
             );
    *(data_gay_berne->f1) += f1;
    *(data_gay_berne->f2) += f2;
    *(data_gay_berne->t1) += t1;
    *(data_gay_berne->t2) += t2;
    }
  }// gay-berne

  else if(int_type==8){
    double en;
    int sz = data_mb->sz;
    VECTOR* r;  r = new VECTOR[sz];
    VECTOR* g;  g = new VECTOR[sz];
    VECTOR* m;  m = new VECTOR[sz];
    VECTOR* f;  f = new VECTOR[sz];
    double* q;  q = new double[sz];
    double* epsilon; epsilon = new double[sz];
    double* sigma;   sigma = new double[sz];
    double* dr2;dr2 = new double[sz];
    for(int i=0;i<sz;i++){
      r[i] = *(data_mb->r[i]);
      g[i] = *(data_mb->g[i]);
      m[i] = *(data_mb->m[i]);
      q[i] = *(data_mb->q[i]);
      epsilon[i] = *(data_mb->epsilon[i]);
      sigma[i] = *(data_mb->sigma[i]);
      f[i] = 0.0;
      dr2[i] = *(data_mb->displr_2[i]);
    }

    int rec_deg, pbc_deg;
    pbc_deg = 1;
    rec_deg = 1;
    int is_cutoff = 1;
    double R_on, R_off; R_on = 8.0; R_off = 10.0;
    if(data_mb->is_cutoff){ R_on = data_mb->R_on; R_off = data_mb->R_off; }
    int is_update = 0;
    MATRIX3x3 at_st, fr_st,ml_st;

    if(functional==0){
//    en = Vdw_LJ1(r,g,m,f,at_st,fr_st,ml_st,sz,epsilon,sigma,data_mb->nexcl,data_mb->excl1,data_mb->excl2,data_mb->scale,Box,rec_deg,pbc_deg,data_mb->elec_etha,is_cutoff,R_on,R_off,data_mb->time,data_mb->at_neib,data_mb->central_translation,dr2,*(data_mb->displT_2),is_update);

//    cout<<"data_mb->excl_scales.size = "<<data_mb->excl_scales.size()<<endl;
      en = Vdw_LJ2_excl(r,g,m,f,at_st,fr_st,ml_st,sz,epsilon,sigma,data_mb->nexcl,data_mb->excl1,data_mb->excl2,data_mb->scale,Box,rec_deg,pbc_deg,data_mb->elec_etha,is_cutoff,R_on,R_off,data_mb->time,data_mb->excl_scales);
      is_update = 1;

    }

    energy += en;
    stress_at += at_st;
    stress_fr += fr_st;
    stress_ml += ml_st;

    if(is_update){ update_displ2 = 1; }

    for(i=0;i<sz;i++){   *(data_mb->f[i]) += f[i];   }

    delete [] r;
    delete [] g;
    delete [] m;
    delete [] f;
    delete [] q;
    delete [] epsilon;
    delete [] sigma;
    delete [] dr2;

  }// mb exclusions


  }// if call_type==int_type
  return energy;
}

