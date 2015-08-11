#include "System.h"

using namespace libmmath::librandom;

namespace libchemobjects{
namespace libchemsys{




void System::cool_atoms(){
  for(int i=0;i<Number_of_atoms;i++){
    Atoms[i].Atom_RB.scale_angular_(0.0);
    Atoms[i].Atom_RB.scale_linear_(0.0);
  }
  zero_atom_forces();
}

void System::cool_fragments(){
  for(int i=0;i<Number_of_fragments;i++){
    Fragments[i].Group_RB.scale_angular_(0.0);
    Fragments[i].Group_RB.scale_linear_(0.0);
  }
  zero_forces_and_torques();
}

void System::cool(){
  cool_atoms();
  cool_fragments();
}


void System::zero_atom_forces(){
  is_stress_at = 0;
  is_stress_fr = 0;
  is_stress_ml = 0;
  for(int i=0;i<Number_of_atoms;i++){   Atoms[i].Atom_RB.rb_force = 0.0;  }
}

void System::zero_fragment_forces(){
  is_stress_at = 0;
  is_stress_fr = 0;
  is_stress_ml = 0;
  for(int i=0;i<Number_of_fragments;i++){ Fragments[i].Group_RB.rb_force = 0.0; } 
}

void System::zero_fragment_torques(){
  for(int i=0;i<Number_of_fragments;i++){ Fragments[i].Group_RB.rb_torque_e = 0.0; }
}

void System::zero_forces(){
  zero_atom_forces();
  zero_fragment_forces();
}

void System::zero_forces_and_torques(){
  zero_forces();
  zero_fragment_torques();
}

void System::update_fragment_forces(){
  for(int i=0;i<Number_of_fragments;i++){
    for(int j=0;j<Fragments[i].Group_Size;j++){
      int at_indx = Fragments[i].globAtom_Index[j];
      Fragments[i].Group_RB.rb_force += Atoms[at_indx].Atom_RB.rb_force;
    }
  }
}

void System::update_fragment_torques(){
  for(int i=0;i<Number_of_fragments;i++){
    for(int j=0;j<Fragments[i].Group_Size;j++){
      int at_indx = Fragments[i].globAtom_Index[j];
      Fragments[i].Group_RB.rb_torque_e += cross(1.0,Fragments[i].Group_RB.rb_centers[j],Fragments[i].Group_RB.rb_A_I_to_e*Atoms[at_indx].Atom_RB.rb_force);
    }
  }
}

void System::update_fragment_forces_and_torques(){
  for(int i=0;i<Number_of_fragments;i++){
    for(int j=0;j<Fragments[i].Group_Size;j++){
      int at_indx = Fragments[i].globAtom_Index[j];
      Fragments[i].Group_RB.rb_force += Atoms[at_indx].Atom_RB.rb_force;
      Fragments[i].Group_RB.rb_torque_e += cross(1.0,Fragments[i].Group_RB.rb_centers[j],Fragments[i].Group_RB.rb_A_I_to_e*Atoms[at_indx].Atom_RB.rb_force);
    }
  }

}

void System::save_forces(vector<VECTOR>& frcs){
// This function extracts current forces on fragments and saves them
// in external vector <frcs>
  VECTOR x(0.0,0.0,0.0);
  frcs = std::vector<VECTOR>(Number_of_fragments,x);
  for(int i=0;i<Number_of_fragments;i++){  frcs[i] = Fragments[i].Group_RB.rb_force; }
}

void System::save_torques(vector<VECTOR>& trcs){
// This function extracts current forces on fragments and saves them
// in external vector <trcs>
  VECTOR x(0.0,0.0,0.0);
  trcs = std::vector<VECTOR>(Number_of_fragments,x);
  for(int i=0;i<Number_of_fragments;i++){  trcs[i] = Fragments[i].Group_RB.rb_torque_e; }
}

void System::load_forces(vector<VECTOR>& frcs){
// This function assignes external forces <frcs> to
// current forces on fragments 
  if(frcs.size()==Number_of_fragments){
    for(int i=0;i<Number_of_fragments;i++){  Fragments[i].Group_RB.rb_force = frcs[i]; }
  }
}

void System::load_torques(vector<VECTOR>& trcs){
// This function assignes external torques <trcs> to
// current torques on fragments
  if(trcs.size()==Number_of_fragments){
    for(int i=0;i<Number_of_fragments;i++){  Fragments[i].Group_RB.rb_torque_e = trcs[i]; }
  }
}

void System::save_stress(MATRIX3x3& strs){
  strs = stress_fr;
}

void System::increment_stress(MATRIX3x3& strs){
  stress_fr += strs;
}

void System::save_respa_state(std::string respa_type){

// Former version: memory leak, memory-intensive
  VECTOR x(0.0,0.0,0.0);
  std::vector<VECTOR> frcs(Number_of_fragments,x);
  std::vector<VECTOR> trcs(Number_of_fragments,x);


  for(int i=0;i<Number_of_fragments;i++){
    frcs[i] = Fragments[i].Group_RB.rb_force;
    trcs[i] = Fragments[i].Group_RB.rb_torque_e;
  }
  if(respa_type=="fast"){ respa_f_fast = frcs; respa_t_fast = trcs; }
  else if(respa_type=="medium"){ respa_f_medium = frcs; respa_t_medium = trcs; }


/*
  if(respa_type=="fast"){ 
    for(int i=0;i<Number_of_fragments;i++){
      respa_f_fast[i] = Fragments[i].Group_RB.rb_force; 
      respa_t_fast[i] = Fragments[i].Group_RB.rb_torque_e;
    }
  }
  else if(respa_type=="medium"){
    for(int i=0;i<Number_of_fragments;i++){
      respa_f_medium[i] = Fragments[i].Group_RB.rb_force; 
      respa_t_medium[i] = Fragments[i].Group_RB.rb_torque_e;
    }
  }
 */
}

void System::load_respa_state(std::string respa_type){

  VECTOR x(0.0,0.0,0.0);
  std::vector<VECTOR> frcs(Number_of_fragments,x);
  std::vector<VECTOR> trcs(Number_of_fragments,x);

  if(respa_type=="fast"){
    if(respa_f_fast.size()==Number_of_fragments){ frcs = respa_f_fast; }
    if(respa_t_fast.size()==Number_of_fragments){ trcs = respa_t_fast; }
  }
  else if(respa_type=="medium"){
    if(respa_f_medium.size()==Number_of_fragments){ frcs = respa_f_medium; }
    if(respa_t_medium.size()==Number_of_fragments){ trcs = respa_t_medium; }
  }

  for(int i=0;i<Number_of_fragments;i++){
    Fragments[i].Group_RB.rb_force = frcs[i];
    Fragments[i].Group_RB.rb_torque_e = trcs[i];
  }


/*
  if(respa_type=="fast"){
    if(respa_f_fast.size()==Number_of_fragments){ 
      for(int i=0;i<Number_of_fragments;i++){
        Fragments[i].Group_RB.rb_force = respa_f_fast[i];
      }
    }
    if(respa_t_fast.size()==Number_of_fragments){ 
      for(int i=0;i<Number_of_fragments;i++){
        Fragments[i].Group_RB.rb_torque_e = respa_t_fast[i];
      }
    }
  }// fast

  else if(respa_type=="medium"){
    if(respa_f_fast.size()==Number_of_fragments){ 
      for(int i=0;i<Number_of_fragments;i++){
        Fragments[i].Group_RB.rb_force = respa_f_medium[i];
      }
    }
    if(respa_t_fast.size()==Number_of_fragments){ 
      for(int i=0;i<Number_of_fragments;i++){
        Fragments[i].Group_RB.rb_torque_e = respa_t_medium[i];
      }
    }

  }// medium

*/

}


void System::init_fragments(){
  Nf_t = 0;     is_Nf_t = 1;
  Nf_r = 0;     is_Nf_r = 1;
  for(int i=0;i<Number_of_fragments;i++){
    int sz = Fragments[i].Group_Size;
    double* masses; masses = new double[sz];
    VECTOR* coords; coords = new VECTOR[sz];
    for(int j=0;j<sz;j++){
      int at_indx = Fragments[i].globAtom_Index[j];
      masses[j] = Atoms[at_indx].Atom_RB.rb_mass;
      coords[j] = Atoms[at_indx].Atom_RB.rb_cm;
    }

    Fragments[i].Group_RB.init(sz,masses,coords);
    Fragments[i].is_Group_RB = 1;

    Nf_t += Fragments[i].Group_RB.get_Nf_t();
    Nf_r += Fragments[i].Group_RB.get_Nf_r();

    delete [] masses;
    delete [] coords;
  }
}

void System::init_molecules(){

  init_fragments();

  for(int i=0;i<Number_of_molecules;i++){
    int sz = Molecules[i].globAtom_Index.size(); //Molecule_Size;
    double* masses; masses = new double[sz];
    VECTOR* coords; coords = new VECTOR[sz];

    for(int j=0;j<sz;j++){
      int at_indx = Molecules[i].globAtom_Index[j];
      masses[j] = Atoms[at_indx].Atom_RB.rb_mass;
      coords[j] = Atoms[at_indx].Atom_RB.rb_cm;
    }

    Molecules[i].Molecule_RB.init(sz,masses,coords);
    Molecules[i].is_Molecule_RB = 1;
    delete [] masses;
    delete [] coords;
  }// for i
    
}

void System::init_box_origin(){
  Box_origin = Atoms[0].Atom_RB.rb_cm;
  for(int i=0;i<Number_of_atoms;i++){
    if(Atoms[i].Atom_RB.rb_cm.x<Box_origin.x){ Box_origin.x = Atoms[i].Atom_RB.rb_cm.x; }
    if(Atoms[i].Atom_RB.rb_cm.y<Box_origin.y){ Box_origin.y = Atoms[i].Atom_RB.rb_cm.y; }
    if(Atoms[i].Atom_RB.rb_cm.z<Box_origin.z){ Box_origin.z = Atoms[i].Atom_RB.rb_cm.z; }
  }
  is_Box_origin = 1;
}

void System::init_box(){
/******************************************************************
  This function automatically determines the size of the box
  based on the system size. The box will be cubical or parallelepiped
******************************************************************/
  double maxx, maxy,maxz;
  VECTOR rij;
  maxx = maxy = maxz = 0.0;
  for(int i=0;i<Number_of_atoms;i++){
    for(int j=i+1;j<Number_of_atoms;j++){
      rij = Atoms[i].Atom_RB.rb_cm - Atoms[j].Atom_RB.rb_cm;
      if(fabs(rij.x)>maxx) { maxx = fabs(rij.x); }
      if(fabs(rij.y)>maxy) { maxy = fabs(rij.y); }
      if(fabs(rij.z)>maxz) { maxz = fabs(rij.z); }      
    }// for j
  }// for i
  // Add some space so the atoms of neighbor cells do not overlap
  // or they are not too close to each other
  maxx += 5.0;
  maxy += 5.0;
  maxz += 5.0;

/* If we want cubic Box
  double maxL = maxx;
  maxL = (maxL>maxy)?maxL:maxy;
  maxL = (maxL>maxz)?maxL:maxz;
  maxx = maxy = maxz = maxL;
*/
  VECTOR tv1,tv2,tv3; 
  tv1.x = maxx; tv1.y = 0.0;  tv1.z = 0.0;   
  tv2.x = 0.0;  tv2.y = maxy; tv2.z = 0.0; 
  tv3.x = 0.0;  tv3.y = 0.0;  tv3.z = maxz; 

  Box.init(tv1,tv2,tv3); is_Box = 1;
  Boxold.init(tv1,tv2,tv3); is_Boxold = 1;

  if(!is_Box_origin){ init_box_origin(); }

//  apply_frag_pbc("abc");

/*
  //------------- Update interactions ----------------
  int n_inter = interactions.size();
  for(i=0;i<n_inter;i++){
    interactions[i].set_pbc(&Box,0,0,0);
  }
*/

  //cout<<"is_Box = "<<is_Box<<endl;
}

void System::init_box(VECTOR tv1,VECTOR tv2,VECTOR tv3){
/******************************************************************
  This function creats the box based on argument vectors
******************************************************************/
  Box.init(tv1,tv2,tv3); is_Box = 1;
  Boxold.init(tv1,tv2,tv3); is_Boxold = 1;
  if(!is_Box_origin){ init_box_origin(); }

//  apply_frag_pbc("abc");

  //------------- Update interactions ----------------
//  int n_inter = interactions.size();
//  for(int i=0;i<n_inter;i++){
//    interactions[i].set_pbc(&Box,0,0,0);
//  }

}

void System::init_box(double maxx,double maxy,double maxz){
/******************************************************************
  This function creates the parralelepiped box with dimensions
  given by argument
******************************************************************/
  VECTOR tv1,tv2,tv3;
  tv1.x = maxx; tv1.y = 0.0;  tv1.z = 0.0;
  tv2.x = 0.0;  tv2.y = maxy; tv2.z = 0.0;
  tv3.x = 0.0;  tv3.y = 0.0;  tv3.z = maxz;

  Box.init(tv1,tv2,tv3); is_Box = 1;
  Boxold.init(tv1,tv2,tv3); is_Boxold = 1;
  if(!is_Box_origin){ init_box_origin(); }

//  apply_frag_pbc("abc");

  //------------- Update interactions ----------------
//  int n_inter = interactions.size();
//  for(int i=0;i<n_inter;i++){
//    interactions[i].set_pbc(&Box,0,0,0);
//  }

}

/*
void System::apply_atom_pbc(std::string pbc_type){
****************************************************************
 This function shifts the centers of mass of all atoms such
 that they are inside of the box [0,X] x [0,Y] x [0,Z]
*****************************************************************
  if(is_Box){
  MATRIX invBox;
  Box.Inverse(invBox);
  if(!is_Box_origin) {init_box_origin(); is_Box_origin = 1;}

  for(int i=0;i<Number_of_atoms;i++){
    RigidBody& top = Atoms[i].Atom_RB;
    VECTOR r = invBox*(top.rb_cm - Box_origin);

    if((pbc_type=="a")||(pbc_type=="ab")||(pbc_type=="ac")||(pbc_type=="abc")) {
      r.x = r.x - floor(r.x);
    }
    if((pbc_type=="b")||(pbc_type=="ab")||(pbc_type=="bc")||(pbc_type=="abc")) {
      r.y = r.y - floor(r.y);
    }
    if((pbc_type=="c")||(pbc_type=="ac")||(pbc_type=="bc")||(pbc_type=="abc")) {
      r.z = r.z - floor(r.z);
    }
    top.rb_cm = Box * r + Box_origin;
  }// for i
  }
}
*/

void System::apply_frag_pbc(std::string pbc_type){
/****************************************************************
 This function shifts the centers of mass of all fragments such
 that they are inside of the box [0,X] x [0,Y] x [0,Z]
 Because the basic unit of all simulations is rigid body(fragment)
 the PBC will be applied to it. This will cure the problem of several
 molecules grouped into one fragment. However, if the molecule is
 build up of several fragments, the special care should be taken for
 bonded interactions (bonds, angles, dihedrals,exclusions, etc.)!
*****************************************************************/
  if(is_Box){
  MATRIX3x3 invBox;
  invBox = Box.inverse();
  if(!is_Box_origin) {init_box_origin(); is_Box_origin = 1;}

  for(int i=0;i<Number_of_fragments;i++){
    RigidBody& top = Fragments[i].Group_RB;
    VECTOR r = invBox*(top.rb_cm - Box_origin);

    if((pbc_type=="a")||(pbc_type=="ab")||(pbc_type=="ac")||(pbc_type=="abc")) {
      r.x = r.x - floor(r.x);
    }
    if((pbc_type=="b")||(pbc_type=="ab")||(pbc_type=="bc")||(pbc_type=="abc")) {
      r.y = r.y - floor(r.y);
    }
    if((pbc_type=="c")||(pbc_type=="ac")||(pbc_type=="bc")||(pbc_type=="abc")) {
      r.z = r.z - floor(r.z);
    }
    top.rb_cm = Box * r + Box_origin;
  }// for i
  }
}

/*
void System::apply_mol_pbc(std::string pbc_type){
****************************************************************
 This function shifts the centers of mass of all molecules such
 that they are inside of the box [0,X] x [0,Y] x [0,Z]
*****************************************************************
  if(is_Box){
  MATRIX invBox;
  Box.Inverse(invBox);
  if(!is_Box_origin) {init_box_origin(); is_Box_origin = 1;}

  for(int i=0;i<Number_of_molecules;i++){
    RigidBody& top = Molecules[i].Molecule_RB;
    VECTOR r = invBox*(top.rb_cm - Box_origin);

    if((pbc_type=="a")||(pbc_type=="ab")||(pbc_type=="ac")||(pbc_type=="abc")) {
      r.x = r.x - floor(r.x);
    }
    if((pbc_type=="b")||(pbc_type=="ab")||(pbc_type=="bc")||(pbc_type=="abc")) {
      r.y = r.y - floor(r.y);
    }
    if((pbc_type=="c")||(pbc_type=="ac")||(pbc_type=="bc")||(pbc_type=="abc")) {
      r.z = r.z - floor(r.z);
    }
    top.rb_cm = Box * r + Box_origin;
  }// for i

  }// if is_Box
}
*/

/*
void System::apply_mol_pbc(std::string pbc_type){
**********************************************************
 This function shifts the centers of mass of molecules
 (and hence fragments) such that they are inside of
 the box [0,X] x [0,Y] x [0,Z]
************************************************************
  MATRIX invBox;
  Box.Inverse(invBox);
  if(!is_Box_origin) {init_box_origin(); }

  double x,y,z;
  int is_translated = 0;
  double M = 0.0;
  VECTOR R; R = 0.0;
  VECTOR O; O = 0.0;

  // Update center of mass of the molecules and masses
  for(int i=0;i<Number_of_molecules;i++){
    RigidBody& mol = ObjectSpace_Molecules[i].Molecule_RB;
    double& mM = mol.rb_mass;
    VECTOR& mR = mol.rb_cm
    mM = 0.0;
    mR = 0.0;
    for(int j=0;j<ObjectSpace_Molecules[i].globGroup_Index.size();j++){
      int grp_indx = os.ObjectSpace_Molecules[i].globGroup_Index[j];
      RigidBody& top = os.ObjectSpace_Fragments[grp_indx].Fragment_RB;
      double m =   top.rb_mass;

      mR += (m*top.rb_cm);
      mM += m;
    }// for j

    R += mR;
    M += mM;
    mR = mR/mM;
  }// for i

  R = R/M;
  // Apply PBC to centers of mass of the molecules
  // and modify centers of mass of corresponding groups
  for(i=0;i<ObjectSpace_Number_of_molecules;i++){
    RigidBody& mol = ObjectSpace_Molecules[i].Molecule_RB;
    double& mM = mol.rb_mass;
    VECTOR& mR = mol.rb_cm;
    VECTOR s,s0;

    s.x = ((mR-ObjectSpace_box_origin) * ObjectSpace_g1);
    s.y = ((mR-ObjectSpace_box_origin) * ObjectSpace_g2);
    s.z = ((mR-ObjectSpace_box_origin) * ObjectSpace_g3);
    s0 = s;

    if((pbc_type=="a")||(pbc_type=="ab")||(pbc_type=="ac")||(pbc_type=="abc")) {
      s.x = s.x - floor(s.x);
    }
    if((pbc_type=="b")||(pbc_type=="ab")||(pbc_type=="bc")||(pbc_type=="abc")) {
      s.y = s.y - floor(s.y);
    }
    if((pbc_type=="c")||(pbc_type=="ac")||(pbc_type=="bc")||(pbc_type=="abc")) {
      s.z = s.z - floor(s.z);
    }

    if(s0==s){  is_translated = 0; }
    else{       is_translated = 1; }

    VECTOR dmR; dmR = mR;
    mR = s.x * ObjectSpace_tv1 + s.y * ObjectSpace_tv2 + s.z * ObjectSpace_tv3 + ObjectSpace_box_origin;
    dmR = mR - dmR;// this is molecule's displacement due to PDB
    for(int j=0;j<os.ObjectSpace_Molecules[i].globGroup_Index.size();j++){
      int grp_indx = os.ObjectSpace_Molecules[i].globGroup_Index[j];
      ObjectSpace_Fragments[grp_indx].Fragment_RB.rb_cm += dmR;
    }// for j
  }// for i
}
*/


void System::fix_fragment_translation(int gr_id){
  int indx = get_fragment_index_by_fragment_id(gr_id);
  if(indx>-1){ 
    cout<<"Fixing fragment with indx = "<<indx<<". Translational DOFs\n";
    Fragments[indx].Group_RB.fix_translation();
    Nf_t -= 3;
  }
}

void System::fix_fragment_rotation(int gr_id){
  int indx = get_fragment_index_by_fragment_id(gr_id);
  if(indx>-1){
    cout<<"Fixing fragment with indx = "<<indx<<". Rotational DOFs\n";
    Fragments[indx].Group_RB.fix_rotation();
    Nf_r -= 3;
  }
}

void System::fix_fragment(int gr_id){
  int indx = get_fragment_index_by_fragment_id(gr_id);
  if(indx>-1){ 
    cout<<"Fixing fragment with indx = "<<indx<<". Translational DOFs\n";
    Fragments[indx].Group_RB.fix_translation();
    Nf_t -= 3;
  }
  if(indx>-1){
    cout<<"Fixing fragment with indx = "<<indx<<". Rotational DOFs\n";
    Fragments[indx].Group_RB.fix_rotation();
    Nf_r -= 3;
  }
}


/*
double System::energy(){
//  int sz = active_interactions.size();
  int sz = interactions.size();
  double res = 0.0;
  int is_update = 0;
  int update_displ2;
  for(int i=0;i<sz;i++){
//    res += interactions[active_interactions[i]].calculate();
    res += interactions[i].calculate(update_displ2);
    if(update_displ2){ is_update = 1; }
  }

  if(is_update){
    // Update displacements:
    Boxold = Box;
    dT_2 = 0.0;
    for(i=0;i<Number_of_atoms;i++){ 
      Atoms[i].Atom_displ2 = 0.0;
      Atoms[i].Atom_RB_old.rb_cm = Atoms[i].Atom_RB.rb_cm;
    }
  }
  return res;
}

double System::energy(std::string s_int_type){
//  int sz = active_interactions.size();
  int sz = interactions.size();
  // The following conversion should be consistent with function:
  // void Interaction::set_interaction_type_and_functional(std::string t,std::string f) 
  int int_type = -1;
  if(s_int_type=="bond")         { int_type = 0; }
  else if(s_int_type=="angle")   { int_type = 1; }
  else if(s_int_type=="dihedral"){ int_type = 2; }
  else if(s_int_type=="oop")     { int_type = 3; }
  else if(s_int_type=="vdw")     { int_type = 4; }
  else if(s_int_type=="elec")    { int_type = 5; }
  else if(s_int_type=="mb")      { int_type = 6; }
  else if(s_int_type=="cg")      { int_type = 7; }

  double res = 0.0;
  int is_update = 0;
  int update_displ2;
  for(int i=0;i<sz;i++){
//    res += interactions[active_interactions[i]].calculate(int_type);
    res += interactions[i].calculate(int_type,update_displ2);
    if(update_displ2){ is_update = 1; }
  }

  if(is_update){
    // Update displacements:
    Boxold = Box;
    dT_2 = 0.0;
    for(i=0;i<Number_of_atoms;i++){ 
      Atoms[i].Atom_displ2 = 0.0;
      Atoms[i].Atom_RB_old.rb_cm = Atoms[i].Atom_RB.rb_cm;
    }
  }

  return res;
}

double System::energy_respa(std::string s_respa_type){
  double res,de; de = res = 0.0; 
  int respa_type = 0;
  int sz = interactions.size();

  if(s_respa_type=="fast"){  respa_type = 0; }
  else if(s_respa_type=="medium"){ respa_type = 1; }
  else if(s_respa_type=="slow") { respa_type = 2; }

  int update_displ2;
//  int count = 0;
  for(int i=0;i<sz;i++){
    if(interactions[i].get_respa_type()==respa_type){
      de = interactions[i].calculate(update_displ2);
      res += de;
      //cout<<"respa_type = "<<respa_type<<" type = "<<interactions[i].get_type()<<"  de = "<<de<<endl;
      //count++;
    }
  }
//  cout<<"count = "<<count<<endl;
  return res;
}

*/

double System::ekin_tr(){
  double res = 0.0;
  for(int i=0;i<Number_of_fragments;i++){
    res += Fragments[i].Group_RB.ekin_tr();
  }
  return res;
}

double System::ekin_tr_int(){
  double res = 0.0;
  VECTOR p_tot; p_tot = 0.0; // total momentum of the system
  double m_tot; m_tot = 0.0; // total mass of the system
  for(int i=0;i<Number_of_fragments;i++){
    res += Fragments[i].Group_RB.ekin_tr();
    p_tot +=  Fragments[i].Group_RB.rb_p; // sum up momenta
    m_tot +=  Fragments[i].Group_RB.rb_mass;
  }
  return ( res - (0.5/m_tot)*p_tot.length2() );
}


double System::ekin_rot(){
  double res = 0.0;
  for(int i=0;i<Number_of_fragments;i++){
    res += Fragments[i].Group_RB.ekin_rot();
  }
  return res;
}

double System::volume(){
  double res = 1e+10;
  if(is_Box){ res = Box.Determinant(); }
// else{ init_box(); res = Box.Determinant(); } // This is wrong!
  return res;
}


MATRIX3x3 System::pressure_tensor(){  

  MATRIX3x3 P_tens,Pid,Pvir,tmp;

  //------------ Kinetic energy part of pressure (ideal gas) -------------------
  Pid = 0.0;
  VECTOR p_tot; p_tot = 0.0;
  double m_tot; m_tot = 0.0;

  if(stress_opt=="at"){
    for(int i=0;i<Number_of_atoms;i++){
      RigidBody& top = Atoms[i].Atom_RB;  tmp.tensor_product(top.rb_p,top.rb_p);  Pid += top.rb_iM*tmp;
      p_tot += top.rb_p; m_tot += top.rb_mass;
    }// for i - all fragments
  }
  else if(stress_opt=="fr"){
    for(int i=0;i<Number_of_fragments;i++){
      RigidBody& top = Fragments[i].Group_RB; tmp.tensor_product(top.rb_p,top.rb_p); Pid += top.rb_iM*tmp;
      p_tot += top.rb_p; m_tot += top.rb_mass;
    }// for i - all fragments
  }
  else if(stress_opt=="ml"){
    for(int i=0;i<Number_of_molecules;i++){
      RigidBody& top = Molecules[i].Molecule_RB; tmp.tensor_product(top.rb_p,top.rb_p); Pid += top.rb_iM*tmp;
      p_tot += top.rb_p; m_tot += top.rb_mass;
    }// for i - all fragments
  }

  tmp.tensor_product(p_tot,p_tot); Pid -= (1.0/m_tot)*tmp;

  //----------------------------------------------------------------------------
  
  //------------ Virial part of the pessure (due to interactions )  ---------------
  int stat = 0; // status, if we need to recalculate stress
  if((stress_opt=="at")&&(is_stress_at==0)){ stat  = 1; }
  if((stress_opt=="fr")&&(is_stress_fr==0)){ stat  = 1; }
  if((stress_opt=="ml")&&(is_stress_ml==0)){ stat  = 1; }

  stat = 1;
  if(stat){

/*  Revise to get stresses from outside
  // In this case we recalculate (update) required stess tensor
  int nint = active_interactions.size();
  stress_at = 0.0;
  stress_fr = 0.0;
  stress_ml = 0.0;
  for(int i=0;i<nint;i++){
    if(interactions[active_interactions[i]].get_status()){
      if(stress_opt=="at"){ stress_at += interactions[active_interactions[i]].stress_at; } 
      else if(stress_opt=="fr"){ stress_fr += interactions[active_interactions[i]].stress_fr; }
      else if(stress_opt=="ml"){ stress_ml += interactions[active_interactions[i]].stress_ml; }
    }
  }
*/
  }// if stat

  if(stress_opt=="at"){ Pvir = stress_at; }
  else if(stress_opt=="fr"){ Pvir = stress_fr; }
  else if(stress_opt=="ml"){ Pvir = stress_ml; }

  //-------------------------------------------------------------------------------
  if(is_Box){  P_tens = (Pid + Pvir)/(volume());  }
  else{ P_tens = 0.0; }

  // Simmetrize pressure tensor:
  P_tens = 0.5*(P_tens + P_tens.T());

  return P_tens;

}


void System::init_fragment_velocities(double Temp){
  VECTOR TOT_P,TOT_L; TOT_P = 0.0; TOT_L = 0.0;
  init_fragment_velocities(Temp,TOT_P,TOT_L);
}

void System::init_fragment_velocities(double Temp,VECTOR TOT_P,VECTOR TOT_L){
  // This approach sets linear and angular momenta independently
  // and makes sure that the total linear and angular momenta are
  // proportional to corresponding arguments (or equal if this is consistent)
  // with the Temperature parameters
  // In any case the temperature is enforced to be constrained, while the
  // amplitude of the total momenta may be scaled if needed

//  VECTOR L_tot; L_tot = 0.0;
//  VECTOR P_tot; P_tot = 0.0;

  // Initialize random number generator
  srand( (unsigned)time(NULL) );
  int i;

  VECTOR* temp_p;
  VECTOR* temp_l;
  temp_p = new VECTOR[Number_of_fragments];
  temp_l = new VECTOR[Number_of_fragments];

  VECTOR tot_p; tot_p = 0.0;
  VECTOR tot_l; tot_l = 0.0;

  Random rnd;
  for(i=0;i<Number_of_fragments;i++){
    temp_p[i].x  = rnd.uniform(-0.5,0.5);
    temp_p[i].y  = rnd.uniform(-0.5,0.5);
    temp_p[i].z  = rnd.uniform(-0.5,0.5);
    temp_l[i].x  = rnd.uniform(-0.5,0.5);
    temp_l[i].y  = rnd.uniform(-0.5,0.5);
    temp_l[i].z  = rnd.uniform(-0.5,0.5);
    tot_p += temp_p[i];
  }

  double size = Number_of_fragments;
  tot_p = (TOT_P - tot_p)/size;
  for(i=0;i<Number_of_fragments;i++){   temp_p[i] +=  tot_p; }

  // Required temperature value scaling
  double temp_tr  = 0.0;
  double temp_rot = 0.0;
  for(i=0;i<Number_of_fragments;i++){
    RigidBody& top = Fragments[i].Group_RB;
    top.set_momentum(temp_p[i]);
    temp_tr += top.ekin_tr();
  }// for i

  // Rescale linear momenta to satisfy the translational kinetic energy
  double target_ekin_tr = 0.5*((double)Nf_t)*(boltzmann*Temp)/hartree;  // in a.u. of energy
  double scaling_factor_tr = (temp_tr==0.0)?0.0:sqrt(target_ekin_tr/temp_tr);
  for(i=0;i<Number_of_fragments;i++){  temp_p[i] *= scaling_factor_tr; }

  // Angular momenta
  for(i=0;i<Number_of_fragments;i++){
    RigidBody& top = Fragments[i].Group_RB;
    VECTOR tmp; tmp.cross(top.rb_cm,temp_p[i]);
    tot_l += tmp;
  }
  tot_l = (TOT_L - tot_l)/size;
  for(i=0;i<Number_of_fragments;i++){  temp_l[i] =  Fragments[i].Group_RB.rb_A_I_to_e * tot_l;   }

  temp_tr = 0.0;
  for(i=0;i<Number_of_fragments;i++){
    RigidBody& top = Fragments[i].Group_RB;
    top.set_angular_momentum(temp_l[i]);
    top.set_momentum(temp_p[i]);
    temp_rot += top.ekin_rot();
    temp_tr += top.ekin_tr();
  }

  double target_ekin = 0.5*((double)(Nf_r + Nf_t))*(boltzmann*Temp)/hartree;
  double scaling_factor = (temp_rot==0.0)?0.0:sqrt(target_ekin/(temp_rot+temp_tr));


  for(i=0;i<Number_of_fragments;i++){
      temp_p[i] *= scaling_factor;
      temp_l[i] *= scaling_factor;
  }

  // Set corresponding variables
  double E_kin = 0.0;
  VECTOR L_tot; L_tot = 0.0;
  VECTOR P_tot; P_tot = 0.0;
  for(i=0;i<Number_of_fragments;i++){
    RigidBody& top = Fragments[i].Group_RB; 
    top.set_momentum(temp_p[i]);
    top.set_angular_momentum(temp_l[i]);
    E_kin += (top.ekin_tr() + top.ekin_rot());

    VECTOR tmp; tmp.cross(top.rb_cm,top.rb_p);
    L_tot += top.rb_A_I_to_e_T * top.rb_l_e + tmp;
    P_tot += top.rb_p;

  }// for i - all fragments
  
  double curr_T =  2.0*(E_kin*hartree)/(((double)(Nf_t + Nf_r))*boltzmann);
  cout<<"in init_velocities...\n";
  cout<<"P_tot = "<<P_tot<<endl;
  cout<<"L_tot = "<<L_tot<<endl;
  cout<<"cutt_T = "<<curr_T<<endl;
  delete [] temp_p;
  delete [] temp_l;

}





}// namespace libchemsys
}// namespace libchemobjects
