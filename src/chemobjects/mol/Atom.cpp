#include "Atom.h"

namespace libchemobjects{
namespace libmol{


void Atom::init_variables(){
  is_Atom_id      = 0;

  globAtom_Index  = 0;
  locAtom_Index   = 0;
  globGroup_Index = 0;
  globMolecule_Index = 0;

  is_Atom_RB      = 0;
  is_Atom_RB_old  = 0;
  Atom_displ2     = 0.0;   is_Atom_displ2 = 1;

  is_Atom_element = 0;
  is_Atom_atomic_radius = 0;
  is_Atom_charge  = 0;
  is_Atom_electronegativity = 0;

  is_Atom_formal_charge = 0;
  is_Atom_coordination = 0;
  is_Atom_functional_group = 0;
  is_Atom_ring_sizes = 0;
  Atom_min_ring_size = 1;  is_Atom_min_ring_size = 1;

  is_Atom_ff_type     = 0;

//  is_Atom_ff_int_type = 0;
//  Atom_is_surface_atom = 0;  is_Atom_is_surface_atom = 1;
//  is_Atom_surface_index = 0;
//  Atom_is_basis_atom = 0;    is_Atom_is_basis_atom = 1;
//  Atom_is_C60_CT = 0;        is_Atom_is_C60_CT = 1;

}

void Atom::copy_content(const Atom& at){

  universe = at.universe;

  // Basic information is copied without respect to if it has been defined
  // in destination object
  globAtom_Index     = at.globAtom_Index;
  locAtom_Index      = at.locAtom_Index;
  globGroup_Index    = at.globGroup_Index;
  globMolecule_Index = at.globMolecule_Index;

//  globAtom_Adjacent_Atoms = at.globAtom_Adjacent_Atoms;
  // Explicit copy
  if(at.globAtom_Adjacent_Atoms.size()){
    if(globAtom_Adjacent_Atoms.size()>0){ globAtom_Adjacent_Atoms.clear(); }
    for(int i=0;i<at.globAtom_Adjacent_Atoms.size();i++){ 
      globAtom_Adjacent_Atoms.push_back(at.globAtom_Adjacent_Atoms[i]);
    }
  }

  // This operator copies data from the source if and only if the corresponding data
  // is defined in the source
  if(at.is_Atom_id){ Atom_id = at.Atom_id; is_Atom_id = 1;}

  if(at.is_Atom_RB){ Atom_RB = at.Atom_RB; is_Atom_RB = 1; }
  if(at.is_Atom_RB_old){ Atom_RB_old = at.Atom_RB_old; is_Atom_RB_old = 1; }
  if(at.is_Atom_displ2){ Atom_displ2 = at.Atom_displ2; is_Atom_displ2 = 1; }

  if(at.is_Atom_element) { Atom_element = at.Atom_element;  is_Atom_element = 1;}
  if(at.is_Atom_atomic_radius){Atom_atomic_radius = at.Atom_atomic_radius;  is_Atom_atomic_radius = 1;} 
  if(at.is_Atom_charge){ Atom_charge = at.Atom_charge;  is_Atom_charge = 1;}
  if(at.is_Atom_electronegativity){Atom_electronegativity = at.Atom_electronegativity;  is_Atom_electronegativity = 1;}

  if(at.is_Atom_formal_charge){ Atom_formal_charge = at.Atom_formal_charge;  is_Atom_formal_charge = 1;}
  if(at.is_Atom_coordination){Atom_coordination = at.Atom_coordination;  is_Atom_coordination = 1;}
  if(at.is_Atom_functional_group){ Atom_functional_group = at.Atom_functional_group; is_Atom_functional_group = 1; }
  if(at.is_Atom_ring_sizes){Atom_ring_sizes = at.Atom_ring_sizes; is_Atom_ring_sizes = 1; }
  if(at.is_Atom_min_ring_size){Atom_min_ring_size = at.Atom_min_ring_size; is_Atom_min_ring_size = 1; }

  if(at.is_Atom_ff_type) { Atom_ff_type = at.Atom_ff_type;  is_Atom_ff_type = 1;}

//  if(at.is_Atom_ff_int_type) { Atom_ff_int_type = at.Atom_ff_int_type;  is_Atom_ff_int_type = 1;}
//  if(at.is_Atom_is_surface_atom){ Atom_is_surface_atom = at.Atom_is_surface_atom; is_Atom_is_surface_atom = 1;}
//  if(at.is_Atom_surface_index){ Atom_surface_index = at.Atom_surface_index; is_Atom_surface_index = 1;}
//  if(at.is_Atom_is_basis_atom){ Atom_is_basis_atom = at.Atom_is_basis_atom; is_Atom_is_basis_atom = 1;}
//  if(at.is_Atom_is_C60_CT){ Atom_is_C60_CT = at.Atom_is_C60_CT; is_Atom_is_C60_CT = 1; }

}

Atom::Atom(Universe& u){
  /****************
     Constructor
  ******************/
  universe = &u;
  // Initialize variables to default values
  init_variables();
}

Atom::Atom(Universe& u, boost::python::dict at){
  /****************
     Constructor
  ******************/
  universe = &u;
  // Initialize variables to default values
  init_variables();
  std::string key;
  for(int i=0;i<len(at.values());i++){
    key = extract<std::string>(at.keys()[i]);

    if(key=="Atom_id") { Atom_id = extract<int>(at.values()[i]);  is_Atom_id = 1; }
    else if(key=="Atom_element") { 
      Atom_element = extract<std::string>(at.values()[i]); is_Atom_element = 1; 
      Element elt = universe->Get_Element(Atom_element);
      // Get basic atomic properties
      if(elt.is_Elt_mass){  Atom_RB.set_mass(elt.Elt_mass); }
    }
    else if(key=="Atom_atomic_radius") { Atom_atomic_radius = extract<double>(at.values()[i]); is_Atom_atomic_radius = 1; }
    else if(key=="Atom_charge") { Atom_charge = extract<double>(at.values()[i]); is_Atom_charge = 1; }
    else if(key=="Atom_electronegativity") { Atom_electronegativity = extract<double>(at.values()[i]); is_Atom_electronegativity = 1; }
    else if(key=="Atom_formal_charge") { Atom_formal_charge = extract<double>(at.values()[i]); is_Atom_formal_charge = 1; }
    else if(key=="Atom_ff_type") { Atom_ff_type = extract<std::string>(at.values()[i]); is_Atom_ff_type = 1; }

    else if(key=="Atom_x") {
      Atom_RB.rb_cm.x = extract<double>(at.values()[i]); is_Atom_RB = 1; Atom_RB.is_rb_cm = 1;
      Atom_RB_old.rb_cm.x = Atom_RB.rb_cm.x; is_Atom_RB_old = 1; Atom_RB.is_rb_cm = 1;
    }
    else if(key=="Atom_y") {
      Atom_RB.rb_cm.y = extract<double>(at.values()[i]); is_Atom_RB = 1; Atom_RB.is_rb_cm = 1;
      Atom_RB_old.rb_cm.y = Atom_RB.rb_cm.y; is_Atom_RB_old = 1; Atom_RB.is_rb_cm = 1;
    }
    else if(key=="Atom_z") { 
      Atom_RB.rb_cm.z = extract<double>(at.values()[i]); is_Atom_RB = 1; Atom_RB.is_rb_cm = 1;
      Atom_RB_old.rb_cm.z = Atom_RB.rb_cm.z; is_Atom_RB_old = 1; Atom_RB.is_rb_cm = 1;
    }


  }// for i
}



Atom::Atom(const Atom& at){
  /********************
    Copy constructor
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of at object which is defined
  copy_content(at);
}

Atom::~Atom(){
  if(globAtom_Adjacent_Atoms.size()>0) { globAtom_Adjacent_Atoms.clear(); }
}

Atom& Atom::operator=(const Atom& at){
  /********************
    Assignment operator
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of at object which is defined
  copy_content(at);
  return *this;
}

void Atom::set(object at){
  set_value(is_Atom_id,      Atom_id,      at,"Atom_id");

  set_value(is_Atom_element, Atom_element, at,"Atom_element");
  set_value(is_Atom_atomic_radius,        Atom_atomic_radius,       at,"Atom_atomic_radius");
  set_value(is_Atom_charge,  Atom_charge,  at,"Atom_charge");
  set_value(is_Atom_electronegativity,    Atom_electronegativity,    at,"Atom_electronegativity");

  set_value(is_Atom_formal_charge,  Atom_formal_charge,  at,"Atom_formal_charge");
  set_value(is_Atom_coordination,   Atom_coordination,   at,"Atom_coordination");

  set_value(is_Atom_ff_type, Atom_ff_type, at, "Atom_ff_type");

//  set_value(is_Atom_ff_int_type, Atom_ff_int_type, at, "Atom_ff_int_type");
//  set_value(is_Atom_is_surface_atom, Atom_is_surface_atom, at, "Atom_is_surface_atom");
//  set_value(is_Atom_surface_index, Atom_surface_index, at, "Atom_surface_index");
//  set_value(is_Atom_is_basis_atom,  Atom_is_basis_atom, at, "Atom_is_basis_atom");
//  set_value(is_Atom_is_C60_CT, Atom_is_C60_CT, at, "Atom_is_C60_CT");

}

void Atom::show_info(){

//   std::cout<<"Atom "<<globAtom_Index<<" topology:"<<std::endl;
//   std::cout<<"globAtom_Index = "<<globAtom_Index<<std::endl;
//   std::cout<<"locAtom_Index = "<<locAtom_Index<<std::endl;
//   std::cout<<"globGroup_Index = "<<globGroup_Index<<std::endl;
//   std::cout<<"globMolecule_Index = "<<globMolecule_Index<<std::endl;
 
//   std::cout<<"Atom "<<globAtom_Index<<" properties:"<<std::endl;
   if(is_Atom_id)     {std::cout<<"Atom_id = "<<Atom_id<<std::endl;     }
   if(is_Atom_displ2) {std::cout<<"Atom_displ2 = "<<Atom_displ2<<std::endl; }
//   if(is_Atom_RB){  Atom_RB.show_info(); }

   if(is_Atom_element){std::cout<<"Atom_element = "<<Atom_element<<std::endl;}
   if(is_Atom_atomic_radius) {std::cout<<"Atom_atomic_radius = "<<Atom_atomic_radius<<std::endl;}
   if(is_Atom_charge) {std::cout<<"Atom_charge = "<<Atom_charge<<std::endl;}
   if(is_Atom_electronegativity) {std::cout<<"Atom_electronegativity = "<<Atom_electronegativity<<std::endl;}

   if(is_Atom_formal_charge) {std::cout<<"Atom_formal_charge = "<<Atom_formal_charge<<std::endl;}
   if(is_Atom_coordination) {std::cout<<"Atom_coordination = "<<Atom_coordination<<std::endl;}
   if(is_Atom_functional_group) {std::cout<<"Atom_functional_group = "<<Atom_functional_group<<std::endl;}
   if(is_Atom_ring_sizes){ std::cout<<"Atom_ring_sizes = "; 
       for(int i=0;i<Atom_ring_sizes.size();i++){ std::cout<<Atom_ring_sizes[i]<<"  "; }
       std::cout<<std::endl;
   }
   if(is_Atom_min_ring_size){ std::cout<<"Atom_min_ring_size = "<<Atom_min_ring_size<<std::endl;}

   if(is_Atom_ff_type) {std::cout<<"Atom_ff_type = "<<Atom_ff_type<<std::endl;}

//   if(is_Atom_ff_int_type) {std::cout<<"Atom_ff_int_type = "<<Atom_ff_int_type<<std::endl;}
//   if(is_Atom_is_surface_atom){std::cout<<"Atom_is_surface_atom = "<<Atom_is_surface_atom<<std::endl;}
//   if(is_Atom_surface_index){std::cout<<"Atom_surface_index = "<<Atom_surface_index<<std::endl;}
//   if(is_Atom_is_basis_atom){std::cout<<"Atom_is_basis_atom = "<<Atom_is_basis_atom<<std::endl;}
//   if(is_Atom_is_C60_CT){ std::cout<<"Atom_is_C60_CT = "<<Atom_is_C60_CT<<std::endl; }
   std::cout<<std::endl;
}

void Atom::save(boost::property_tree::ptree& pt,std::string path){

  if(is_Atom_id){  ::save(pt,path+".Atom_id",Atom_id);    }

  ::save(pt,path+".globAtom_Index",globAtom_Index);
  ::save(pt,path+".locAtom_Index",locAtom_Index);
  ::save(pt,path+".globGroup_Index",globGroup_Index);
  ::save(pt,path+".globMolecule_Index",globMolecule_Index);
  ::save(pt,path+".globAtom_Adjacent_Atoms",globAtom_Adjacent_Atoms);

  if(is_Atom_RB){  Atom_RB.save(pt,path+".Atom_RB");    }
  if(is_Atom_RB_old){  Atom_RB_old.save(pt,path+".Atom_RB_old");    }
  //if(is_Atom_RB){  ::save(pt,path+".Atom_RB",Atom_RB);    }
  //if(is_Atom_RB_old){  ::save(pt,path+".Atom_RB_old",Atom_RB_old);    }

  if(is_Atom_displ2) { ::save(pt,path+".Atom_displ2",Atom_displ2); }

  if(is_Atom_element){  ::save(pt,path+".Atom_element",Atom_element);    }
  if(is_Atom_atomic_radius){  ::save(pt,path+".Atom_atomic_radius",Atom_atomic_radius);    }
  if(is_Atom_charge){  ::save(pt,path+".Atom_charge",Atom_charge);    }
  if(is_Atom_electronegativity){  ::save(pt,path+".Atom_electronegativity",Atom_electronegativity);    }
  if(is_Atom_formal_charge){  ::save(pt,path+".Atom_formal_charge",Atom_formal_charge);    }
  if(is_Atom_coordination){  ::save(pt,path+".Atom_coordination",Atom_coordination);    }
  if(is_Atom_functional_group){  ::save(pt,path+".Atom_functional_group",Atom_functional_group);    }
  if(is_Atom_ring_sizes){  ::save(pt,path+".Atom_ring_sizes",Atom_ring_sizes);    }
  if(is_Atom_min_ring_size){  ::save(pt,path+".Atom_min_ring_size",Atom_min_ring_size);    }
  if(is_Atom_ff_type){  ::save(pt,path+".Atom_ff_type",Atom_ff_type);    }

}

void save(boost::property_tree::ptree& pt,std::string path,vector<Atom>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".Atom"+rt);
  }
}


void Atom::load(boost::property_tree::ptree& pt,std::string path,int& status){

  int st;
  status = 0;

  ::load(pt,path+".Atom_id",Atom_id,is_Atom_id); if(is_Atom_id==1) { status=1;}

  ::load(pt,path+".globAtom_Index",globAtom_Index,st); if(st==1) { status=1;}
  ::load(pt,path+".locAtom_Index",locAtom_Index,st); if(st==1) { status=1;}
  ::load(pt,path+".globGroup_Index",globGroup_Index,st); if(st==1) { status=1;}
  ::load(pt,path+".globMolecule_Index",globMolecule_Index,st); if(st==1) { status=1;}
  ::load(pt,path+".globAtom_Adjacent_Atoms",globAtom_Adjacent_Atoms,st); if(st==1) { status=1;}

  Atom_RB.load(pt,path+".Atom_RB",is_Atom_RB); if(is_Atom_RB==1) { status=1;}
  Atom_RB_old.load(pt,path+".Atom_RB_old",is_Atom_RB_old); if(is_Atom_RB_old==1) { status=1;}
//  ::load(pt,path+".Atom_RB",Atom_RB, is_Atom_RB); if(is_Atom_RB==1) { status=1;}
//  ::load(pt,path+".Atom_RB_old",Atom_RB_old,is_Atom_RB_old); if(is_Atom_RB_old==1) { status=1;}

  ::load(pt,path+".Atom_displ2",Atom_displ2,is_Atom_displ2); if(is_Atom_displ2==1) { status=1;}

  ::load(pt,path+".Atom_element",Atom_element,is_Atom_element); if(is_Atom_element==1) { status=1;}
  ::load(pt,path+".Atom_atomic_radius",Atom_atomic_radius,is_Atom_atomic_radius); if(is_Atom_atomic_radius==1) { status=1;}
  ::load(pt,path+".Atom_charge",Atom_charge,is_Atom_charge); if(is_Atom_charge==1) { status=1;}
  ::load(pt,path+".Atom_electronegativity",Atom_electronegativity,is_Atom_electronegativity); if(is_Atom_electronegativity==1) { status=1;}
  ::load(pt,path+".Atom_formal_charge",Atom_formal_charge,is_Atom_formal_charge); if(is_Atom_formal_charge==1) { status=1;}
  ::load(pt,path+".Atom_coordination",Atom_coordination,is_Atom_coordination); if(is_Atom_coordination==1) { status=1;}
  ::load(pt,path+".Atom_functional_group",Atom_functional_group,is_Atom_functional_group); if(is_Atom_functional_group==1) { status=1;}
  ::load(pt,path+".Atom_ring_sizes",Atom_ring_sizes,is_Atom_ring_sizes); if(is_Atom_ring_sizes==1) { status=1;}
  ::load(pt,path+".Atom_min_ring_size",Atom_min_ring_size,is_Atom_min_ring_size); if(is_Atom_min_ring_size==1) { status=1;}
  ::load(pt,path+".Atom_ff_type",Atom_ff_type,is_Atom_ff_type); if(is_Atom_ff_type==1) { status=1;}

}

void load(boost::property_tree::ptree& pt,std::string path,vector<Atom>& vt,Universe& u,int& status){
  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      Atom x(u); x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}


}// namespace libmol
}// namespace libchemobjects


