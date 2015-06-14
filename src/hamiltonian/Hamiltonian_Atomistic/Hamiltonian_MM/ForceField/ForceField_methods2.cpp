#include "../ForceField.h"

void ForceField::set_functionals(boost::python::dict d){
  std::string key,value;
  for(int i=0;i<len(d.values());i++){    
    key = extract<std::string>(d.keys()[i]);
    value = extract<std::string>(d.values()[i]);
    std::cout<<"key = "<<key<<"  value = "<<value<<std::endl;
    if(key=="bond"){ bond_functional = value; is_bond_functional = 1;}
    else if(key=="angle"){ angle_functional = value; is_angle_functional = 1;}
    else if(key=="dihedral"){ dihedral_functional = value; is_dihedral_functional = 1;}
    else if(key=="oop"){ oop_functional = value; is_oop_functional = 1;}
    else if(key=="vdw"){ vdw_functional = value; is_vdw_functional = 1;}
    else if(key=="elec"){ elec_functional = value; is_elec_functional = 1;}
  }// for i
}

