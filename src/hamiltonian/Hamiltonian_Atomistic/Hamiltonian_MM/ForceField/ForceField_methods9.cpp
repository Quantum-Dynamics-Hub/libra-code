#include "ForceField.h"

namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_mm{
namespace libforcefield{




int ForceField::set_ff_epsilon_and_sigma(int sz, vector<string> types,double** epsilon, double** sigma){

  int res = 1;
  for(int i=0;i<sz;i++){
    //-------------- Start with looking the whole record ---------------------
    // Find index of Atom_Record corresponding to the force field type types[i]
    int ff_atom1_indx = Atom_Record_Index(types[i]);
    if(ff_atom1_indx==-1){
      *epsilon[i] = 0.0; *sigma[i] = 3.0 * Angst;
      cout<<"Warning: In ForceField::set_ff_epsilon_and_sigma : Can not fing the atom type ff_atom type in force field "<<ForceField_Name<<"\n";
      cout<<"Setting epsilon to zero and sigma to 3.0\n";
    }
    else{
      map<string,double> d;
      get_vdw_parameters(types[i],types[i],"none",d);
      *epsilon[i] = d["epsilon"];
      *sigma[i] = d["sigma"];
    }
  }// for i

  return res;

}


}// namespace libforcefield
}// namespace libhamiltonian_mm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian

