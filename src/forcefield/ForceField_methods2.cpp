/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "ForceField.h"

/// liblibra namespace
namespace liblibra{

//using namespace liblinalg;
//using namespace libio;

namespace libforcefield{



void ForceField::set_functionals(boost::python::dict d){

  std::string key,value;
  for(int i=0;i<len(d.values());i++){    

    key = extract<std::string>(d.keys()[i]);
//    value = extract<std::string>(d.values()[i]);
//    std::cout<<"key = "<<key<<"  value = "<<value<<std::endl;

    if(key=="ForceField_Name"){ ForceField_Name = extract<std::string>(d.values()[i]); is_ForceField_Name = 1; }
    else if(key=="bond" || key=="bond_functional"){ bond_functional = extract<std::string>(d.values()[i]); is_bond_functional = 1;}
    else if(key=="angle" || key=="angle_functional"){ angle_functional = extract<std::string>(d.values()[i]); is_angle_functional = 1;}
    else if(key=="dihedral" || key=="dihedral_functional"){ dihedral_functional = extract<std::string>(d.values()[i]); is_dihedral_functional = 1;}
    else if(key=="oop" || key=="oop_functional"){ oop_functional = extract<std::string>(d.values()[i]); is_oop_functional = 1;}
    else if(key=="vdw" || key=="vdw_functional"){ vdw_functional = extract<std::string>(d.values()[i]); is_vdw_functional = 1;}
    else if(key=="elec" ||  key=="elec_functional"){ elec_functional = extract<std::string>(d.values()[i]); is_elec_functional = 1;}
    else if(key=="mb" || key=="mb_functional"){ mb_functional = extract<std::string>(d.values()[i]); is_mb_functional = 1;}
    else if(key=="cg" || key=="cg_functional"){ cg_functional = extract<std::string>(d.values()[i]); is_cg_functional = 1;}
    else if(key=="mb_excl" || key=="mb_excl_functional"){ mb_excl_functional = extract<std::string>(d.values()[i]); is_mb_excl_functional = 1;}

/*
    if(key=="ForceField_Name"){ ForceField_Name = value; is_ForceField_Name = 1; }
    else if(key=="bond" || key=="bond_functional"){ bond_functional = value; is_bond_functional = 1;}
    else if(key=="angle" || key=="angle_functional"){ angle_functional = value; is_angle_functional = 1;}
    else if(key=="dihedral" || key=="dihedral_functional"){ dihedral_functional = value; is_dihedral_functional = 1;}
    else if(key=="oop" || key=="oop_functional"){ oop_functional = value; is_oop_functional = 1;}
    else if(key=="vdw" || key=="vdw_functional"){ vdw_functional = value; is_vdw_functional = 1;}
    else if(key=="elec" || key=="angle_functional"){ elec_functional = value; is_elec_functional = 1;}
    else if(key=="mb" || key=="mb_functional"){ mb_functional = value; is_mb_functional = 1;}
    else if(key=="cg" || key=="cg_functional"){ cg_functional = value; is_cg_functional = 1;}
    else if(key=="mb_excl" || key=="mb_excl_functional"){ mb_excl_functional = value; is_mb_excl_functional = 1;}
*/
    else if(key=="system_pbc"){ system_pbc = extract<std::string>(d.values()[i]); is_system_pbc = 1; }
    else if(key=="pbc_degree_x"){ pbc_degree_x = extract<int>(d.values()[i]); is_pbc_degree_x = 1; }
    else if(key=="pbc_degree_y"){ pbc_degree_y = extract<int>(d.values()[i]); is_pbc_degree_y = 1; }
    else if(key=="pbc_degree_z"){ pbc_degree_z = extract<int>(d.values()[i]); is_pbc_degree_z = 1; }
    else if(key=="reciprocal_degree_x"){ reciprocal_degree_x = extract<int>(d.values()[i]); is_reciprocal_degree_x = 1; }
    else if(key=="reciprocal_degree_y"){ reciprocal_degree_y = extract<int>(d.values()[i]); is_reciprocal_degree_y = 1; }
    else if(key=="reciprocal_degree_z"){ reciprocal_degree_z = extract<int>(d.values()[i]); is_reciprocal_degree_z = 1; }
    else if(key=="R_vdw_on"){ R_vdw_on = extract<double>(d.values()[i]); is_R_vdw_on = 1; }
    else if(key=="R_vdw_off"){ R_vdw_off = extract<double>(d.values()[i]); is_R_vdw_off = 1; }
    else if(key=="R_elec_on"){ R_elec_on = extract<double>(d.values()[i]); is_R_elec_on = 1; }
    else if(key=="R_elec_off"){ R_elec_off = extract<double>(d.values()[i]); is_R_elec_off = 1; }
    else if(key=="R_vlist"){ R_vlist = extract<double>(d.values()[i]); is_R_vlist = 1; }
    else if(key=="elec_etha"){ elec_etha = extract<double>(d.values()[i]); is_elec_etha = 1; }

    else if(key=="sigma_comb_rule"){ sigma_comb_rule = extract<std::string>(d.values()[i]); is_sigma_comb_rule = 1; }
    else if(key=="epsilon_comb_rule"){ epsilon_comb_rule = extract<std::string>(d.values()[i]); is_epsilon_comb_rule = 1; }
    else if(key=="vdw_scale12"){ vdw_scale12 = extract<double>(d.values()[i]); is_vdw_scale12 = 1; }
    else if(key=="vdw_scale13"){ vdw_scale13 = extract<double>(d.values()[i]); is_vdw_scale13 = 1; }
    else if(key=="vdw_scale14"){ vdw_scale14 = extract<double>(d.values()[i]); is_vdw_scale14 = 1; }
    else if(key=="elec_scale12"){ elec_scale12 = extract<double>(d.values()[i]); is_elec_scale12 = 1; }
    else if(key=="elec_scale13"){ elec_scale13 = extract<double>(d.values()[i]); is_elec_scale13 = 1; }
    else if(key=="elec_scale14"){ elec_scale14 = extract<double>(d.values()[i]); is_elec_scale14 = 1; }


  }// for i

}


}// namespace libforcefield
}// namespace liblibra

