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

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libforcefield.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;


namespace libforcefield{




//---------------------------------------------------------------------
int (ForceField::*Add_Angle_Record1)(Angle_Record)       = &ForceField::Add_Angle_Record;
int (ForceField::*Add_Angle_Record2)(Angle_Record,int)   = &ForceField::Add_Angle_Record;

int (ForceField::*Add_Dihedral_Record1)(Dihedral_Record)       = &ForceField::Add_Dihedral_Record;
int (ForceField::*Add_Dihedral_Record2)(Dihedral_Record,int)   = &ForceField::Add_Dihedral_Record;

int (ForceField::*Add_Improper_Record1)(Dihedral_Record)       = &ForceField::Add_Improper_Record;



void export_forcefield_objects(){

    class_<Atom_Record>("Atom_Record",init<>())      
        .def(init<const Atom_Record&>())
        .def("__copy__", &generic__copy__<Atom_Record>)
        .def("__deepcopy__", &generic__deepcopy__<Atom_Record>)

        .def_readwrite("Atom_ff_type",&Atom_Record::Atom_ff_type)
        .def_readwrite("Atom_ff_int_type",&Atom_Record::Atom_ff_int_type)
        .def_readwrite("Atom_ff_type_H",&Atom_Record::Atom_ff_type_H)
        .def_readwrite("Atom_element",&Atom_Record::Atom_element)
        .def_readwrite("Atom_atomic_number",&Atom_Record::Atom_atomic_number)
        .def_readwrite("Atom_electronegativity",&Atom_Record::Atom_electronegativity)
        .def_readwrite("Atom_partial_charge",&Atom_Record::Atom_partial_charge)        
        .def_readwrite("Atom_ff_eq_int_type2",&Atom_Record::Atom_ff_eq_int_type2)
        .def_readwrite("Atom_ff_eq_int_type3",&Atom_Record::Atom_ff_eq_int_type3)
        .def_readwrite("Atom_ff_eq_int_type4",&Atom_Record::Atom_ff_eq_int_type4)
        .def_readwrite("Atom_ff_eq_int_type5",&Atom_Record::Atom_ff_eq_int_type5)
        .def_readwrite("Atom_radius",&Atom_Record::Atom_radius)
        .def_readwrite("Atom_Z_star",&Atom_Record::Atom_Z_star)
        .def_readwrite("Atom_theta",&Atom_Record::Atom_theta)
        .def_readwrite("Atom_sigma",&Atom_Record::Atom_sigma)
        .def_readwrite("Atom_epsilon",&Atom_Record::Atom_epsilon)
        .def_readwrite("Atom_GMP",&Atom_Record::Atom_GMP)
        .def_readwrite("Atom_crd",&Atom_Record::Atom_crd)
        .def_readwrite("Atom_val",&Atom_Record::Atom_val)
        .def_readwrite("Atom_pilp",&Atom_Record::Atom_pilp)
        .def_readwrite("Atom_mltb",&Atom_Record::Atom_mltb)
        .def_readwrite("Atom_arom",&Atom_Record::Atom_arom)
        .def_readwrite("Atom_lin",&Atom_Record::Atom_lin)
        .def_readwrite("Atom_sbmb",&Atom_Record::Atom_sbmb)
        .def_readwrite("Atom_alpha",&Atom_Record::Atom_alpha)
        .def_readwrite("Atom_N_eff",&Atom_Record::Atom_N_eff)
        .def_readwrite("Atom_A_scale",&Atom_Record::Atom_A_scale)
        .def_readwrite("Atom_G_scale",&Atom_Record::Atom_G_scale)
        .def_readwrite("Atom_DAN",&Atom_Record::Atom_DAN)
        .def_readwrite("Atom_pbci",&Atom_Record::Atom_pbci)
        .def_readwrite("Atom_fcadj",&Atom_Record::Atom_fcadj)
        .def_readwrite("Atom_dative",&Atom_Record::Atom_dative)
        .def_readwrite("Atom_brdr1",&Atom_Record::Atom_brdr1)
        .def_readwrite("Atom_brdr2",&Atom_Record::Atom_brdr2)
        .def_readwrite("Atom_brdr3",&Atom_Record::Atom_brdr3)
        .def_readwrite("Atom_such_n",&Atom_Record::Atom_such_n)
        .def_readwrite("Atom_such_m",&Atom_Record::Atom_such_m)
        .def_readwrite("Atom_such_a",&Atom_Record::Atom_such_a)
        .def_readwrite("Atom_such_D",&Atom_Record::Atom_such_D)
        .def_readwrite("Atom_such_c",&Atom_Record::Atom_such_c)
        .def("set",&Atom_Record::set)
        .def("show_info",&Atom_Record::show_info)
    ;

    class_< Atom_RecordList >("Atom_RecordList")
        .def(vector_indexing_suite< Atom_RecordList >())
    ;


    class_<Bond_Record>("Bond_Record",init<>())
        .def_readwrite("Atom1_ff_type",&Bond_Record::Atom1_ff_type)
        .def_readwrite("Atom1_ff_int_type",&Bond_Record::Atom1_ff_int_type)
        .def_readwrite("Atom2_ff_type",&Bond_Record::Atom2_ff_type)
        .def_readwrite("Atom2_ff_int_type",&Bond_Record::Atom2_ff_int_type)
        .def_readwrite("Atom1_atomic_number",&Bond_Record::Atom1_atomic_number)
        .def_readwrite("Atom2_atomic_number",&Bond_Record::Atom2_atomic_number)
        .def_readwrite("Atom1_element",&Bond_Record::Atom1_element)
        .def_readwrite("Atom2_element",&Bond_Record::Atom2_element)
        .def_readwrite("Bond_type_index",&Bond_Record::Bond_type_index)
        .def_readwrite("Bond_r_eq",&Bond_Record::Bond_r_eq)
        .def_readwrite("Bond_k_bond",&Bond_Record::Bond_k_bond)
        .def_readwrite("Bond_D_bond",&Bond_Record::Bond_D_bond)
        .def_readwrite("Bond_alpha",&Bond_Record::Bond_alpha)
        .def_readwrite("Bond_r_eq_ref",&Bond_Record::Bond_r_eq_ref)
        .def_readwrite("Bond_k_bond_ref",&Bond_Record::Bond_k_bond_ref)
        .def_readwrite("Bond_shift_elec",&Bond_Record::Bond_shift_elec)
        .def_readwrite("Bond_bci",&Bond_Record::Bond_bci)
        .def_readwrite("Bond_wij",&Bond_Record::Bond_wij)
        .def_readwrite("Bond_wij_1",&Bond_Record::Bond_wij_1)
        .def_readwrite("Bond_wij_2",&Bond_Record::Bond_wij_2)
        .def_readwrite("Bond_alpij",&Bond_Record::Bond_alpij)
        .def_readwrite("Bond_alpij_1",&Bond_Record::Bond_alpij_1)
        .def_readwrite("Bond_alpij_2",&Bond_Record::Bond_alpij_2)
        .def_readwrite("Bond_rij_1",&Bond_Record::Bond_rij_1)
        .def_readwrite("Bond_rij_2",&Bond_Record::Bond_rij_2)

        .def("set",&Bond_Record::set)
        .def("show_info",&Bond_Record::show_info)
    ;

    class_< Bond_RecordList >("Bond_RecordList")
        .def(vector_indexing_suite< Bond_RecordList >())
    ;

        
    class_<Angle_Record>("Angle_Record",init<>())
        .def_readwrite("Atom1_ff_type",&Angle_Record::Atom1_ff_type)
        .def_readwrite("Atom1_ff_int_type",&Angle_Record::Atom1_ff_int_type)
        .def_readwrite("Atom2_ff_type",&Angle_Record::Atom2_ff_type)
        .def_readwrite("Atom2_ff_int_type",&Angle_Record::Atom2_ff_int_type)
        .def_readwrite("Atom3_ff_type",&Angle_Record::Atom3_ff_type)
        .def_readwrite("Atom3_ff_int_type",&Angle_Record::Atom3_ff_int_type)
        .def_readwrite("Angle_type_index",&Angle_Record::Angle_type_index)
        .def_readwrite("Angle_theta_eq",&Angle_Record::Angle_theta_eq)
        .def_readwrite("Angle_k_angle",&Angle_Record::Angle_k_angle)
        .def_readwrite("Angle_r_eq",&Angle_Record::Angle_r_eq)
        .def_readwrite("Angle_k_ub",&Angle_Record::Angle_k_ub)
        .def("set",&Angle_Record::set)
        .def("show_info",&Angle_Record::show_info)
    ;

    class_< Angle_RecordList >("Angle_RecordList")
        .def(vector_indexing_suite< Angle_RecordList >())
    ;


    class_<Dihedral_Record>("Dihedral_Record",init<>())
        .def_readwrite("Atom1_ff_type",&Dihedral_Record::Atom1_ff_type)
        .def_readwrite("Atom1_ff_int_type",&Dihedral_Record::Atom1_ff_int_type)
        .def_readwrite("Atom2_ff_type",&Dihedral_Record::Atom2_ff_type)
        .def_readwrite("Atom2_ff_int_type",&Dihedral_Record::Atom2_ff_int_type)
        .def_readwrite("Atom3_ff_type",&Dihedral_Record::Atom3_ff_type)
        .def_readwrite("Atom3_ff_int_type",&Dihedral_Record::Atom3_ff_int_type)
        .def_readwrite("Atom4_ff_type",&Dihedral_Record::Atom4_ff_type)
        .def_readwrite("Atom4_ff_int_type",&Dihedral_Record::Atom4_ff_int_type)

        .def_readwrite("Dihedral_vphi",&Dihedral_Record::Dihedral_vphi)
        .def_readwrite("Dihedral_phase",&Dihedral_Record::Dihedral_phase)
        .def_readwrite("Dihedral_mult",&Dihedral_Record::Dihedral_mult)
        .def("set",&Dihedral_Record::set)
        .def("show_info",&Dihedral_Record::show_info)
    ;

    class_< Dihedral_RecordList >("Dihedral_RecordList")
        .def(vector_indexing_suite< Dihedral_RecordList >())
    ;


    class_<Fragment_Record>("Fragment_Record",init<>())
        .def_readwrite("Fragment_ff_type",&Fragment_Record::Fragment_ff_type)
        .def_readwrite("Fragment_ff_int_type",&Fragment_Record::Fragment_ff_int_type)
        .def_readwrite("Fragment_di",&Fragment_Record::Fragment_di)
        .def_readwrite("Fragment_li",&Fragment_Record::Fragment_li)
        .def_readwrite("Fragment_e0",&Fragment_Record::Fragment_e0)
        .def_readwrite("Fragment_rat",&Fragment_Record::Fragment_rat)
        .def_readwrite("Fragment_dw",&Fragment_Record::Fragment_dw)
        .def_readwrite("Fragment_mu",&Fragment_Record::Fragment_mu)
        .def_readwrite("Fragment_nu",&Fragment_Record::Fragment_nu)
        .def("set",&Fragment_Record::set)
        .def("show_info",&Fragment_Record::show_info)
    ;

    class_<ForceField>("ForceField",init<>())
        .def(init<boost::python::dict>())
        .def(init<const ForceField&>())
        .def("__copy__", &generic__copy__<ForceField>)
        .def("__deepcopy__", &generic__deepcopy__<ForceField>)

        .def_readwrite("ForceField_Name",&ForceField::ForceField_Name)
        .def_readwrite("bond_functional",&ForceField::bond_functional)
        .def_readwrite("angle_functional",&ForceField::angle_functional)
        .def_readwrite("dihedral_functional",&ForceField::dihedral_functional)
        .def_readwrite("oop_functional",&ForceField::oop_functional)
        .def_readwrite("vdw_functional",&ForceField::vdw_functional)
        .def_readwrite("elec_functional",&ForceField::elec_functional)
        .def_readwrite("mb_functional",&ForceField::mb_functional)
        .def_readwrite("cg_functional",&ForceField::cg_functional)
        .def_readwrite("mb_excl_functional",&ForceField::mb_excl_functional)
        .def_readwrite("sigma_comb_rule",&ForceField::sigma_comb_rule)
        .def_readwrite("epsilon_comb_rule",&ForceField::epsilon_comb_rule)
        .def_readwrite("vdw_scale12",&ForceField::vdw_scale12)
        .def_readwrite("vdw_scale13",&ForceField::vdw_scale13)
        .def_readwrite("vdw_scale14",&ForceField::vdw_scale14)
        .def_readwrite("elec_scale12",&ForceField::elec_scale12)
        .def_readwrite("elec_scale13",&ForceField::elec_scale13)
        .def_readwrite("elec_scale14",&ForceField::elec_scale14)

        .def_readwrite("Atom_Records",&ForceField::Atom_Records)
        .def_readwrite("Bond_Records",&ForceField::Bond_Records)
        .def_readwrite("Angle_Records",&ForceField::Angle_Records)
        .def_readwrite("Dihedral_Records",&ForceField::Dihedral_Records)

        .def("Add_Atom_Record",&ForceField::Add_Atom_Record)
        .def("Add_Bond_Record",&ForceField::Add_Bond_Record)
        .def("Add_Angle_Record",Add_Angle_Record1)
        .def("Add_Angle_Record",Add_Angle_Record2)
        .def("Add_Dihedral_Record",Add_Dihedral_Record1)
        .def("Add_Dihedral_Record",Add_Dihedral_Record2)
        .def("Add_Improper_Record",Add_Improper_Record1)
        .def("Add_Fragment_Record",&ForceField::Add_Fragment_Record)
        .def("show_atom_records",&ForceField::show_atom_records)
        .def("show_bond_records",&ForceField::show_bond_records)
        .def("show_angle_records",&ForceField::show_angle_records)
        .def("show_dihedral_records",&ForceField::show_dihedral_records)
        .def("show_improper_records",&ForceField::show_improper_records)
        .def("show_fragment_records",&ForceField::show_fragment_records)
        .def("set",&ForceField::set)
        .def("show_info",&ForceField::show_info)
        .def("set_functionals",&ForceField::set_functionals)

    ;


}// export_forcefield_objects


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygforcefield){
#else
BOOST_PYTHON_MODULE(libforcefield){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_forcefield_objects();

}




}// namespace libforcefield
}// namespace liblibra

