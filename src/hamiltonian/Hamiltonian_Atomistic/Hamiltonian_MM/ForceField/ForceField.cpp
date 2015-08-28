/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "ForceField.h"

//#include "../../../../mmath/libmmath.h"
using namespace libmmath;

//#include "../../../../io/libio.h"
using namespace libio;



namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_mm{
namespace libforcefield{



//--------------- Atom Record class members ---------------------------

int Atom_Record::set(object at){

   set_value(is_Atom_ff_type,    Atom_ff_type,    at,"Atom_ff_type");
   set_value(is_Atom_ff_type_H,  Atom_ff_type_H,  at,"Atom_ff_type_H");
   set_value(is_Atom_ff_int_type,Atom_ff_int_type,at,"Atom_ff_int_type");
   set_value(is_Atom_element,    Atom_element,    at,"Atom_element");
   set_value(is_Atom_atomic_number, Atom_atomic_number, at,"Atom_atomic_number");

   set_value(is_Atom_electronegativity, Atom_electronegativity, at,"Atom_electronegativity");
   set_value(is_Atom_partial_charge, Atom_partial_charge, at,"Atom_partial_charge");

   set_value(is_Atom_ff_eq_int_type2,Atom_ff_eq_int_type2, at, "Atom_ff_eq_int_type2");
   set_value(is_Atom_ff_eq_int_type3,Atom_ff_eq_int_type3, at, "Atom_ff_eq_int_type3");
   set_value(is_Atom_ff_eq_int_type4,Atom_ff_eq_int_type4, at, "Atom_ff_eq_int_type4");
   set_value(is_Atom_ff_eq_int_type5,Atom_ff_eq_int_type5, at, "Atom_ff_eq_int_type5");
 
   set_value(is_Atom_radius,Atom_radius,at,"Atom_radius");
   set_value(is_Atom_Z_star,Atom_Z_star,at,"Atom_Z_star");
   set_value(is_Atom_theta,Atom_theta,at,"Atom_theta");
   set_value(is_Atom_sigma,Atom_sigma,at,"Atom_sigma");
   set_value(is_Atom_epsilon,Atom_epsilon,at,"Atom_epsilon");
   set_value(is_Atom_GMP,Atom_GMP,at,"Atom_GMP");

   set_value(is_Atom_crd,Atom_crd,at,"Atom_crd");
   set_value(is_Atom_val,Atom_val,at,"Atom_val");
   set_value(is_Atom_pilp,Atom_pilp,at,"Atom_pilp");
   set_value(is_Atom_mltb,Atom_mltb,at,"Atom_mltb");
   set_value(is_Atom_arom,Atom_arom,at,"Atom_arom");
   set_value(is_Atom_lin,Atom_lin,at,"Atom_lin");
   set_value(is_Atom_sbmb,Atom_sbmb,at,"Atom_sbmb");
   set_value(is_Atom_alpha,Atom_alpha,at,"Atom_alpha");
   set_value(is_Atom_N_eff,Atom_N_eff,at,"Atom_N_eff");
   set_value(is_Atom_A_scale,Atom_A_scale,at,"Atom_A_scale");
   set_value(is_Atom_G_scale,Atom_G_scale,at,"Atom_G_scale");
   set_value(is_Atom_DAN,Atom_DAN,at,"Atom_DAN");
   set_value(is_Atom_pbci,Atom_pbci,at,"Atom_pbci");
   set_value(is_Atom_fcadj,Atom_fcadj,at,"Atom_fcadj");


   set_value(is_Atom_dative,Atom_dative,at,"Atom_dative");
   set_value(is_Atom_brdr1, Atom_brdr1, at,"Atom_brdr1");
   set_value(is_Atom_brdr2, Atom_brdr2, at,"Atom_brdr2");
   set_value(is_Atom_brdr3, Atom_brdr3, at,"Atom_brdr3");

   set_value(is_Atom_such_n,Atom_such_n,at,"Atom_such_n");
   set_value(is_Atom_such_m,Atom_such_m,at,"Atom_such_m");
   set_value(is_Atom_such_a,Atom_such_a,at,"Atom_such_a");
   set_value(is_Atom_such_D,Atom_such_D,at,"Atom_such_D");
   set_value(is_Atom_such_c,Atom_such_c,at,"Atom_such_c");

   return 1;
}



int Atom_Record::show_info(){

   std::cout<<"Atom_Record properties:"<<std::endl;
   if(is_Atom_ff_type)        {std::cout<<"Atom_ff_type = "<<Atom_ff_type<<std::endl;   }
   if(is_Atom_ff_type_H)      {std::cout<<"Atom_ff_type_H = "<<Atom_ff_type_H<<std::endl;   }
   if(is_Atom_ff_int_type)    {std::cout<<"Atom_ff_int_type = "<<Atom_ff_int_type<<std::endl;     }
   if(is_Atom_element)        {std::cout<<"Atom_element = "<<Atom_element<<std::endl;   }
   if(is_Atom_atomic_number)  {std::cout<<"Atom_atomic_number = "<<Atom_atomic_number<<std::endl;   }

   if(is_Atom_electronegativity)  {std::cout<<"Atom_electronegativity = "<<Atom_electronegativity<<std::endl;   }
   if(is_Atom_partial_charge)  {std::cout<<"Atom_partial_charge = "<<Atom_partial_charge<<std::endl;   }

   if(is_Atom_ff_eq_int_type2){std::cout<<"Atom_ff_eq_int_type2 = "<<Atom_ff_eq_int_type2<<std::endl; }
   if(is_Atom_ff_eq_int_type3){std::cout<<"Atom_ff_eq_int_type3 = "<<Atom_ff_eq_int_type3<<std::endl; }
   if(is_Atom_ff_eq_int_type4){std::cout<<"Atom_ff_eq_int_type4 = "<<Atom_ff_eq_int_type4<<std::endl; }
   if(is_Atom_ff_eq_int_type5){std::cout<<"Atom_ff_eq_int_type5 = "<<Atom_ff_eq_int_type5<<std::endl; }

   if(is_Atom_radius)         {std::cout<<"Atom_radius = "<<Atom_radius<<std::endl;   }
   if(is_Atom_Z_star)         {std::cout<<"Atom_Z_star = "<<Atom_Z_star<<std::endl;   }
   if(is_Atom_theta)          {std::cout<<"Atom_theta = "<<Atom_theta<<std::endl;   }
   if(is_Atom_sigma)          {std::cout<<"Atom_sigma = "<<Atom_sigma<<std::endl;   }
   if(is_Atom_epsilon)        {std::cout<<"Atom_epsilon = "<<Atom_epsilon<<std::endl;   }
   if(is_Atom_GMP)            {std::cout<<"Atom_GMP = "<<Atom_GMP<<std::endl; }

   if(is_Atom_crd)            {std::cout<<"Atom_crd = "<<Atom_crd<<std::endl; }
   if(is_Atom_val)            {std::cout<<"Atom_val = "<<Atom_val<<std::endl; }
   if(is_Atom_pilp)           {std::cout<<"Atom_pilp = "<<Atom_pilp<<std::endl; }
   if(is_Atom_mltb)           {std::cout<<"Atom_mltb = "<<Atom_mltb<<std::endl; }
   if(is_Atom_arom)           {std::cout<<"Atom_arom = "<<Atom_arom<<std::endl; }
   if(is_Atom_lin)            {std::cout<<"Atom_lin = "<<Atom_lin<<std::endl; }
   if(is_Atom_sbmb)           {std::cout<<"Atom_sbmb = "<<Atom_sbmb<<std::endl; }
   if(is_Atom_alpha)          {std::cout<<"Atom_alpha = "<<Atom_alpha<<std::endl; }
   if(is_Atom_N_eff)          {std::cout<<"Atom_N_eff = "<<Atom_N_eff<<std::endl; }
   if(is_Atom_A_scale)        {std::cout<<"Atom_A_scale = "<<Atom_A_scale<<std::endl; }
   if(is_Atom_G_scale)        {std::cout<<"Atom_G_scale = "<<Atom_G_scale<<std::endl; }
   if(is_Atom_DAN)            {std::cout<<"Atom_DAN = "<<Atom_DAN<<std::endl; }
   if(is_Atom_pbci)           {std::cout<<"Atom_pbci = "<<Atom_pbci<<std::endl; }
   if(is_Atom_fcadj)          {std::cout<<"Atom_fcadj = "<<Atom_fcadj<<std::endl; }

   if(is_Atom_dative)         {std::cout<<"Atom_dative = "<<Atom_dative<<std::endl; }
   if(is_Atom_brdr1)          {std::cout<<"Atom_brdr1  = "<<Atom_brdr1<<std::endl;  }
   if(is_Atom_brdr2)          {std::cout<<"Atom_brdr2  = "<<Atom_brdr2<<std::endl;  }
   if(is_Atom_brdr3)          {std::cout<<"Atom_brdr3  = "<<Atom_brdr3<<std::endl;  }

   if(is_Atom_such_n)         {std::cout<<"Atom_such_n = "<<Atom_such_n<<std::endl; }
   if(is_Atom_such_m)         {std::cout<<"Atom_such_m = "<<Atom_such_m<<std::endl; }
   if(is_Atom_such_a)         {std::cout<<"Atom_such_a = "<<Atom_such_a<<std::endl; }
   if(is_Atom_such_D)         {std::cout<<"Atom_such_D = "<<Atom_such_D<<std::endl; }
   if(is_Atom_such_c)         {std::cout<<"Atom_such_c = "<<Atom_such_c<<std::endl; }

   std::cout<<std::endl;   

   return 1;

}

Atom_Record& Atom_Record::operator=(const Atom_Record& at){

   is_Atom_ff_type     = 0;
   is_Atom_ff_type_H   = 0;
   is_Atom_ff_int_type = 0;
   is_Atom_element     = 0;
   is_Atom_atomic_number = 0;

   is_Atom_electronegativity = 0;
   is_Atom_partial_charge = 0;

   is_Atom_ff_eq_int_type2 = 0;
   is_Atom_ff_eq_int_type3 = 0;
   is_Atom_ff_eq_int_type4 = 0;
   is_Atom_ff_eq_int_type5 = 0;

   is_Atom_radius      = 0;
   is_Atom_Z_star      = 0;
   is_Atom_theta       = 0;
   is_Atom_sigma       = 0;
   is_Atom_epsilon     = 0;
   is_Atom_GMP         = 0;

   is_Atom_crd         = 0;
   is_Atom_val         = 0;
   is_Atom_pilp        = 0;
   is_Atom_mltb        = 0;
   is_Atom_arom        = 0;
   is_Atom_lin         = 0;
   is_Atom_sbmb        = 0;
   is_Atom_alpha       = 0;
   is_Atom_N_eff       = 0;
   is_Atom_A_scale     = 0;
   is_Atom_G_scale     = 0;
   is_Atom_DAN         = 0;
   is_Atom_pbci        = 0;
   is_Atom_fcadj       = 0;

   is_Atom_dative      = 0;
   is_Atom_brdr1       = 0;
   is_Atom_brdr2       = 0;
   is_Atom_brdr3       = 0;

   is_Atom_such_n      = 0;
   is_Atom_such_m      = 0;
   is_Atom_such_a      = 0;
   is_Atom_such_D      = 0;
   is_Atom_such_c      = 0;



// This assignment is only working if the source object contain some data

   if(at.is_Atom_ff_type)     {Atom_ff_type     = at.Atom_ff_type;  is_Atom_ff_type = 1;}
   if(at.is_Atom_ff_type_H)   {Atom_ff_type_H   = at.Atom_ff_type_H;is_Atom_ff_type_H = 1;}
   if(at.is_Atom_ff_int_type) {Atom_ff_int_type = at.Atom_ff_int_type; is_Atom_ff_int_type = 1;}
   if(at.is_Atom_element)     {Atom_element     = at.Atom_element;  is_Atom_element = 1;}
   if(at.is_Atom_atomic_number){Atom_atomic_number = at.Atom_atomic_number;  is_Atom_atomic_number = 1;}

   if(at.is_Atom_electronegativity){Atom_electronegativity = at.Atom_electronegativity;  is_Atom_electronegativity = 1;}
   if(at.is_Atom_partial_charge){Atom_partial_charge = at.Atom_partial_charge;  is_Atom_partial_charge = 1;}

   if(at.is_Atom_ff_eq_int_type2){ Atom_ff_eq_int_type2 = at.Atom_ff_eq_int_type2; is_Atom_ff_eq_int_type2 = 1; }
   if(at.is_Atom_ff_eq_int_type3){ Atom_ff_eq_int_type3 = at.Atom_ff_eq_int_type3; is_Atom_ff_eq_int_type3 = 1; }
   if(at.is_Atom_ff_eq_int_type4){ Atom_ff_eq_int_type4 = at.Atom_ff_eq_int_type4; is_Atom_ff_eq_int_type4 = 1; }
   if(at.is_Atom_ff_eq_int_type5){ Atom_ff_eq_int_type5 = at.Atom_ff_eq_int_type5; is_Atom_ff_eq_int_type5 = 1; }


   if(at.is_Atom_radius)      {Atom_radius      = at.Atom_radius;   is_Atom_radius = 1;}
   if(at.is_Atom_Z_star)      {Atom_Z_star      = at.Atom_Z_star;   is_Atom_Z_star = 1;}
   if(at.is_Atom_theta)       {Atom_theta       = at.Atom_theta;    is_Atom_theta  = 1;}
   if(at.is_Atom_sigma)       {Atom_sigma       = at.Atom_sigma;    is_Atom_sigma  = 1;}
   if(at.is_Atom_epsilon)     {Atom_epsilon     = at.Atom_epsilon;  is_Atom_epsilon = 1;}
   if(at.is_Atom_GMP)         {Atom_GMP         = at.Atom_GMP;      is_Atom_GMP    = 1;}

   if(at.is_Atom_crd)         {Atom_crd         = at.Atom_crd;      is_Atom_crd    = 1;}
   if(at.is_Atom_val)         {Atom_val         = at.Atom_val;      is_Atom_val    = 1;}
   if(at.is_Atom_pilp)        {Atom_pilp        = at.Atom_pilp;     is_Atom_pilp   = 1;}
   if(at.is_Atom_mltb)        {Atom_mltb        = at.Atom_mltb;     is_Atom_mltb   = 1;}
   if(at.is_Atom_arom)        {Atom_arom        = at.Atom_arom;     is_Atom_arom   = 1;}
   if(at.is_Atom_lin)         {Atom_lin         = at.Atom_lin;      is_Atom_lin    = 1;}
   if(at.is_Atom_sbmb)        {Atom_sbmb        = at.Atom_sbmb;     is_Atom_sbmb   = 1;}
   if(at.is_Atom_alpha)       {Atom_alpha       = at.Atom_alpha;    is_Atom_alpha  = 1;}
   if(at.is_Atom_N_eff)       {Atom_N_eff       = at.Atom_N_eff;    is_Atom_N_eff  = 1;}
   if(at.is_Atom_A_scale)     {Atom_A_scale     = at.Atom_A_scale;  is_Atom_A_scale= 1;}
   if(at.is_Atom_G_scale)     {Atom_G_scale     = at.Atom_G_scale;  is_Atom_G_scale= 1;}
   if(at.is_Atom_DAN)         {Atom_DAN         = at.Atom_DAN;      is_Atom_DAN= 1;}
   if(at.is_Atom_pbci)        {Atom_pbci        = at.Atom_pbci;     is_Atom_pbci= 1;}
   if(at.is_Atom_fcadj)       {Atom_fcadj       = at.Atom_fcadj;    is_Atom_fcadj= 1;}

   if(at.is_Atom_dative)      {Atom_dative      = at.Atom_dative;   is_Atom_dative = 1;}
   if(at.is_Atom_brdr1)       {Atom_brdr1       = at.Atom_brdr1;    is_Atom_brdr1  = 1;}
   if(at.is_Atom_brdr2)       {Atom_brdr2       = at.Atom_brdr2;    is_Atom_brdr2  = 1;}
   if(at.is_Atom_brdr3)       {Atom_brdr3       = at.Atom_brdr3;    is_Atom_brdr3  = 1;}


   if(at.is_Atom_such_n)      {Atom_such_n      = at.Atom_such_n;   is_Atom_such_n = 1;}
   if(at.is_Atom_such_m)      {Atom_such_m      = at.Atom_such_m;   is_Atom_such_m = 1;}
   if(at.is_Atom_such_a)      {Atom_such_a      = at.Atom_such_a;   is_Atom_such_a = 1;}
   if(at.is_Atom_such_D)      {Atom_such_D      = at.Atom_such_D;   is_Atom_such_D = 1;}
   if(at.is_Atom_such_c)      {Atom_such_c      = at.Atom_such_c;   is_Atom_such_c = 1;}

   return *this;

}

void Atom_Record::merge(const Atom_Record& at){

   if(!is_Atom_ff_type && at.is_Atom_ff_type)         {Atom_ff_type     = at.Atom_ff_type;  is_Atom_ff_type = 1;}
   if(!is_Atom_ff_type_H && at.is_Atom_ff_type_H)     {Atom_ff_type_H   = at.Atom_ff_type_H;is_Atom_ff_type_H = 1;}
   if(!is_Atom_ff_int_type && at.is_Atom_ff_int_type) {Atom_ff_int_type = at.Atom_ff_int_type; is_Atom_ff_int_type = 1;}
   if(!is_Atom_element && at.is_Atom_element)         {Atom_element     = at.Atom_element;  is_Atom_element = 1;}
   if(!is_Atom_atomic_number && at.is_Atom_atomic_number){Atom_atomic_number = at.Atom_atomic_number;  is_Atom_atomic_number = 1;}

   if(!is_Atom_electronegativity && at.is_Atom_electronegativity){Atom_electronegativity = at.Atom_electronegativity;  is_Atom_electronegativity = 1;}
   if(!is_Atom_partial_charge && at.is_Atom_partial_charge){Atom_partial_charge = at.Atom_partial_charge;  is_Atom_partial_charge = 1;}

   if(!is_Atom_ff_eq_int_type2 && at.is_Atom_ff_eq_int_type2) { Atom_ff_eq_int_type2 = at.Atom_ff_eq_int_type2; is_Atom_ff_eq_int_type2 = 1; }
   if(!is_Atom_ff_eq_int_type3 && at.is_Atom_ff_eq_int_type3) { Atom_ff_eq_int_type3 = at.Atom_ff_eq_int_type3; is_Atom_ff_eq_int_type3 = 1; }
   if(!is_Atom_ff_eq_int_type4 && at.is_Atom_ff_eq_int_type4) { Atom_ff_eq_int_type4 = at.Atom_ff_eq_int_type4; is_Atom_ff_eq_int_type4 = 1; }
   if(!is_Atom_ff_eq_int_type5 && at.is_Atom_ff_eq_int_type5) { Atom_ff_eq_int_type5 = at.Atom_ff_eq_int_type5; is_Atom_ff_eq_int_type5 = 1; }


   if(!is_Atom_radius && at.is_Atom_radius)      {Atom_radius      = at.Atom_radius;   is_Atom_radius = 1;}
   if(!is_Atom_Z_star && at.is_Atom_Z_star)      {Atom_Z_star      = at.Atom_Z_star;   is_Atom_Z_star = 1;}
   if(!is_Atom_theta && at.is_Atom_theta)        {Atom_theta       = at.Atom_theta;    is_Atom_theta  = 1;}
   if(!is_Atom_sigma && at.is_Atom_sigma)        {Atom_sigma       = at.Atom_sigma;    is_Atom_sigma  = 1;}
   if(!is_Atom_epsilon && at.is_Atom_epsilon)    {Atom_epsilon     = at.Atom_epsilon;  is_Atom_epsilon = 1;}
   if(!is_Atom_GMP && at.is_Atom_GMP)            {Atom_GMP         = at.Atom_GMP;      is_Atom_GMP    = 1;}

   if(!is_Atom_crd && at.is_Atom_crd)            {Atom_crd         = at.Atom_crd;      is_Atom_crd    = 1;}
   if(!is_Atom_val && at.is_Atom_val)            {Atom_val         = at.Atom_val;      is_Atom_val    = 1;}
   if(!is_Atom_pilp && at.is_Atom_pilp)          {Atom_pilp        = at.Atom_pilp;     is_Atom_pilp   = 1;}
   if(!is_Atom_mltb && at.is_Atom_mltb)          {Atom_mltb        = at.Atom_mltb;     is_Atom_mltb   = 1;}
   if(!is_Atom_arom && at.is_Atom_arom)          {Atom_arom        = at.Atom_arom;     is_Atom_arom   = 1;}
   if(!is_Atom_lin && at.is_Atom_lin)            {Atom_lin         = at.Atom_lin;      is_Atom_lin    = 1;}
   if(!is_Atom_sbmb && at.is_Atom_sbmb)          {Atom_sbmb        = at.Atom_sbmb;     is_Atom_sbmb   = 1;}
   if(!is_Atom_alpha && at.is_Atom_alpha)        {Atom_alpha       = at.Atom_alpha;    is_Atom_alpha  = 1;}
   if(!is_Atom_N_eff && at.is_Atom_N_eff)        {Atom_N_eff       = at.Atom_N_eff;    is_Atom_N_eff  = 1;}
   if(!is_Atom_A_scale && at.is_Atom_A_scale)    {Atom_A_scale     = at.Atom_A_scale;  is_Atom_A_scale= 1;}
   if(!is_Atom_G_scale && at.is_Atom_G_scale)    {Atom_G_scale     = at.Atom_G_scale;  is_Atom_G_scale= 1;}
   if(!is_Atom_DAN && at.is_Atom_DAN)            {Atom_DAN         = at.Atom_DAN;      is_Atom_DAN= 1;}
   if(!is_Atom_pbci && at.is_Atom_pbci)          {Atom_pbci        = at.Atom_pbci;     is_Atom_pbci= 1;}
   if(!is_Atom_fcadj && at.is_Atom_fcadj)        {Atom_fcadj       = at.Atom_fcadj;    is_Atom_fcadj= 1;}

   if(!is_Atom_dative && at.is_Atom_dative)      {Atom_dative      = at.Atom_dative;   is_Atom_dative = 1;}
   if(!is_Atom_brdr1 && at.is_Atom_brdr1)        {Atom_brdr1       = at.Atom_brdr1;    is_Atom_brdr1  = 1;}
   if(!is_Atom_brdr2 && at.is_Atom_brdr2)        {Atom_brdr2       = at.Atom_brdr2;    is_Atom_brdr2  = 1;}
   if(!is_Atom_brdr3 && at.is_Atom_brdr3)        {Atom_brdr3       = at.Atom_brdr3;    is_Atom_brdr3  = 1;}


   if(!is_Atom_such_n && at.is_Atom_such_n)      {Atom_such_n      = at.Atom_such_n;   is_Atom_such_n = 1;}
   if(!is_Atom_such_m && at.is_Atom_such_m)      {Atom_such_m      = at.Atom_such_m;   is_Atom_such_m = 1;}
   if(!is_Atom_such_a && at.is_Atom_such_a)      {Atom_such_a      = at.Atom_such_a;   is_Atom_such_a = 1;}
   if(!is_Atom_such_D && at.is_Atom_such_D)      {Atom_such_D      = at.Atom_such_D;   is_Atom_such_D = 1;}
   if(!is_Atom_such_c && at.is_Atom_such_c)      {Atom_such_c      = at.Atom_such_c;   is_Atom_such_c = 1;}

}

void Atom_Record::save(boost::property_tree::ptree& pt,std::string path){

  if(is_Atom_ff_type){  ::save(pt,path+".Atom_ff_type",Atom_ff_type);    }
  if(is_Atom_ff_type_H){  ::save(pt,path+".Atom_ff_type_H",Atom_ff_type_H);    }
  if(is_Atom_ff_int_type){  ::save(pt,path+".Atom_ff_int_type",Atom_ff_int_type);    }
  if(is_Atom_element){  ::save(pt,path+".Atom_element",Atom_element);    }
  if(is_Atom_atomic_number){  ::save(pt,path+".Atom_atomic_number",Atom_atomic_number);    }
  if(is_Atom_electronegativity){  ::save(pt,path+".Atom_electronegativity",Atom_electronegativity);    }
  if(is_Atom_partial_charge){  ::save(pt,path+".Atom_partial_charge",Atom_partial_charge);    }
  if(is_Atom_ff_eq_int_type2){  ::save(pt,path+".Atom_ff_eq_int_type2",Atom_ff_eq_int_type2);    }
  if(is_Atom_ff_eq_int_type3){  ::save(pt,path+".Atom_ff_eq_int_type3",Atom_ff_eq_int_type3);    }
  if(is_Atom_ff_eq_int_type4){  ::save(pt,path+".Atom_ff_eq_int_type4",Atom_ff_eq_int_type4);    }
  if(is_Atom_ff_eq_int_type5){  ::save(pt,path+".Atom_ff_eq_int_type5",Atom_ff_eq_int_type5);    }
  if(is_Atom_radius){  ::save(pt,path+".Atom_radius",Atom_radius);    }
  if(is_Atom_Z_star){  ::save(pt,path+".Atom_Z_star",Atom_Z_star);    }
  if(is_Atom_theta){  ::save(pt,path+".Atom_theta",Atom_theta);    }
  if(is_Atom_sigma){  ::save(pt,path+".Atom_sigma",Atom_sigma);    }
  if(is_Atom_epsilon){  ::save(pt,path+".Atom_epsilon",Atom_epsilon);    }
  if(is_Atom_GMP){  ::save(pt,path+".Atom_GMP",Atom_GMP);    }
  if(is_Atom_crd){  ::save(pt,path+".Atom_crd",Atom_crd);    }
  if(is_Atom_val){  ::save(pt,path+".Atom_val",Atom_val);    }
  if(is_Atom_pilp){  ::save(pt,path+".Atom_pilp",Atom_pilp);    }
  if(is_Atom_mltb){  ::save(pt,path+".Atom_mltb",Atom_mltb);    }
  if(is_Atom_arom){  ::save(pt,path+".Atom_arom",Atom_arom);    }
  if(is_Atom_lin){  ::save(pt,path+".Atom_lin",Atom_lin);    }
  if(is_Atom_sbmb){  ::save(pt,path+".Atom_sbmb",Atom_sbmb);    }
  if(is_Atom_alpha){  ::save(pt,path+".Atom_alpha",Atom_alpha);    }
  if(is_Atom_N_eff){  ::save(pt,path+".Atom_N_eff",Atom_N_eff);    }
  if(is_Atom_A_scale){  ::save(pt,path+".Atom_A_scale",Atom_A_scale);    }
  if(is_Atom_G_scale){  ::save(pt,path+".Atom_G_scale",Atom_G_scale);    }
  if(is_Atom_DAN){  ::save(pt,path+".Atom_DAN",Atom_DAN);    }
  if(is_Atom_pbci){  ::save(pt,path+".Atom_pbci",Atom_pbci);    }
  if(is_Atom_fcadj){  ::save(pt,path+".Atom_fcadj",Atom_fcadj);    }
  if(is_Atom_dative){  ::save(pt,path+".Atom_dative",Atom_dative);    }
  if(is_Atom_brdr1){  ::save(pt,path+".Atom_brdr1",Atom_brdr1);    }
  if(is_Atom_brdr2){  ::save(pt,path+".Atom_brdr2",Atom_brdr2);    }
  if(is_Atom_brdr3){  ::save(pt,path+".Atom_brdr3",Atom_brdr3);    }
  if(is_Atom_such_n){  ::save(pt,path+".Atom_such_n",Atom_such_n);    }
  if(is_Atom_such_m){  ::save(pt,path+".Atom_such_m",Atom_such_m);    }
  if(is_Atom_such_a){  ::save(pt,path+".Atom_such_a",Atom_such_a);    }
  if(is_Atom_such_D){  ::save(pt,path+".Atom_such_D",Atom_such_D);    }
  if(is_Atom_such_c){  ::save(pt,path+".Atom_such_c",Atom_such_c);    }


}

void save(boost::property_tree::ptree& pt,std::string path,vector<Atom_Record>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".Atom_Record"+rt);
  }
}


void Atom_Record::load(boost::property_tree::ptree& pt,std::string path,int& status){
  int st;
  status = 0;

  ::load(pt,path+".Atom_ff_type",Atom_ff_type,is_Atom_ff_type); if(is_Atom_ff_type==1) { status=1;}
  ::load(pt,path+".Atom_ff_type_H",Atom_ff_type_H,is_Atom_ff_type_H); if(is_Atom_ff_type_H==1) { status=1;}
  ::load(pt,path+".Atom_ff_int_type",Atom_ff_int_type,is_Atom_ff_int_type); if(is_Atom_ff_int_type==1) { status=1;}
  ::load(pt,path+".Atom_element",Atom_element,is_Atom_element); if(is_Atom_element==1) { status=1;}
  ::load(pt,path+".Atom_atomic_number",Atom_atomic_number,is_Atom_atomic_number); if(is_Atom_atomic_number==1) { status=1;}

  ::load(pt,path+".Atom_electronegativity",Atom_electronegativity,is_Atom_electronegativity); if(is_Atom_electronegativity==1) { status=1;}
  ::load(pt,path+".Atom_partial_charge",Atom_partial_charge,is_Atom_partial_charge); if(is_Atom_partial_charge==1) { status=1;}
  ::load(pt,path+".Atom_ff_eq_int_type2",Atom_ff_eq_int_type2,is_Atom_ff_eq_int_type2); if(is_Atom_ff_eq_int_type2==1) { status=1;}
  ::load(pt,path+".Atom_ff_eq_int_type3",Atom_ff_eq_int_type3,is_Atom_ff_eq_int_type3); if(is_Atom_ff_eq_int_type3==1) { status=1;}
  ::load(pt,path+".Atom_ff_eq_int_type4",Atom_ff_eq_int_type4,is_Atom_ff_eq_int_type4); if(is_Atom_ff_eq_int_type4==1) { status=1;}
  ::load(pt,path+".Atom_ff_eq_int_type5",Atom_ff_eq_int_type5,is_Atom_ff_eq_int_type5); if(is_Atom_ff_eq_int_type5==1) { status=1;}

  ::load(pt,path+".Atom_radius",Atom_radius,is_Atom_radius); if(is_Atom_radius==1) { status=1;}
  ::load(pt,path+".Atom_Z_star",Atom_Z_star,is_Atom_Z_star); if(is_Atom_Z_star==1) { status=1;}
  ::load(pt,path+".Atom_theta",Atom_theta,is_Atom_theta); if(is_Atom_theta==1) { status=1;}
  ::load(pt,path+".Atom_sigma",Atom_sigma,is_Atom_sigma); if(is_Atom_sigma==1) { status=1;}
  ::load(pt,path+".Atom_epsilon",Atom_epsilon,is_Atom_epsilon); if(is_Atom_epsilon==1) { status=1;}
  ::load(pt,path+".Atom_GMP",Atom_GMP,is_Atom_GMP); if(is_Atom_GMP==1) { status=1;}

  ::load(pt,path+".Atom_crd",Atom_crd,is_Atom_crd); if(is_Atom_crd==1) { status=1;}
  ::load(pt,path+".Atom_val",Atom_val,is_Atom_val); if(is_Atom_val==1) { status=1;}
  ::load(pt,path+".Atom_pilp",Atom_pilp,is_Atom_pilp); if(is_Atom_pilp==1) { status=1;}
  ::load(pt,path+".Atom_mltb",Atom_mltb,is_Atom_mltb); if(is_Atom_mltb==1) { status=1;}
  ::load(pt,path+".Atom_arom",Atom_arom,is_Atom_arom); if(is_Atom_arom==1) { status=1;}
  ::load(pt,path+".Atom_lin",Atom_lin,is_Atom_lin); if(is_Atom_lin==1) { status=1;}
  ::load(pt,path+".Atom_sbmb",Atom_sbmb,is_Atom_sbmb); if(is_Atom_sbmb==1) { status=1;}
  ::load(pt,path+".Atom_alpha",Atom_alpha,is_Atom_alpha); if(is_Atom_alpha==1) { status=1;}
  ::load(pt,path+".Atom_N_eff",Atom_N_eff,is_Atom_N_eff); if(is_Atom_N_eff==1) { status=1;}
  ::load(pt,path+".Atom_A_scale",Atom_A_scale,is_Atom_A_scale); if(is_Atom_A_scale==1) { status=1;}
  ::load(pt,path+".Atom_G_scale",Atom_G_scale,is_Atom_G_scale); if(is_Atom_G_scale==1) { status=1;}
  ::load(pt,path+".Atom_DAN",Atom_DAN,is_Atom_DAN); if(is_Atom_DAN==1) { status=1;}
  ::load(pt,path+".Atom_pbci",Atom_pbci,is_Atom_pbci); if(is_Atom_pbci==1) { status=1;}
  ::load(pt,path+".Atom_fcadj",Atom_fcadj,is_Atom_fcadj); if(is_Atom_fcadj==1) { status=1;}

  ::load(pt,path+".Atom_dative",Atom_dative,is_Atom_dative); if(is_Atom_dative==1) { status=1;}
  ::load(pt,path+".Atom_brdr1",Atom_brdr1,is_Atom_brdr1); if(is_Atom_brdr1==1) { status=1;}
  ::load(pt,path+".Atom_brdr2",Atom_brdr2,is_Atom_brdr2); if(is_Atom_brdr2==1) { status=1;}
  ::load(pt,path+".Atom_brdr3",Atom_brdr3,is_Atom_brdr3); if(is_Atom_brdr3==1) { status=1;}

  ::load(pt,path+".Atom_such_n",Atom_such_n,is_Atom_such_n); if(is_Atom_such_n==1) { status=1;}
  ::load(pt,path+".Atom_such_m",Atom_such_m,is_Atom_such_m); if(is_Atom_such_m==1) { status=1;}
  ::load(pt,path+".Atom_such_a",Atom_such_a,is_Atom_such_a); if(is_Atom_such_a==1) { status=1;}
  ::load(pt,path+".Atom_such_D",Atom_such_D,is_Atom_such_D); if(is_Atom_such_D==1) { status=1;}
  ::load(pt,path+".Atom_such_c",Atom_such_c,is_Atom_such_c); if(is_Atom_such_c==1) { status=1;}


}

void load(boost::property_tree::ptree& pt,std::string path,vector<Atom_Record>& vt,int& status){
  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      Atom_Record x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}



//-------------------------- Bond Record Class members -------------------------
int Bond_Record::set(object at){

   set_value(is_Atom1_ff_type,    Atom1_ff_type,    at,"Atom1_ff_type");
   set_value(is_Atom1_ff_int_type,Atom1_ff_int_type,at,"Atom1_ff_int_type");
   set_value(is_Atom2_ff_type,    Atom2_ff_type,    at,"Atom2_ff_type");
   set_value(is_Atom2_ff_int_type,Atom2_ff_int_type,at,"Atom2_ff_int_type");
   set_value(is_Atom1_atomic_number, Atom1_atomic_number, at,"Atom1_atomic_number");
   set_value(is_Atom2_atomic_number, Atom2_atomic_number, at,"Atom2_atomic_number");
   set_value(is_Atom1_element, Atom1_element, at,"Atom1_element");
   set_value(is_Atom2_element, Atom2_element, at,"Atom2_element");
   set_value(is_Bond_type_index,  Bond_type_index,  at,"Bond_type_index");

   set_value(is_Bond_r_eq,   Bond_r_eq,   at, "Bond_r_eq");
   set_value(is_Bond_k_bond, Bond_k_bond, at, "Bond_k_bond"); 
   set_value(is_Bond_D_bond, Bond_D_bond, at, "Bond_D_bond");
   set_value(is_Bond_alpha,  Bond_alpha,  at, "Bond_alpha");  
   set_value(is_Bond_r_eq_ref,  Bond_r_eq_ref,  at, "Bond_r_eq_ref");
   set_value(is_Bond_k_bond_ref,Bond_k_bond_ref,at, "Bond_k_bond_ref");

   set_value(is_Bond_shift_elec,Bond_shift_elec,at, "Bond_shift_elec");

   set_value(is_Bond_bci,    Bond_bci,  at, "Bond_bci");

   set_value(is_Bond_wij,    Bond_wij,  at, "Bond_wij");
   set_value(is_Bond_wij_1,  Bond_wij_1,at, "Bond_wij_1");
   set_value(is_Bond_wij_2,  Bond_wij_2,at, "Bond_wij_2");
   set_value(is_Bond_alpij,  Bond_alpij,at, "Bond_alpij");
   set_value(is_Bond_alpij_1,Bond_alpij_1,at,"Bond_alpij_1");
   set_value(is_Bond_alpij_2,Bond_alpij_2,at,"Bond_alpij_2");
   set_value(is_Bond_rij_1,  Bond_rij_1,at, "Bond_rij_1");
   set_value(is_Bond_rij_2,  Bond_rij_2,at, "Bond_rij_2");
 
  
   return 1;

}
Bond_Record& Bond_Record::operator=(const Bond_Record& at){

   is_Atom1_ff_type = 0;
   is_Atom2_ff_type = 0;
   is_Atom1_ff_int_type = 0;
   is_Atom2_ff_int_type = 0;
   is_Atom1_atomic_number = 0;
   is_Atom2_atomic_number = 0;
   is_Atom1_element = 0;
   is_Atom2_element = 0;

   is_Bond_type_index = 0;

   is_Bond_r_eq = 0;
   is_Bond_k_bond = 0;
   is_Bond_D_bond = 0;
   is_Bond_alpha = 0;
   is_Bond_r_eq_ref = 0;
   is_Bond_k_bond_ref = 0;
   is_Bond_shift_elec = 0;
   is_Bond_bci = 0;
   is_Bond_wij = 0;
   is_Bond_wij_1 = 0;
   is_Bond_wij_2 = 0;
   is_Bond_alpij = 0;
   is_Bond_alpij_1 = 0;
   is_Bond_alpij_2 = 0;
   is_Bond_rij_1 = 0;
   is_Bond_rij_2 = 0;


   if(at.is_Atom1_ff_type)     {Atom1_ff_type     = at.Atom1_ff_type;     is_Atom1_ff_type = 1;}
   if(at.is_Atom1_ff_int_type) {Atom1_ff_int_type = at.Atom1_ff_int_type; is_Atom1_ff_int_type = 1;}
   if(at.is_Atom2_ff_type)     {Atom2_ff_type     = at.Atom2_ff_type;     is_Atom2_ff_type = 1;}
   if(at.is_Atom2_ff_int_type) {Atom2_ff_int_type = at.Atom2_ff_int_type; is_Atom2_ff_int_type = 1;}
   if(at.is_Atom1_atomic_number) {Atom1_atomic_number = at.Atom1_atomic_number; is_Atom1_atomic_number = 1;}
   if(at.is_Atom2_atomic_number) {Atom2_atomic_number = at.Atom2_atomic_number; is_Atom2_atomic_number = 1;}
   if(at.is_Atom1_element) {Atom1_element = at.Atom1_element; is_Atom1_element = 1;}
   if(at.is_Atom2_element) {Atom2_element = at.Atom2_element; is_Atom2_element = 1;}


   if(at.is_Bond_type_index)   {Bond_type_index   = at.Bond_type_index;   is_Bond_type_index = 1; }

   if(at.is_Bond_r_eq)      { Bond_r_eq     = at.Bond_r_eq;     is_Bond_r_eq = 1; }
   if(at.is_Bond_k_bond)    { Bond_k_bond   = at.Bond_k_bond;   is_Bond_k_bond = 1; }
   if(at.is_Bond_D_bond)    { Bond_D_bond   = at.Bond_D_bond;   is_Bond_D_bond = 1; }
   if(at.is_Bond_alpha)     { Bond_alpha    = at.Bond_alpha;    is_Bond_alpha = 1; }
   if(at.is_Bond_r_eq_ref)  { Bond_r_eq_ref   = at.Bond_r_eq_ref;   is_Bond_r_eq_ref   = 1; }
   if(at.is_Bond_k_bond_ref){ Bond_k_bond_ref = at.Bond_k_bond_ref; is_Bond_k_bond_ref = 1; }
   if(at.is_Bond_shift_elec){ Bond_shift_elec = at.Bond_shift_elec; is_Bond_shift_elec = 1; }
   if(at.is_Bond_bci)       { Bond_bci      = at.Bond_bci;      is_Bond_bci = 1; }

   if(at.is_Bond_wij)       { Bond_wij      = at.Bond_wij;      is_Bond_wij = 1; }
   if(at.is_Bond_wij_1)     { Bond_wij_1    = at.Bond_wij_1;    is_Bond_wij_1 = 1; }
   if(at.is_Bond_wij_2)     { Bond_wij_2    = at.Bond_wij_2;    is_Bond_wij_2 = 1; }
   if(at.is_Bond_alpij)     { Bond_alpij    = at.Bond_alpij;    is_Bond_alpij = 1; }
   if(at.is_Bond_alpij_1)   { Bond_alpij_1  = at.Bond_alpij_1;  is_Bond_alpij_1 = 1; }
   if(at.is_Bond_alpij_2)   { Bond_alpij_2  = at.Bond_alpij_2;  is_Bond_alpij_2 = 1; }
   if(at.is_Bond_rij_1)     { Bond_rij_1    = at.Bond_rij_1;    is_Bond_rij_1 = 1; }
   if(at.is_Bond_rij_2)     { Bond_rij_2    = at.Bond_rij_2;    is_Bond_rij_2 = 1; }


   return *this;
}

void Bond_Record::merge(const Bond_Record& at){

   if(!is_Atom1_ff_type && at.is_Atom1_ff_type)         {Atom1_ff_type     = at.Atom1_ff_type;     is_Atom1_ff_type = 1;}
   if(!is_Atom1_ff_int_type && at.is_Atom1_ff_int_type) {Atom1_ff_int_type = at.Atom1_ff_int_type; is_Atom1_ff_int_type = 1;}
   if(!is_Atom2_ff_type && at.is_Atom2_ff_type)         {Atom2_ff_type     = at.Atom2_ff_type;     is_Atom2_ff_type = 1;}
   if(!is_Atom2_ff_int_type && at.is_Atom2_ff_int_type) {Atom2_ff_int_type = at.Atom2_ff_int_type; is_Atom2_ff_int_type = 1;}
   if(!is_Atom1_atomic_number && at.is_Atom1_atomic_number) {Atom1_atomic_number = at.Atom1_atomic_number; is_Atom1_atomic_number = 1;}
   if(!is_Atom2_atomic_number && at.is_Atom2_atomic_number) {Atom2_atomic_number = at.Atom2_atomic_number; is_Atom2_atomic_number = 1;}
   if(!is_Atom1_element && at.is_Atom1_element) {Atom1_element = at.Atom1_element; is_Atom1_element = 1;}
   if(!is_Atom2_element && at.is_Atom2_element) {Atom2_element = at.Atom2_element; is_Atom2_element = 1;}


   if(!is_Bond_type_index && at.is_Bond_type_index)     {Bond_type_index   = at.Bond_type_index;   is_Bond_type_index = 1; }

   if(!is_Bond_r_eq && at.is_Bond_r_eq)        { Bond_r_eq     = at.Bond_r_eq;     is_Bond_r_eq = 1; }
   if(!is_Bond_k_bond && at.is_Bond_k_bond)    { Bond_k_bond   = at.Bond_k_bond;   is_Bond_k_bond = 1; }
   if(!is_Bond_D_bond && at.is_Bond_D_bond)    { Bond_D_bond   = at.Bond_D_bond;   is_Bond_D_bond = 1; }
   if(!is_Bond_alpha && at.is_Bond_alpha)      { Bond_alpha    = at.Bond_alpha;    is_Bond_alpha = 1; }

   if(!is_Bond_r_eq_ref && at.is_Bond_r_eq_ref)  { Bond_r_eq_ref   = at.Bond_r_eq_ref;   is_Bond_r_eq_ref   = 1; }
   if(!is_Bond_k_bond_ref && at.is_Bond_k_bond_ref){ Bond_k_bond_ref = at.Bond_k_bond_ref; is_Bond_k_bond_ref = 1; }

   if(!is_Bond_shift_elec && at.is_Bond_shift_elec){ Bond_shift_elec = at.Bond_shift_elec; is_Bond_shift_elec = 1; }

   if(!is_Bond_bci && at.is_Bond_bci)          { Bond_bci      = at.Bond_bci;      is_Bond_bci = 1; }

   if(!is_Bond_wij && at.is_Bond_wij)         { Bond_wij      = at.Bond_wij;      is_Bond_wij = 1; }
   if(!is_Bond_wij_1 && at.is_Bond_wij_1)     { Bond_wij_1    = at.Bond_wij_1;    is_Bond_wij_1 = 1; }
   if(!is_Bond_wij_2 && at.is_Bond_wij_2)     { Bond_wij_2    = at.Bond_wij_2;    is_Bond_wij_2 = 1; }
   if(!is_Bond_alpij && at.is_Bond_alpij)     { Bond_alpij    = at.Bond_alpij;    is_Bond_alpij = 1; }
   if(!is_Bond_alpij_1 && at.is_Bond_alpij_1) { Bond_alpij_1  = at.Bond_alpij_1;  is_Bond_alpij_1 = 1; }
   if(!is_Bond_alpij_2 && at.is_Bond_alpij_2) { Bond_alpij_2  = at.Bond_alpij_2;  is_Bond_alpij_2 = 1; }
   if(!is_Bond_rij_1 && at.is_Bond_rij_1)     { Bond_rij_1    = at.Bond_rij_1;    is_Bond_rij_1 = 1; }
   if(!is_Bond_rij_2 && at.is_Bond_rij_2)     { Bond_rij_2    = at.Bond_rij_2;    is_Bond_rij_2 = 1; }



}


int Bond_Record::show_info(){

   std::cout<<"Bond_Record properties:"<<std::endl;
   if(is_Atom1_ff_type)        {std::cout<<"Atom1_ff_type = "<<Atom1_ff_type<<std::endl;   }
   if(is_Atom1_ff_int_type)    {std::cout<<"Atom1_ff_int_type = "<<Atom1_ff_int_type<<std::endl;     }
   if(is_Atom2_ff_type)        {std::cout<<"Atom2_ff_type = "<<Atom2_ff_type<<std::endl;   }
   if(is_Atom2_ff_int_type)    {std::cout<<"Atom2_ff_int_type = "<<Atom2_ff_int_type<<std::endl;     }
   if(is_Atom1_atomic_number)  {std::cout<<"Atom1_atomic_number = "<<Atom1_atomic_number<<std::endl;   }
   if(is_Atom2_atomic_number)  {std::cout<<"Atom2_atomic_number = "<<Atom2_atomic_number<<std::endl;   }
   if(is_Atom1_element)        {std::cout<<"Atom1_element = "<<Atom1_element<<std::endl;   }
   if(is_Atom2_element)        {std::cout<<"Atom2_element = "<<Atom2_element<<std::endl;   }


   if(is_Bond_type_index)      {std::cout<<"Bond_type_index = "<<Bond_type_index<<std::endl; }

   if(is_Bond_r_eq)            {std::cout<<"Bond_r_eq = "<<Bond_r_eq<<" Angstrom"<<std::endl; }
   if(is_Bond_k_bond)          {std::cout<<"Bond_k_bond = "<<Bond_k_bond<<" kcal/(mol*Angst^2)"<<std::endl; }
   if(is_Bond_D_bond)          {std::cout<<"Bond_D_bond = "<<Bond_D_bond<<" kcal/mol"<<std::endl; }
   if(is_Bond_alpha)           {std::cout<<"Bond_alpha = "<<Bond_alpha<<" Angst^-1"<<std::endl; }
   if(is_Bond_r_eq_ref)        {std::cout<<"Bond_r_eq_ref = "<<Bond_r_eq_ref<<std::endl; }
   if(is_Bond_k_bond_ref)      {std::cout<<"Bond_k_bond_ref = "<<Bond_k_bond_ref<<std::endl; }

   if(is_Bond_shift_elec)      {std::cout<<"Bond_shift_elec = "<<Bond_shift_elec<<std::endl; }

   if(is_Bond_bci)             {std::cout<<"Bond_bci = "<<Bond_bci<<std::endl; }

   if(is_Bond_wij)             {std::cout<<"Bond_wij = "<<Bond_wij<<std::endl; }
   if(is_Bond_wij_1)           {std::cout<<"Bond_wij_1 = "<<Bond_wij_1<<std::endl; }
   if(is_Bond_wij_2)           {std::cout<<"Bond_wij_2 = "<<Bond_wij_2<<std::endl; }
   if(is_Bond_alpij)           {std::cout<<"Bond_alpij = "<<Bond_alpij<<std::endl; }
   if(is_Bond_alpij_1)         {std::cout<<"Bond_alpij_1 = "<<Bond_alpij_1<<std::endl; }
   if(is_Bond_alpij_2)         {std::cout<<"Bond_alpij_2 = "<<Bond_alpij_2<<std::endl; }
   if(is_Bond_rij_1)           {std::cout<<"Bond_rij_1 = "<<Bond_rij_1<<std::endl; }
   if(is_Bond_rij_2)           {std::cout<<"Bond_rij_2 = "<<Bond_rij_2<<std::endl; }


   std::cout<<std::endl;

   return 1;

}

void Bond_Record::save(boost::property_tree::ptree& pt,std::string path){

  if(is_Atom1_ff_type){  ::save(pt,path+".Atom1_ff_type",Atom1_ff_type);    }
  if(is_Atom2_ff_type){  ::save(pt,path+".Atom2_ff_type",Atom2_ff_type);    }
  if(is_Atom1_ff_int_type){  ::save(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type);    }
  if(is_Atom2_ff_int_type){  ::save(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type);    }
  if(is_Atom1_atomic_number){  ::save(pt,path+".Atom1_atomic_number",Atom1_atomic_number);    }
  if(is_Atom2_atomic_number){  ::save(pt,path+".Atom2_atomic_number",Atom2_atomic_number);    }
  if(is_Atom1_element){  ::save(pt,path+".Atom1_element",Atom1_element);    }
  if(is_Atom2_element){  ::save(pt,path+".Atom2_element",Atom2_element);    }

  if(is_Bond_type_index){  ::save(pt,path+".Bond_type_index",Bond_type_index);    }
  if(is_Bond_r_eq){  ::save(pt,path+".Bond_r_eq",Bond_r_eq);    }
  if(is_Bond_k_bond){  ::save(pt,path+".Bond_k_bond",Bond_k_bond);    }
  if(is_Bond_D_bond){  ::save(pt,path+".Bond_D_bond",Bond_D_bond);    }
  if(is_Bond_alpha){  ::save(pt,path+".Bond_alpha",Bond_alpha);    }
  if(is_Bond_r_eq_ref){  ::save(pt,path+".Bond_r_eq_ref",Bond_r_eq_ref);    }
  if(is_Bond_k_bond_ref){  ::save(pt,path+".Bond_k_bond_ref",Bond_k_bond_ref);    }
  if(is_Bond_shift_elec){  ::save(pt,path+".Bond_shift_elec",Bond_shift_elec);    }
  if(is_Bond_bci){  ::save(pt,path+".Bond_bci",Bond_bci);    }

  if(is_Bond_wij){  ::save(pt,path+".Bond_wij",Bond_wij);    }
  if(is_Bond_wij_1){  ::save(pt,path+".Bond_wij_1",Bond_wij_1);    }
  if(is_Bond_wij_2){  ::save(pt,path+".Bond_wij_2",Bond_wij_2);    }
  if(is_Bond_alpij){  ::save(pt,path+".Bond_alpij",Bond_alpij);    }
  if(is_Bond_alpij_1){  ::save(pt,path+".Bond_alpij_1",Bond_alpij_1);    }
  if(is_Bond_alpij_2){  ::save(pt,path+".Bond_alpij_2",Bond_alpij_2);    }
  if(is_Bond_rij_1){  ::save(pt,path+".Bond_rij_1",Bond_rij_1);    }
  if(is_Bond_rij_2){  ::save(pt,path+".Bond_rij_2",Bond_rij_2);    }

}

void save(boost::property_tree::ptree& pt,std::string path,vector<Bond_Record>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".Bond_Record"+rt);
  }
}

void Bond_Record::load(boost::property_tree::ptree& pt,std::string path,int& status){
  int st;
  status = 0;

  ::load(pt,path+".Atom1_ff_type",Atom1_ff_type,is_Atom1_ff_type); if(is_Atom1_ff_type==1) { status=1;}
  ::load(pt,path+".Atom2_ff_type",Atom2_ff_type,is_Atom2_ff_type); if(is_Atom2_ff_type==1) { status=1;}
  ::load(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type,is_Atom1_ff_int_type); if(is_Atom1_ff_int_type==1) { status=1;}
  ::load(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type,is_Atom2_ff_int_type); if(is_Atom2_ff_int_type==1) { status=1;}
  ::load(pt,path+".Atom1_atomic_number",Atom1_atomic_number,is_Atom1_atomic_number); if(is_Atom1_atomic_number==1) { status=1;}
  ::load(pt,path+".Atom2_atomic_number",Atom2_atomic_number,is_Atom2_atomic_number); if(is_Atom2_atomic_number==1) { status=1;}
  ::load(pt,path+".Atom1_element",Atom1_element,is_Atom1_element); if(is_Atom1_element==1) { status=1;}
  ::load(pt,path+".Atom2_element",Atom2_element,is_Atom2_element); if(is_Atom2_element==1) { status=1;}

  ::load(pt,path+".Bond_type_index",Bond_type_index,is_Bond_type_index); if(is_Bond_type_index==1) { status=1;}
  ::load(pt,path+".Bond_r_eq",Bond_r_eq,is_Bond_r_eq); if(is_Bond_r_eq==1) { status=1;}
  ::load(pt,path+".Bond_k_bond",Bond_k_bond,is_Bond_k_bond); if(is_Bond_k_bond==1) { status=1;}
  ::load(pt,path+".Bond_D_bond",Bond_D_bond,is_Bond_D_bond); if(is_Bond_D_bond==1) { status=1;}
  ::load(pt,path+".Bond_alpha",Bond_alpha,is_Bond_alpha); if(is_Bond_alpha==1) { status=1;}
  ::load(pt,path+".Bond_r_eq_ref",Bond_r_eq_ref,is_Bond_r_eq_ref); if(is_Bond_r_eq_ref==1) { status=1;}
  ::load(pt,path+".Bond_k_bond_ref",Bond_k_bond_ref,is_Bond_k_bond_ref); if(is_Bond_k_bond_ref==1) { status=1;}
  ::load(pt,path+".Bond_shift_elec",Bond_shift_elec,is_Bond_shift_elec); if(is_Bond_shift_elec==1) { status=1;}
  ::load(pt,path+".Bond_bci",Bond_bci,is_Bond_bci); if(is_Bond_bci==1) { status=1;}

  ::load(pt,path+".Bond_wij",Bond_wij,is_Bond_wij); if(is_Bond_wij==1) { status=1;}
  ::load(pt,path+".Bond_wij_1",Bond_wij_1,is_Bond_wij_1); if(is_Bond_wij_1==1) { status=1;}
  ::load(pt,path+".Bond_wij_2",Bond_wij_2,is_Bond_wij_2); if(is_Bond_wij_2==1) { status=1;}
  ::load(pt,path+".Bond_alpij",Bond_alpij,is_Bond_alpij); if(is_Bond_alpij==1) { status=1;}
  ::load(pt,path+".Bond_alpij_1",Bond_alpij_1,is_Bond_alpij_1); if(is_Bond_alpij_1==1) { status=1;}
  ::load(pt,path+".Bond_alpij_2",Bond_alpij_2,is_Bond_alpij_2); if(is_Bond_alpij_2==1) { status=1;}
  ::load(pt,path+".Bond_rij_1",Bond_rij_1,is_Bond_rij_1); if(is_Bond_rij_1==1) { status=1;}
  ::load(pt,path+".Bond_rij_2",Bond_rij_2,is_Bond_rij_2); if(is_Bond_rij_2==1) { status=1;}

}


void load(boost::property_tree::ptree& pt,std::string path,vector<Bond_Record>& vt,int& status){
  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      Bond_Record x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}


//----------------------- Angle Record Class members ------------------------------
int Angle_Record::set(object at){

   set_value(is_Atom1_ff_type,    Atom1_ff_type,    at,"Atom1_ff_type");
   set_value(is_Atom1_ff_int_type,Atom1_ff_int_type,at,"Atom1_ff_int_type");
   set_value(is_Atom2_ff_type,    Atom2_ff_type,    at,"Atom2_ff_type");
   set_value(is_Atom2_ff_int_type,Atom2_ff_int_type,at,"Atom2_ff_int_type");
   set_value(is_Atom3_ff_type,    Atom3_ff_type,    at,"Atom3_ff_type");
   set_value(is_Atom3_ff_int_type,Atom3_ff_int_type,at,"Atom3_ff_int_type");
   set_value(is_Angle_type_index, Angle_type_index, at,"Angle_type_index");


   set_value(is_Angle_theta_eq,   Angle_theta_eq,   at, "Angle_theta_eq");
   set_value(is_Angle_k_angle,    Angle_k_angle,    at, "Angle_k_angle");
   set_value(is_Angle_r_eq,       Angle_r_eq,       at, "Angle_r_eq");
   set_value(is_Angle_k_ub,       Angle_k_ub,       at, "Angle_k_ub");
   set_value(is_Angle_kijk_sb,    Angle_kijk_sb,    at, "Angle_kijk_sb");
   set_value(is_Angle_kkji_sb,    Angle_kkji_sb,    at, "Angle_kkji_sb");


   return 1;

}
Angle_Record& Angle_Record::operator=(const Angle_Record& at){

   is_Atom1_ff_type = 0;
   is_Atom2_ff_type = 0;
   is_Atom3_ff_type = 0;

   is_Atom1_ff_int_type = 0;
   is_Atom2_ff_int_type = 0;
   is_Atom3_ff_int_type = 0;

   is_Angle_type_index = 0;

   is_Angle_theta_eq = 0;
   is_Angle_k_angle  = 0;
   is_Angle_r_eq     = 0;
   is_Angle_k_ub     = 0;
   is_Angle_kijk_sb  = 0;
   is_Angle_kkji_sb  = 0;


   if(at.is_Atom1_ff_type)     {Atom1_ff_type     = at.Atom1_ff_type;     is_Atom1_ff_type = 1;}
   if(at.is_Atom1_ff_int_type) {Atom1_ff_int_type = at.Atom1_ff_int_type; is_Atom1_ff_int_type = 1;}
   if(at.is_Atom2_ff_type)     {Atom2_ff_type     = at.Atom2_ff_type;     is_Atom2_ff_type = 1;}
   if(at.is_Atom2_ff_int_type) {Atom2_ff_int_type = at.Atom2_ff_int_type; is_Atom2_ff_int_type = 1;}
   if(at.is_Atom3_ff_type)     {Atom3_ff_type     = at.Atom3_ff_type;     is_Atom3_ff_type = 1;}
   if(at.is_Atom3_ff_int_type) {Atom3_ff_int_type = at.Atom3_ff_int_type; is_Atom3_ff_int_type = 1;}
   if(at.is_Angle_type_index)  {Angle_type_index  = at.Angle_type_index;  is_Angle_type_index = 1;}


   if(at.is_Angle_theta_eq)    { Angle_theta_eq     = at.Angle_theta_eq;     is_Angle_theta_eq = 1; }
   if(at.is_Angle_k_angle)     { Angle_k_angle      = at.Angle_k_angle;      is_Angle_k_angle = 1; }
   if(at.is_Angle_r_eq)        { Angle_r_eq         = at.Angle_r_eq;         is_Angle_r_eq = 1; }
   if(at.is_Angle_k_ub)        { Angle_k_ub         = at.Angle_k_ub;         is_Angle_k_ub = 1;} 
   if(at.is_Angle_kijk_sb)     { Angle_kijk_sb      = at.Angle_kijk_sb;      is_Angle_kijk_sb = 1;}
   if(at.is_Angle_kkji_sb)     { Angle_kkji_sb      = at.Angle_kkji_sb;      is_Angle_kkji_sb = 1;}

   return *this;
}

void Angle_Record::merge(const Angle_Record& at){

   if(!is_Atom1_ff_type && at.is_Atom1_ff_type)         {Atom1_ff_type     = at.Atom1_ff_type;     is_Atom1_ff_type = 1;}
   if(!is_Atom1_ff_int_type && at.is_Atom1_ff_int_type) {Atom1_ff_int_type = at.Atom1_ff_int_type; is_Atom1_ff_int_type = 1;}
   if(!is_Atom2_ff_type && at.is_Atom2_ff_type)         {Atom2_ff_type     = at.Atom2_ff_type;     is_Atom2_ff_type = 1;}
   if(!is_Atom2_ff_int_type && at.is_Atom2_ff_int_type) {Atom2_ff_int_type = at.Atom2_ff_int_type; is_Atom2_ff_int_type = 1;}
   if(!is_Atom3_ff_type && at.is_Atom3_ff_type)         {Atom3_ff_type     = at.Atom3_ff_type;     is_Atom3_ff_type = 1;}
   if(!is_Atom3_ff_int_type && at.is_Atom3_ff_int_type) {Atom3_ff_int_type = at.Atom3_ff_int_type; is_Atom3_ff_int_type = 1;}
   if(!is_Angle_type_index && at.is_Angle_type_index)   {Angle_type_index  = at.Angle_type_index;  is_Angle_type_index = 1;}


   if(!is_Angle_theta_eq && at.is_Angle_theta_eq)   { Angle_theta_eq     = at.Angle_theta_eq;     is_Angle_theta_eq = 1; }
   if(!is_Angle_k_angle && at.is_Angle_k_angle)     { Angle_k_angle      = at.Angle_k_angle;      is_Angle_k_angle = 1; }
   if(!is_Angle_r_eq && at.is_Angle_r_eq)           { Angle_r_eq         = at.Angle_r_eq;         is_Angle_r_eq = 1; }
   if(!is_Angle_k_ub && at.is_Angle_k_ub)           { Angle_k_ub         = at.Angle_k_ub;         is_Angle_k_ub = 1;}
   if(!is_Angle_kijk_sb && at.is_Angle_kijk_sb)     { Angle_kijk_sb      = at.Angle_kijk_sb;      is_Angle_kijk_sb = 1;}
   if(!is_Angle_kkji_sb && at.is_Angle_kkji_sb)     { Angle_kkji_sb      = at.Angle_kkji_sb;      is_Angle_kkji_sb = 1;}

}


int Angle_Record::show_info(){

   std::cout<<"Angle_Record properties:"<<std::endl;
   if(is_Atom1_ff_type)        {std::cout<<"Atom1_ff_type = "<<Atom1_ff_type<<std::endl;   }
   if(is_Atom1_ff_int_type)    {std::cout<<"Atom1_ff_int_type = "<<Atom1_ff_int_type<<std::endl;     }
   if(is_Atom2_ff_type)        {std::cout<<"Atom2_ff_type = "<<Atom2_ff_type<<std::endl;   }
   if(is_Atom2_ff_int_type)    {std::cout<<"Atom2_ff_int_type = "<<Atom2_ff_int_type<<std::endl;     }
   if(is_Atom3_ff_type)        {std::cout<<"Atom3_ff_type = "<<Atom3_ff_type<<std::endl;   }
   if(is_Atom3_ff_int_type)    {std::cout<<"Atom3_ff_int_type = "<<Atom3_ff_int_type<<std::endl;     }
   if(is_Angle_type_index)     {std::cout<<"Angle_type_index = "<<Angle_type_index<<std::endl; }

   if(is_Angle_theta_eq)       {std::cout<<"Angle_theta_eq = "<<Angle_theta_eq<<" deg"<<std::endl; }
   if(is_Angle_k_angle)        {std::cout<<"Angle_k_angle = "<<Angle_k_angle<<" kcal/(mol*rad^2)"<<std::endl; }
   if(is_Angle_r_eq)           {std::cout<<"Angle_r_eq = "<<Angle_r_eq<<" Angstrom"<<std::endl; }
   if(is_Angle_k_ub)           {std::cout<<"Angle_k_ub = "<<Angle_k_ub<<" kcal/(mol*Angstrom^2)"<<std::endl; }
   if(is_Angle_kijk_sb)        {std::cout<<"Angle_kijk_sb = "<<Angle_kijk_sb<<" some units"<<std::endl; }
   if(is_Angle_kkji_sb)        {std::cout<<"Angle_kkji_sb = "<<Angle_kkji_sb<<" some units"<<std::endl; }


   std::cout<<std::endl;

   return 1;

}

void Angle_Record::save(boost::property_tree::ptree& pt,std::string path){

  if(is_Atom1_ff_type){  ::save(pt,path+".Atom1_ff_type",Atom1_ff_type);    }
  if(is_Atom2_ff_type){  ::save(pt,path+".Atom2_ff_type",Atom2_ff_type);    }
  if(is_Atom3_ff_type){  ::save(pt,path+".Atom3_ff_type",Atom3_ff_type);    }
  if(is_Atom1_ff_int_type){  ::save(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type);    }
  if(is_Atom2_ff_int_type){  ::save(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type);    }
  if(is_Atom3_ff_int_type){  ::save(pt,path+".Atom3_ff_int_type",Atom3_ff_int_type);    }

  if(is_Angle_type_index){  ::save(pt,path+".Angle_type_index",Angle_type_index);    }
  if(is_Angle_theta_eq){  ::save(pt,path+".Angle_theta_eq",Angle_theta_eq);    }
  if(is_Angle_k_angle){  ::save(pt,path+".Angle_k_angle",Angle_k_angle);    }
  if(is_Angle_r_eq){  ::save(pt,path+".Angle_r_eq",Angle_r_eq);    }
  if(is_Angle_k_ub){  ::save(pt,path+".Angle_k_ub",Angle_k_ub);    }
  if(is_Angle_kijk_sb){  ::save(pt,path+".Angle_kijk_sb",Angle_kijk_sb);    }
  if(is_Angle_kkji_sb){  ::save(pt,path+".Angle_kkji_sb",Angle_kkji_sb);    }

}

void save(boost::property_tree::ptree& pt,std::string path,vector<Angle_Record>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".Angle_Record"+rt);
  }
}

void Angle_Record::load(boost::property_tree::ptree& pt,std::string path,int& status){
  int st;
  status = 0;

  ::load(pt,path+".Atom1_ff_type",Atom1_ff_type,is_Atom1_ff_type); if(is_Atom1_ff_type==1) { status=1;}
  ::load(pt,path+".Atom2_ff_type",Atom2_ff_type,is_Atom2_ff_type); if(is_Atom2_ff_type==1) { status=1;}
  ::load(pt,path+".Atom3_ff_type",Atom3_ff_type,is_Atom3_ff_type); if(is_Atom3_ff_type==1) { status=1;}
  ::load(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type,is_Atom1_ff_int_type); if(is_Atom1_ff_int_type==1) { status=1;}
  ::load(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type,is_Atom2_ff_int_type); if(is_Atom2_ff_int_type==1) { status=1;}
  ::load(pt,path+".Atom3_ff_int_type",Atom3_ff_int_type,is_Atom3_ff_int_type); if(is_Atom3_ff_int_type==1) { status=1;}

  ::load(pt,path+".Angle_type_index",Angle_type_index,is_Angle_type_index); if(is_Angle_type_index==1) { status=1;}
  ::load(pt,path+".Angle_theta_eq",Angle_theta_eq,is_Angle_theta_eq); if(is_Angle_theta_eq==1) { status=1;}
  ::load(pt,path+".Angle_k_angle",Angle_k_angle,is_Angle_k_angle); if(is_Angle_k_angle==1) { status=1;}
  ::load(pt,path+".Angle_r_eq",Angle_r_eq,is_Angle_r_eq); if(is_Angle_r_eq==1) { status=1;}
  ::load(pt,path+".Angle_k_ub",Angle_k_ub,is_Angle_k_ub); if(is_Angle_k_ub==1) { status=1;}
  ::load(pt,path+".Angle_kijk_sb",Angle_kijk_sb,is_Angle_kijk_sb); if(is_Angle_kijk_sb==1) { status=1;}
  ::load(pt,path+".Angle_kkji_sb",Angle_kkji_sb,is_Angle_kkji_sb); if(is_Angle_kkji_sb==1) { status=1;}

}

void load(boost::property_tree::ptree& pt,std::string path,vector<Angle_Record>& vt,int& status){
  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      Angle_Record x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}


//----------------------- Dihedral Record Class members ------------------------------
int Dihedral_Record::set(object at){

   set_value(is_Atom1_ff_type,    Atom1_ff_type,    at,"Atom1_ff_type");
   set_value(is_Atom1_ff_int_type,Atom1_ff_int_type,at,"Atom1_ff_int_type");
   set_value(is_Atom2_ff_type,    Atom2_ff_type,    at,"Atom2_ff_type");
   set_value(is_Atom2_ff_int_type,Atom2_ff_int_type,at,"Atom2_ff_int_type");
   set_value(is_Atom3_ff_type,    Atom3_ff_type,    at,"Atom3_ff_type");
   set_value(is_Atom3_ff_int_type,Atom3_ff_int_type,at,"Atom3_ff_int_type");
   set_value(is_Atom4_ff_type,    Atom4_ff_type,    at,"Atom4_ff_type");
   set_value(is_Atom4_ff_int_type,Atom4_ff_int_type,at,"Atom4_ff_int_type");
   set_value(is_Dihedral_type_index,Dihedral_type_index,at,"Dihedral_type_index");


   set_value(is_Dihedral_vphi,    Dihedral_vphi,    at, "Dihedral_vphi");
   set_value(is_Dihedral_vphi1,   Dihedral_vphi1,   at, "Dihedral_vphi1");
   set_value(is_Dihedral_vphi2,   Dihedral_vphi2,   at, "Dihedral_vphi2");
   set_value(is_Dihedral_vphi3,   Dihedral_vphi3,   at, "Dihedral_vphi3");
   set_value(is_Dihedral_phase,   Dihedral_phase,   at, "Dihedral_phase");
   set_value(is_Dihedral_mult,    Dihedral_mult,    at, "Dihedral_mult");

   return 1;

}

Dihedral_Record& Dihedral_Record::operator=(const Dihedral_Record& at){

   is_Atom1_ff_type = 0;
   is_Atom2_ff_type = 0;
   is_Atom3_ff_type = 0;
   is_Atom4_ff_type = 0;

   is_Atom1_ff_int_type = 0;
   is_Atom2_ff_int_type = 0;
   is_Atom3_ff_int_type = 0;
   is_Atom4_ff_int_type = 0;

   is_Dihedral_type_index = 0;

   is_Dihedral_vphi  = 0;
   is_Dihedral_vphi1 = 0;
   is_Dihedral_vphi2 = 0;
   is_Dihedral_vphi3 = 0;
   is_Dihedral_phase = 0;
   is_Dihedral_mult  = 0;


   if(at.is_Atom1_ff_type)     {Atom1_ff_type     = at.Atom1_ff_type;     is_Atom1_ff_type = 1;}
   if(at.is_Atom1_ff_int_type) {Atom1_ff_int_type = at.Atom1_ff_int_type; is_Atom1_ff_int_type = 1;}
   if(at.is_Atom2_ff_type)     {Atom2_ff_type     = at.Atom2_ff_type;     is_Atom2_ff_type = 1;}
   if(at.is_Atom2_ff_int_type) {Atom2_ff_int_type = at.Atom2_ff_int_type; is_Atom2_ff_int_type = 1;}
   if(at.is_Atom3_ff_type)     {Atom3_ff_type     = at.Atom3_ff_type;     is_Atom3_ff_type = 1;}
   if(at.is_Atom3_ff_int_type) {Atom3_ff_int_type = at.Atom3_ff_int_type; is_Atom3_ff_int_type = 1;}
   if(at.is_Atom4_ff_type)     {Atom4_ff_type     = at.Atom4_ff_type;     is_Atom4_ff_type = 1;}
   if(at.is_Atom4_ff_int_type) {Atom4_ff_int_type = at.Atom4_ff_int_type; is_Atom4_ff_int_type = 1;}
   if(at.is_Dihedral_type_index){Dihedral_type_index = at.Dihedral_type_index; is_Dihedral_type_index = 1;}

   if(at.is_Dihedral_vphi)     { Dihedral_vphi      = at.Dihedral_vphi;      is_Dihedral_vphi  = 1; }
   if(at.is_Dihedral_vphi1)    { Dihedral_vphi1     = at.Dihedral_vphi1;     is_Dihedral_vphi1 = 1; }
   if(at.is_Dihedral_vphi2)    { Dihedral_vphi2     = at.Dihedral_vphi2;     is_Dihedral_vphi2 = 1; }
   if(at.is_Dihedral_vphi3)    { Dihedral_vphi3     = at.Dihedral_vphi3;     is_Dihedral_vphi3 = 1; }
   if(at.is_Dihedral_phase)    { Dihedral_phase     = at.Dihedral_phase;     is_Dihedral_phase = 1; }
   if(at.is_Dihedral_mult)     { Dihedral_mult      = at.Dihedral_mult;      is_Dihedral_mult  = 1; }   

   return *this;
}

void Dihedral_Record::merge(const Dihedral_Record& at){

   if(!is_Atom1_ff_type && at.is_Atom1_ff_type)            {Atom1_ff_type     = at.Atom1_ff_type;     is_Atom1_ff_type = 1;}
   if(!is_Atom1_ff_int_type && at.is_Atom1_ff_int_type)    {Atom1_ff_int_type = at.Atom1_ff_int_type; is_Atom1_ff_int_type = 1;}
   if(!is_Atom2_ff_type && at.is_Atom2_ff_type)            {Atom2_ff_type     = at.Atom2_ff_type;     is_Atom2_ff_type = 1;}
   if(!is_Atom2_ff_int_type && at.is_Atom2_ff_int_type)    {Atom2_ff_int_type = at.Atom2_ff_int_type; is_Atom2_ff_int_type = 1;}
   if(!is_Atom3_ff_type && at.is_Atom3_ff_type)            {Atom3_ff_type     = at.Atom3_ff_type;     is_Atom3_ff_type = 1;}
   if(!is_Atom3_ff_int_type && at.is_Atom3_ff_int_type)    {Atom3_ff_int_type = at.Atom3_ff_int_type; is_Atom3_ff_int_type = 1;}
   if(!is_Atom4_ff_type && at.is_Atom4_ff_type)            {Atom4_ff_type     = at.Atom4_ff_type;     is_Atom4_ff_type = 1;}
   if(!is_Atom4_ff_int_type && at.is_Atom4_ff_int_type)    {Atom4_ff_int_type = at.Atom4_ff_int_type; is_Atom4_ff_int_type = 1;}
   if(!is_Dihedral_type_index && at.is_Dihedral_type_index){Dihedral_type_index = at.Dihedral_type_index; is_Dihedral_type_index = 1;}

   if(!is_Dihedral_vphi && at.is_Dihedral_vphi)     { Dihedral_vphi      = at.Dihedral_vphi;      is_Dihedral_vphi  = 1; }
   if(!is_Dihedral_vphi1 && at.is_Dihedral_vphi1)   { Dihedral_vphi1     = at.Dihedral_vphi1;     is_Dihedral_vphi1 = 1; }
   if(!is_Dihedral_vphi2 && at.is_Dihedral_vphi2)   { Dihedral_vphi2     = at.Dihedral_vphi2;     is_Dihedral_vphi2 = 1; }
   if(!is_Dihedral_vphi3 && at.is_Dihedral_vphi3)   { Dihedral_vphi3     = at.Dihedral_vphi3;     is_Dihedral_vphi3 = 1; }
   if(!is_Dihedral_phase && at.is_Dihedral_phase)   { Dihedral_phase     = at.Dihedral_phase;     is_Dihedral_phase = 1; }
   if(!is_Dihedral_mult && at.is_Dihedral_mult)     { Dihedral_mult      = at.Dihedral_mult;      is_Dihedral_mult  = 1; }

}


int Dihedral_Record::show_info(){

   std::cout<<"Dihedral_Record properties:"<<std::endl;
   if(is_Atom1_ff_type)        {std::cout<<"Atom1_ff_type = "<<Atom1_ff_type<<std::endl;   }
   if(is_Atom1_ff_int_type)    {std::cout<<"Atom1_ff_int_type = "<<Atom1_ff_int_type<<std::endl;     }
   if(is_Atom2_ff_type)        {std::cout<<"Atom2_ff_type = "<<Atom2_ff_type<<std::endl;   }
   if(is_Atom2_ff_int_type)    {std::cout<<"Atom2_ff_int_type = "<<Atom2_ff_int_type<<std::endl;     }
   if(is_Atom3_ff_type)        {std::cout<<"Atom3_ff_type = "<<Atom3_ff_type<<std::endl;   }
   if(is_Atom3_ff_int_type)    {std::cout<<"Atom3_ff_int_type = "<<Atom3_ff_int_type<<std::endl;     }
   if(is_Atom4_ff_type)        {std::cout<<"Atom4_ff_type = "<<Atom4_ff_type<<std::endl;   }
   if(is_Atom4_ff_int_type)    {std::cout<<"Atom4_ff_int_type = "<<Atom4_ff_int_type<<std::endl;     }
   if(is_Dihedral_type_index)  {std::cout<<"Dihedral_type_index = "<<Dihedral_type_index<<std::endl; }

   if(is_Dihedral_vphi)        {std::cout<<"Dihedral_vphi = "<<Dihedral_vphi<<" kcal/mol"<<std::endl; }
   if(is_Dihedral_vphi1)       {std::cout<<"Dihedral_vphi1 = "<<Dihedral_vphi1<<" kcal/mol"<<std::endl; }
   if(is_Dihedral_vphi2)       {std::cout<<"Dihedral_vphi2 = "<<Dihedral_vphi2<<" kcal/mol"<<std::endl; }
   if(is_Dihedral_vphi3)       {std::cout<<"Dihedral_vphi3 = "<<Dihedral_vphi3<<" kcal/mol"<<std::endl; }
   if(is_Dihedral_phase)       {std::cout<<"Dihedral_phase = "<<Dihedral_phase<<" deg"<<std::endl; }
   if(is_Dihedral_mult)        {std::cout<<"Dihedral_mult = "<<Dihedral_mult<<" unitless"<<std::endl; }

   std::cout<<std::endl;

   return 1;

}

void Dihedral_Record::save(boost::property_tree::ptree& pt,std::string path){

  if(is_Atom1_ff_type){  ::save(pt,path+".Atom1_ff_type",Atom1_ff_type);    }
  if(is_Atom2_ff_type){  ::save(pt,path+".Atom2_ff_type",Atom2_ff_type);    }
  if(is_Atom3_ff_type){  ::save(pt,path+".Atom3_ff_type",Atom3_ff_type);    }
  if(is_Atom4_ff_type){  ::save(pt,path+".Atom4_ff_type",Atom4_ff_type);    }
  if(is_Atom1_ff_int_type){  ::save(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type);    }
  if(is_Atom2_ff_int_type){  ::save(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type);    }
  if(is_Atom3_ff_int_type){  ::save(pt,path+".Atom3_ff_int_type",Atom3_ff_int_type);    }
  if(is_Atom4_ff_int_type){  ::save(pt,path+".Atom4_ff_int_type",Atom4_ff_int_type);    }

  if(is_Dihedral_type_index){  ::save(pt,path+".Dihedral_type_index",Dihedral_type_index);    }
  if(is_Dihedral_vphi){  ::save(pt,path+".Dihedral_vphi",Dihedral_vphi);    }
  if(is_Dihedral_vphi1){  ::save(pt,path+".Dihedral_vphi1",Dihedral_vphi1);    }
  if(is_Dihedral_vphi2){  ::save(pt,path+".Dihedral_vphi2",Dihedral_vphi2);    }
  if(is_Dihedral_vphi3){  ::save(pt,path+".Dihedral_vphi3",Dihedral_vphi3);    }
  if(is_Dihedral_phase){  ::save(pt,path+".Dihedral_phase",Dihedral_phase);    }
  if(is_Dihedral_mult){  ::save(pt,path+".Dihedral_mult",Dihedral_mult);    }

}

void save(boost::property_tree::ptree& pt,std::string path,vector<Dihedral_Record>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".Dihedral_Record"+rt);
  }
}

void Dihedral_Record::load(boost::property_tree::ptree& pt,std::string path,int& status){
  int st;
  status = 0;

  ::load(pt,path+".Atom1_ff_type",Atom1_ff_type,is_Atom1_ff_type); if(is_Atom1_ff_type==1) { status=1;}
  ::load(pt,path+".Atom2_ff_type",Atom2_ff_type,is_Atom2_ff_type); if(is_Atom2_ff_type==1) { status=1;}
  ::load(pt,path+".Atom3_ff_type",Atom3_ff_type,is_Atom3_ff_type); if(is_Atom3_ff_type==1) { status=1;}
  ::load(pt,path+".Atom4_ff_type",Atom4_ff_type,is_Atom4_ff_type); if(is_Atom4_ff_type==1) { status=1;}
  ::load(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type,is_Atom1_ff_int_type); if(is_Atom1_ff_int_type==1) { status=1;}
  ::load(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type,is_Atom2_ff_int_type); if(is_Atom2_ff_int_type==1) { status=1;}
  ::load(pt,path+".Atom3_ff_int_type",Atom3_ff_int_type,is_Atom3_ff_int_type); if(is_Atom3_ff_int_type==1) { status=1;}
  ::load(pt,path+".Atom4_ff_int_type",Atom4_ff_int_type,is_Atom4_ff_int_type); if(is_Atom4_ff_int_type==1) { status=1;}

  ::load(pt,path+".Dihedral_type_index",Dihedral_type_index,is_Dihedral_type_index); if(is_Dihedral_type_index==1) { status=1;}
  ::load(pt,path+".Dihedral_vphi",Dihedral_vphi,is_Dihedral_vphi); if(is_Dihedral_vphi==1) { status=1;}
  ::load(pt,path+".Dihedral_vphi1",Dihedral_vphi1,is_Dihedral_vphi1); if(is_Dihedral_vphi1==1) { status=1;}
  ::load(pt,path+".Dihedral_vphi2",Dihedral_vphi2,is_Dihedral_vphi2); if(is_Dihedral_vphi2==1) { status=1;}
  ::load(pt,path+".Dihedral_vphi3",Dihedral_vphi3,is_Dihedral_vphi3); if(is_Dihedral_vphi3==1) { status=1;}
  ::load(pt,path+".Dihedral_phase",Dihedral_phase,is_Dihedral_phase); if(is_Dihedral_phase==1) { status=1;}
  ::load(pt,path+".Dihedral_mult",Dihedral_mult,is_Dihedral_mult); if(is_Dihedral_mult==1) { status=1;}

}

void load(boost::property_tree::ptree& pt,std::string path,vector<Dihedral_Record>& vt,int& status){
  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      Dihedral_Record x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}


//--------------- Fragment Record class members ---------------------------

int Fragment_Record::set(object at){

   set_value(is_Fragment_ff_type,    Fragment_ff_type,    at,"Fragment_ff_type");
   set_value(is_Fragment_ff_int_type,Fragment_ff_int_type,at,"Fragment_ff_int_type");

   set_value(is_Fragment_di,Fragment_di,at,"Fragment_di");
   set_value(is_Fragment_li,Fragment_li,at,"Fragment_li");
   set_value(is_Fragment_e0,Fragment_e0,at,"Fragment_e0");
   set_value(is_Fragment_rat,Fragment_rat,at,"Fragment_rat");
   set_value(is_Fragment_dw,Fragment_dw,at,"Fragment_dw");
   set_value(is_Fragment_mu,Fragment_mu,at,"Fragment_mu");
   set_value(is_Fragment_nu,Fragment_nu,at,"Fragment_nu");

  return 1;
}


int Fragment_Record::show_info(){

   std::cout<<"Fragment_Record properties:"<<std::endl;
   if(is_Fragment_ff_type)        {std::cout<<"Fragment_ff_type = "<<Fragment_ff_type<<std::endl;   }
   if(is_Fragment_ff_int_type)    {std::cout<<"Fragment_ff_int_type = "<<Fragment_ff_int_type<<std::endl;     }

   if(is_Fragment_di)    {std::cout<<"Fragment_di = "<<Fragment_di<<std::endl;     }
   if(is_Fragment_li)    {std::cout<<"Fragment_li = "<<Fragment_li<<std::endl;     }
   if(is_Fragment_e0)    {std::cout<<"Fragment_e0 = "<<Fragment_e0<<std::endl;     }
   if(is_Fragment_rat)    {std::cout<<"Fragment_rat = "<<Fragment_rat<<std::endl;     }
   if(is_Fragment_dw)    {std::cout<<"Fragment_dw = "<<Fragment_dw<<std::endl;     }
   if(is_Fragment_mu)    {std::cout<<"Fragment_mu = "<<Fragment_mu<<std::endl;     }
   if(is_Fragment_nu)    {std::cout<<"Fragment_nu = "<<Fragment_nu<<std::endl;     }

   std::cout<<std::endl;

   return 1;

}

Fragment_Record& Fragment_Record::operator=(const Fragment_Record& at){

   is_Fragment_ff_type     = 0;
   is_Fragment_ff_int_type = 0;
   is_Fragment_di = 0;
   is_Fragment_li = 0;
   is_Fragment_e0 = 0;
   is_Fragment_rat = 0;
   is_Fragment_dw = 0;
   is_Fragment_mu = 0;
   is_Fragment_nu = 0;


// This assignment is only working if the source object contain some data

   if(at.is_Fragment_ff_type)     {Fragment_ff_type     = at.Fragment_ff_type;  is_Fragment_ff_type = 1;}
   if(at.is_Fragment_ff_int_type) {Fragment_ff_int_type = at.Fragment_ff_int_type; is_Fragment_ff_int_type = 1;}

   if(at.is_Fragment_di) {Fragment_di = at.Fragment_di; is_Fragment_di = 1;}
   if(at.is_Fragment_li) {Fragment_li = at.Fragment_li; is_Fragment_li = 1;}
   if(at.is_Fragment_e0) {Fragment_e0 = at.Fragment_e0; is_Fragment_e0 = 1;}
   if(at.is_Fragment_rat) {Fragment_rat = at.Fragment_rat; is_Fragment_rat = 1;}
   if(at.is_Fragment_dw) {Fragment_dw = at.Fragment_dw; is_Fragment_dw = 1;}
   if(at.is_Fragment_mu) {Fragment_mu = at.Fragment_mu; is_Fragment_mu = 1;}
   if(at.is_Fragment_nu) {Fragment_nu = at.Fragment_nu; is_Fragment_nu = 1;}

   return *this;

}

void Fragment_Record::merge(const Fragment_Record& at){

   if(!is_Fragment_ff_type && at.is_Fragment_ff_type)         {Fragment_ff_type     = at.Fragment_ff_type;  is_Fragment_ff_type = 1;}
   if(!is_Fragment_ff_int_type && at.is_Fragment_ff_int_type) {Fragment_ff_int_type = at.Fragment_ff_int_type; is_Fragment_ff_int_type = 1;}

   if(!is_Fragment_di && at.is_Fragment_di) {Fragment_di = at.Fragment_di; is_Fragment_di = 1;}
   if(!is_Fragment_li && at.is_Fragment_li) {Fragment_li = at.Fragment_li; is_Fragment_li = 1;}
   if(!is_Fragment_e0 && at.is_Fragment_e0) {Fragment_e0 = at.Fragment_e0; is_Fragment_e0 = 1;}
   if(!is_Fragment_rat && at.is_Fragment_rat) {Fragment_rat = at.Fragment_rat; is_Fragment_rat = 1;}
   if(!is_Fragment_dw && at.is_Fragment_dw) {Fragment_dw = at.Fragment_dw; is_Fragment_dw = 1;}
   if(!is_Fragment_mu && at.is_Fragment_mu) {Fragment_mu = at.Fragment_mu; is_Fragment_mu = 1;}
   if(!is_Fragment_nu && at.is_Fragment_nu) {Fragment_nu = at.Fragment_nu; is_Fragment_nu = 1;}

}

void Fragment_Record::save(boost::property_tree::ptree& pt,std::string path){

  if(is_Fragment_ff_type){  ::save(pt,path+".Fragment_ff_type",Fragment_ff_type);    }
  if(is_Fragment_ff_int_type){  ::save(pt,path+".Fragment_ff_int_type",Fragment_ff_int_type);    }

  if(is_Fragment_di){  ::save(pt,path+".Fragment_di",Fragment_di);    }
  if(is_Fragment_li){  ::save(pt,path+".Fragment_li",Fragment_li);    }
  if(is_Fragment_e0){  ::save(pt,path+".Fragment_e0",Fragment_e0);    }
  if(is_Fragment_rat){  ::save(pt,path+".Fragment_rat",Fragment_rat);    }
  if(is_Fragment_dw){  ::save(pt,path+".Fragment_dw",Fragment_dw);    }
  if(is_Fragment_mu){  ::save(pt,path+".Fragment_mu",Fragment_mu);    }
  if(is_Fragment_nu){  ::save(pt,path+".Fragment_nu",Fragment_nu);    }

}

void save(boost::property_tree::ptree& pt,std::string path,vector<Fragment_Record>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    vt[i].save(pt,path+".Fragment_Record"+rt);
  }
}


void Fragment_Record::load(boost::property_tree::ptree& pt,std::string path,int& status){
  int st;
  status = 0;

  ::load(pt,path+".Fragment_ff_type",Fragment_ff_type,is_Fragment_ff_type); if(is_Fragment_ff_type==1) { status=1;}
  ::load(pt,path+".Fragment_ff_int_type",Fragment_ff_int_type,is_Fragment_ff_int_type); if(is_Fragment_ff_int_type==1) { status=1;}

  ::load(pt,path+".Fragment_di",Fragment_di,is_Fragment_di); if(is_Fragment_di==1) { status=1;}
  ::load(pt,path+".Fragment_li",Fragment_li,is_Fragment_li); if(is_Fragment_li==1) { status=1;}
  ::load(pt,path+".Fragment_e0",Fragment_e0,is_Fragment_e0); if(is_Fragment_e0==1) { status=1;}
  ::load(pt,path+".Fragment_rat",Fragment_rat,is_Fragment_rat); if(is_Fragment_rat==1) { status=1;}
  ::load(pt,path+".Fragment_dw",Fragment_dw,is_Fragment_dw); if(is_Fragment_dw==1) { status=1;}
  ::load(pt,path+".Fragment_mu",Fragment_mu,is_Fragment_mu); if(is_Fragment_mu==1) { status=1;}
  ::load(pt,path+".Fragment_nu",Fragment_nu,is_Fragment_nu); if(is_Fragment_nu==1) { status=1;}

}

void load(boost::property_tree::ptree& pt,std::string path,vector<Fragment_Record>& vt,int& status){
  int st;
  status = 0;
  try{
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(path)){
      Fragment_Record x; x.load(pt,path+"."+v.first,st); if(st==1){ vt.push_back(x); status = 1; }
    }
  }catch(std::exception& e){ }
}

/*
void ForceField::set(object at){

   set_value(is_ForceField_Name,       ForceField_Name,      at,"ForceField_Name");
   set_value(is_sigma_comb_rule,       sigma_comb_rule,      at,"sigma_comb_rule");
   set_value(is_epsilon_comb_rule,     epsilon_comb_rule,    at,"epsilon_comb_rule");

   set_value(is_system_pbc,            system_pbc,           at,"system_pbc");
   set_value(is_pbc_degree_x,          pbc_degree_x,         at,"pbc_degree_x");
   set_value(is_pbc_degree_y,          pbc_degree_y,         at,"pbc_degree_y");
   set_value(is_pbc_degree_z,          pbc_degree_z,         at,"pbc_degree_z");
   set_value(is_reciprocal_degree_x,   reciprocal_degree_x,  at,"reciprocal_degree_x");
   set_value(is_reciprocal_degree_y,   reciprocal_degree_y,  at,"reciprocal_degree_y");
   set_value(is_reciprocal_degree_z,   reciprocal_degree_z,  at,"reciprocal_degree_z");
   set_value(is_R_vdw_on,              R_vdw_on,             at,"R_vdw_on");
   set_value(is_R_vdw_off,             R_vdw_off,            at,"R_vdw_off");
   set_value(is_R_elec_on,             R_elec_on,            at,"R_elec_on");
   set_value(is_R_elec_off,            R_elec_off,           at,"R_elec_off");
   set_value(is_R_vlist,               R_vlist,              at,"R_vlist");
   set_value(is_elec_etha,             R_elec_etha,          at,"R_elec_etha");

   set_value(is_bond_functional,       bond_functional,      at,"bond_functional");
   set_value(is_angle_functional,      angle_functional,     at,"angle_functional");
   set_value(is_dihedral_functional,   dihedral_functional,  at,"dihedral_functional");
   set_value(is_oop_functional,        oop_functional,       at,"oop_functional");
   set_value(is_vdw_functional,        vdw_functional,       at,"vdw_functional");
   set_value(is_elec_functional,       elec_functional,      at,"elec_functional");

   set_value(is_vdw_scale12,           vdw_scale12,          at,"vdw_scale12");
   set_value(is_vdw_scale13,           vdw_scale13,          at,"vdw_scale13");
   set_value(is_vdw_scale14,           vdw_scale14,          at,"vdw_scale14");
   set_value(is_elec_scale12,          elec_scale12,         at,"elec_scale12");
   set_value(is_elec_scale13,          elec_scale13,         at,"elec_scale13");
   set_value(is_elec_scale14,          elec_scale14,         at,"elec_scale14");

}
*/

void ForceField::set(boost::python::dict d){
  extract_dictionary(d);
}

void ForceField::extract_dictionary(boost::python::dict d){
  std::string key;
  for(int i=0;i<len(d.values());i++){
    key = extract<std::string>(d.keys()[i]);

    if(key=="ForceField_Name"){ ForceField_Name = extract<std::string>(d.values()[i]); is_ForceField_Name = 1; }
    else if(key=="bond_functional"){ bond_functional = extract<std::string>(d.values()[i]); is_bond_functional = 1;}
    else if(key=="angle_functional"){ angle_functional = extract<std::string>(d.values()[i]); is_angle_functional = 1;}
    else if(key=="dihedral_functional"){ dihedral_functional = extract<std::string>(d.values()[i]); is_dihedral_functional = 1;}
    else if(key=="oop_functional"){ oop_functional = extract<std::string>(d.values()[i]); is_oop_functional = 1;}
    else if(key=="vdw_functional"){ vdw_functional = extract<std::string>(d.values()[i]); is_vdw_functional = 1;}
    else if(key=="elec_functional"){ elec_functional = extract<std::string>(d.values()[i]); is_elec_functional = 1;}
    else if(key=="mb_functional"){ mb_functional = extract<std::string>(d.values()[i]); is_mb_functional = 1;}
    else if(key=="cg_functional"){ cg_functional = extract<std::string>(d.values()[i]); is_cg_functional = 1;}
    else if(key=="mb_excl_functional"){ mb_excl_functional = extract<std::string>(d.values()[i]); is_mb_excl_functional = 1;}


//    else if(key=="stress_opt"){ stress_opt = extract<std::string>(d.values()[i]); is_stress_opt = 1;}

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

void ForceField::init_variables(){
  is_ForceField_Name               = 0;

  is_system_pbc                    = 0;
  is_pbc_degree_x                  = 0;
  is_pbc_degree_y                  = 0;
  is_pbc_degree_z                  = 0;
  is_reciprocal_degree_x           = 0;
  is_reciprocal_degree_y           = 0;
  is_reciprocal_degree_z           = 0;
  is_R_vdw_on                      = 0;
  is_R_vdw_off                     = 0;
  is_R_elec_on                     = 0;
  is_R_elec_off                    = 0;
  is_R_vlist                       = 0;
  is_elec_etha                     = 0;

  is_bond_functional               = 0;
  is_angle_functional              = 0;
  is_dihedral_functional           = 0;
  is_oop_functional                = 0;
  is_vdw_functional                = 0;
  is_elec_functional               = 0;
  is_mb_functional                 = 0;
  is_cg_functional                 = 0;
  is_mb_excl_functional            = 0;


//  stress_opt = "fr";               is_stress_opt = 1;

  is_sigma_comb_rule    = 0;
  is_epsilon_comb_rule  = 0;
  is_vdw_scale12        = 0;
  is_vdw_scale13        = 0;
  is_vdw_scale14        = 0;
  is_elec_scale12       = 0;
  is_elec_scale13       = 0;
  is_elec_scale14       = 0;

}

void ForceField::copy_content(const ForceField& ff){

  if(ff.is_ForceField_Name){ ForceField_Name = ff.ForceField_Name; is_ForceField_Name = 1; }

  if(ff.is_bond_functional){ bond_functional = ff.bond_functional; is_bond_functional = 1; }
  if(ff.is_angle_functional){ angle_functional = ff.angle_functional; is_angle_functional = 1; }
  if(ff.is_dihedral_functional){ dihedral_functional = ff.dihedral_functional; is_dihedral_functional = 1; }
  if(ff.is_oop_functional){ oop_functional = ff.oop_functional; is_oop_functional = 1; }
  if(ff.is_vdw_functional){ vdw_functional = ff.vdw_functional; is_vdw_functional = 1; }
  if(ff.is_elec_functional){ elec_functional = ff.elec_functional; is_elec_functional = 1; }
  if(ff.is_mb_functional){ mb_functional = ff.mb_functional; is_mb_functional = 1; }
  if(ff.is_cg_functional){ cg_functional = ff.cg_functional; is_cg_functional = 1; }
  if(ff.is_mb_excl_functional){ mb_excl_functional = ff.mb_excl_functional; is_mb_excl_functional = 1; }


//  if(ff.is_stress_opt){ stress_opt = ff.stress_opt; is_stress_opt = 1; }

  if(ff.is_system_pbc){ system_pbc = ff.system_pbc; is_system_pbc = 1; }
  if(ff.is_pbc_degree_x){ pbc_degree_x = ff.pbc_degree_x; is_pbc_degree_x = 1; }
  if(ff.is_pbc_degree_y){ pbc_degree_y = ff.pbc_degree_y; is_pbc_degree_y = 1; }
  if(ff.is_pbc_degree_z){ pbc_degree_z = ff.pbc_degree_z; is_pbc_degree_z = 1; }
  if(ff.is_reciprocal_degree_x){ reciprocal_degree_x = ff.reciprocal_degree_x; is_reciprocal_degree_x = 1; }
  if(ff.is_reciprocal_degree_y){ reciprocal_degree_y = ff.reciprocal_degree_y; is_reciprocal_degree_y = 1; }
  if(ff.is_reciprocal_degree_z){ reciprocal_degree_z = ff.reciprocal_degree_z; is_reciprocal_degree_z = 1; }
  if(ff.is_R_vdw_on){ R_vdw_on = ff.R_vdw_on; is_R_vdw_on = 1; }
  if(ff.is_R_vdw_off){ R_vdw_off = ff.R_vdw_off; is_R_vdw_off = 1; }
  if(ff.is_R_elec_on){ R_elec_on = ff.R_elec_on; is_R_elec_on = 1; }
  if(ff.is_R_elec_off){ R_elec_off = ff.R_elec_off; is_R_elec_off = 1; }
  if(ff.is_R_vlist){ R_vlist = ff.R_vlist; is_R_vlist = 1; }
  if(ff.is_elec_etha){ elec_etha = ff.elec_etha; is_elec_etha = 1; }

  if(ff.is_sigma_comb_rule){ sigma_comb_rule = ff.sigma_comb_rule; is_sigma_comb_rule = 1; }
  if(ff.is_epsilon_comb_rule){ epsilon_comb_rule = ff.epsilon_comb_rule; is_epsilon_comb_rule = 1; }
  if(ff.is_vdw_scale12){ vdw_scale12 = ff.vdw_scale12; is_vdw_scale12 = 1; }
  if(ff.is_vdw_scale13){ vdw_scale13 = ff.vdw_scale13; is_vdw_scale13 = 1; }
  if(ff.is_vdw_scale14){ vdw_scale14 = ff.vdw_scale14; is_vdw_scale14 = 1; }
  if(ff.is_elec_scale12){ elec_scale12 = ff.elec_scale12; is_elec_scale12 = 1; }
  if(ff.is_elec_scale13){ elec_scale13 = ff.elec_scale13; is_elec_scale13 = 1; }
  if(ff.is_elec_scale14){ elec_scale14 = ff.elec_scale14; is_elec_scale14 = 1; }

  Atom_Records = ff.Atom_Records;
  Bond_Records = ff.Bond_Records;
  Angle_Records = ff.Angle_Records;
  Dihedral_Records = ff.Dihedral_Records;
  Improper_Records = ff.Improper_Records;
  Fragment_Records = ff.Fragment_Records;
}

ForceField::ForceField(){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
}

ForceField::ForceField(boost::python::dict d){
  /****************
     Constructor
  ******************/
  // Initialize variables to default values
  init_variables();
  extract_dictionary(d);
}

ForceField::ForceField(const ForceField& ff){
  /********************
    Copy constructor
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(ff);
}

ForceField& ForceField::operator=(const ForceField& ff){
  /********************
    Assignment operator
  *********************/
  // Initialize variables to default values
  init_variables();
  // Copy content of th object which is defined
  copy_content(ff);
  return *this;
}

ForceField::~ForceField(){
  if(Atom_Records.size()>0) { Atom_Records.clear(); }
  if(Bond_Records.size()>0) { Bond_Records.clear(); }
  if(Angle_Records.size()>0) { Angle_Records.clear(); }
  if(Dihedral_Records.size()>0) { Dihedral_Records.clear(); }
  if(Improper_Records.size()>0) { Improper_Records.clear(); }
  if(Fragment_Records.size()>0) { Fragment_Records.clear(); }
}

void ForceField::show_info(){

  std::cout<<"ForceField properties:"<<std::endl;
  if(is_ForceField_Name)      {std::cout<<"ForceField_Name = "<<ForceField_Name<<std::endl;   }
  if(is_sigma_comb_rule)      {std::cout<<"sigma_comb_rule = "<<sigma_comb_rule<<std::endl;   }
  if(is_epsilon_comb_rule)    {std::cout<<"epsilon_comb_rule = "<<epsilon_comb_rule<<std::endl;   }

  if(is_system_pbc)           {std::cout<<"system_pbc = "<<system_pbc<<std::endl; }
  if(is_pbc_degree_x)         {std::cout<<"pbc_degree_x = "<<pbc_degree_x<<std::endl; }
  if(is_pbc_degree_y)         {std::cout<<"pbc_degree_y = "<<pbc_degree_y<<std::endl; }
  if(is_pbc_degree_z)         {std::cout<<"pbc_degree_z = "<<pbc_degree_z<<std::endl; }
  if(is_reciprocal_degree_x)  {std::cout<<"reciprocal_degree_x = "<<reciprocal_degree_x<<std::endl; }
  if(is_reciprocal_degree_y)  {std::cout<<"reciprocal_degree_y = "<<reciprocal_degree_y<<std::endl; }
  if(is_reciprocal_degree_z)  {std::cout<<"reciprocal_degree_z = "<<reciprocal_degree_z<<std::endl; }
  if(is_R_vdw_on)             {std::cout<<"R_vdw_on = "<<R_vdw_on<<std::endl; }
  if(is_R_vdw_off)            {std::cout<<"R_vdw_off = "<<R_vdw_off<<std::endl; }
  if(is_R_elec_on)            {std::cout<<"R_elec_on = "<<R_elec_on<<std::endl; }
  if(is_R_elec_off)           {std::cout<<"R_elec_off = "<<R_elec_off<<std::endl; }
  if(is_R_vlist)              {std::cout<<"R_vlist = "<<R_vlist<<std::endl; }
  if(is_elec_etha)            {std::cout<<"elec_etha = "<<elec_etha<<std::endl; }

  if(is_bond_functional)      {std::cout<<"bond_functional = "<<bond_functional<<std::endl; }
  if(is_angle_functional)     {std::cout<<"angle_functional = "<<angle_functional<<std::endl; }
  if(is_dihedral_functional)  {std::cout<<"dihedral_functional = "<<dihedral_functional<<std::endl; }
  if(is_oop_functional)       {std::cout<<"oop_functional = "<<oop_functional<<std::endl; }
  if(is_vdw_functional)       {std::cout<<"vdw_functional = "<<vdw_functional<<std::endl; }
  if(is_elec_functional)      {std::cout<<"elec_functional = "<<elec_functional<<std::endl; }
  if(is_mb_functional)        {std::cout<<"mb_functional = "<<mb_functional<<std::endl; }
  if(is_cg_functional)        {std::cout<<"cg_functional = "<<cg_functional<<std::endl; }
  if(is_mb_excl_functional)   {std::cout<<"mb_excl_functional = "<<mb_excl_functional<<std::endl; }

//  if(is_stress_opt)           {std::cout<<"stress_opt = "<<stress_opt<<std::endl; }

  if(is_vdw_scale12)          {std::cout<<"vdw_scale12 = "<<vdw_scale12<<std::endl; }
  if(is_vdw_scale13)          {std::cout<<"vdw_scale13 = "<<vdw_scale13<<std::endl; }
  if(is_vdw_scale14)          {std::cout<<"vdw_scale14 = "<<vdw_scale14<<std::endl; }
  if(is_elec_scale12)         {std::cout<<"elec_scale12 = "<<elec_scale12<<std::endl; }
  if(is_elec_scale13)         {std::cout<<"elec_scale13 = "<<elec_scale13<<std::endl; }
  if(is_elec_scale14)         {std::cout<<"elec_scale14 = "<<elec_scale14<<std::endl; }

}


void ForceField::save(boost::property_tree::ptree& pt,std::string path){

  if(is_ForceField_Name){  ::save(pt,path+".ForceField_Name",ForceField_Name);    }
  if(is_bond_functional){  ::save(pt,path+".bond_functional",bond_functional);    }
  if(is_angle_functional){  ::save(pt,path+".angle_functional",angle_functional);    }
  if(is_dihedral_functional){  ::save(pt,path+".dihedral_functional",dihedral_functional);    }
  if(is_oop_functional){  ::save(pt,path+".oop_functional",oop_functional);    }
  if(is_vdw_functional){  ::save(pt,path+".vdw_functional",vdw_functional);    }
  if(is_elec_functional){  ::save(pt,path+".elec_functional",elec_functional);    }
  if(is_mb_functional){  ::save(pt,path+".mb_functional",mb_functional);    }
  if(is_cg_functional){  ::save(pt,path+".cg_functional",cg_functional);    }
  if(is_mb_excl_functional){  ::save(pt,path+".mb_excl_functional",mb_excl_functional);    }

  if(is_system_pbc){  ::save(pt,path+".system_pbc",system_pbc);    }
  if(is_pbc_degree_x){  ::save(pt,path+".pbc_degree_x",pbc_degree_x);    }
  if(is_pbc_degree_y){  ::save(pt,path+".pbc_degree_y",pbc_degree_y);    }
  if(is_pbc_degree_z){  ::save(pt,path+".pbc_degree_z",pbc_degree_z);    }
  if(is_reciprocal_degree_x){  ::save(pt,path+".reciprocal_degree_x",reciprocal_degree_x);    }
  if(is_reciprocal_degree_y){  ::save(pt,path+".reciprocal_degree_y",reciprocal_degree_y);    }
  if(is_reciprocal_degree_z){  ::save(pt,path+".reciprocal_degree_z",reciprocal_degree_z);    }
  if(is_R_vdw_on){  ::save(pt,path+".R_vdw_on",R_vdw_on);    }
  if(is_R_vdw_off){  ::save(pt,path+".R_vdw_off",R_vdw_off);    }
  if(is_R_elec_on){  ::save(pt,path+".R_elec_on",R_elec_on);    }
  if(is_R_elec_off){  ::save(pt,path+".R_elec_off",R_elec_off);    }
  if(is_R_vlist){  ::save(pt,path+".R_vlist",R_vlist);    }
  if(is_elec_etha){  ::save(pt,path+".elec_etha",elec_etha);    }

  if(is_sigma_comb_rule){  ::save(pt,path+".sigma_comb_rule",sigma_comb_rule);    }
  if(is_epsilon_comb_rule){  ::save(pt,path+".epsilon_comb_rule",epsilon_comb_rule);    }
  if(is_vdw_scale12){  ::save(pt,path+".vdw_scale12",vdw_scale12);    }
  if(is_vdw_scale13){  ::save(pt,path+".vdw_scale13",vdw_scale13);    }
  if(is_vdw_scale14){  ::save(pt,path+".vdw_scale14",vdw_scale14);    }
  if(is_elec_scale12){  ::save(pt,path+".elec_scale12",elec_scale12);    }
  if(is_elec_scale13){  ::save(pt,path+".elec_scale13",elec_scale13);    }
  if(is_elec_scale14){  ::save(pt,path+".elec_scale14",elec_scale14);    }

  namespace here = libhamiltonian::libhamiltonian_atomistic::libhamiltonian_mm::libforcefield;


  if(Atom_Records.size()>0){  here::save(pt,path+".Atom_Records",Atom_Records);    }
  if(Bond_Records.size()>0){  here::save(pt,path+".Bond_Records",Bond_Records);    }
  if(Angle_Records.size()>0){  here::save(pt,path+".Angle_Records",Angle_Records);    }
  if(Dihedral_Records.size()>0){  here::save(pt,path+".Dihedral_Records",Dihedral_Records);    }
  if(Improper_Records.size()>0){  here::save(pt,path+".Improper_Records",Improper_Records);    }
  if(Fragment_Records.size()>0){  here::save(pt,path+".Fragment_Records",Fragment_Records);    }

}

void ForceField::load(boost::property_tree::ptree& pt,std::string path,int& status){
  int st;
  status = 0;

  ::load(pt,path+".ForceField_Name",ForceField_Name,is_ForceField_Name); if(is_ForceField_Name==1) { status=1;}
  ::load(pt,path+".bond_functional",bond_functional,is_bond_functional); if(is_bond_functional==1) { status=1;}
  ::load(pt,path+".angle_functional",angle_functional,is_angle_functional); if(is_angle_functional==1) { status=1;}
  ::load(pt,path+".dihedral_functional",dihedral_functional,is_dihedral_functional); if(is_dihedral_functional==1) { status=1;}
  ::load(pt,path+".oop_functional",oop_functional,is_oop_functional); if(is_oop_functional==1) { status=1;}
  ::load(pt,path+".vdw_functional",vdw_functional,is_vdw_functional); if(is_vdw_functional==1) { status=1;}
  ::load(pt,path+".elec_functional",elec_functional,is_elec_functional); if(is_elec_functional==1) { status=1;}
  ::load(pt,path+".mb_functional",mb_functional,is_mb_functional); if(is_mb_functional==1) { status=1;}
  ::load(pt,path+".cg_functional",cg_functional,is_cg_functional); if(is_cg_functional==1) { status=1;}
  ::load(pt,path+".mb_excl_functional",mb_excl_functional,is_mb_excl_functional); if(is_mb_excl_functional==1) { status=1;}


  ::load(pt,path+".system_pbc",system_pbc,is_system_pbc); if(is_system_pbc==1) { status=1;}
  ::load(pt,path+".pbc_degree_x",pbc_degree_x,is_pbc_degree_x); if(is_pbc_degree_x==1) { status=1;}
  ::load(pt,path+".pbc_degree_y",pbc_degree_y,is_pbc_degree_y); if(is_pbc_degree_y==1) { status=1;}
  ::load(pt,path+".pbc_degree_z",pbc_degree_z,is_pbc_degree_z); if(is_pbc_degree_z==1) { status=1;}
  ::load(pt,path+".reciprocal_degree_x",reciprocal_degree_x,is_reciprocal_degree_x); if(is_reciprocal_degree_x==1) { status=1;}
  ::load(pt,path+".reciprocal_degree_y",reciprocal_degree_y,is_reciprocal_degree_y); if(is_reciprocal_degree_y==1) { status=1;}
  ::load(pt,path+".reciprocal_degree_z",reciprocal_degree_z,is_reciprocal_degree_z); if(is_reciprocal_degree_z==1) { status=1;}
  ::load(pt,path+".R_vdw_on",R_vdw_on,is_R_vdw_on); if(is_R_vdw_on==1) { status=1;}
  ::load(pt,path+".R_vdw_off",R_vdw_off,is_R_vdw_off); if(is_R_vdw_off==1) { status=1;}
  ::load(pt,path+".R_elec_on",R_elec_on,is_R_elec_on); if(is_R_elec_on==1) { status=1;}
  ::load(pt,path+".R_elec_off",R_elec_off,is_R_elec_off); if(is_R_elec_off==1) { status=1;}
  ::load(pt,path+".R_vlist",R_vlist,is_R_vlist); if(is_R_vlist==1) { status=1;}
  ::load(pt,path+".elec_etha",elec_etha,is_elec_etha); if(is_elec_etha==1) { status=1;}
  ::load(pt,path+".sigma_comb_rule",sigma_comb_rule,is_sigma_comb_rule); if(is_sigma_comb_rule==1) { status=1;}
  ::load(pt,path+".epsilon_comb_rule",epsilon_comb_rule,is_epsilon_comb_rule); if(is_epsilon_comb_rule==1) { status=1;}
  ::load(pt,path+".vdw_scale12",vdw_scale12,is_vdw_scale12); if(is_vdw_scale12==1) { status=1;}
  ::load(pt,path+".vdw_scale13",vdw_scale13,is_vdw_scale13); if(is_vdw_scale13==1) { status=1;}
  ::load(pt,path+".vdw_scale14",vdw_scale14,is_vdw_scale14); if(is_vdw_scale14==1) { status=1;}
  ::load(pt,path+".elec_scale12",elec_scale12,is_elec_scale12); if(is_elec_scale12==1) { status=1;}
  ::load(pt,path+".elec_scale13",elec_scale13,is_elec_scale13); if(is_elec_scale13==1) { status=1;}
  ::load(pt,path+".elec_scale14",elec_scale14,is_elec_scale14); if(is_elec_scale14==1) { status=1;}

  namespace here = libhamiltonian::libhamiltonian_atomistic::libhamiltonian_mm::libforcefield;

  here::load(pt,path+".Atom_Records",Atom_Records,st); if(st==1) { status=1;}
  here::load(pt,path+".Bond_Records",Bond_Records,st); if(st==1) { status=1;}
  here::load(pt,path+".Angle_Records",Angle_Records,st); if(st==1) { status=1;}
  here::load(pt,path+".Dihedral_Records",Dihedral_Records,st); if(st==1) { status=1;}
  here::load(pt,path+".Improper_Records",Improper_Records,st); if(st==1) { status=1;}
  here::load(pt,path+".Fragment_Records",Fragment_Records,st); if(st==1) { status=1;}

}




int ForceField::Atom_Record_Index(int Atom_ff_int_type){
/****************************************************************
   This function searches for index of Atom_Records vector in which
   data about atom of force field type "Atom_ff_int_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Atom_Records.size();

   for(int i=0;i<sz;i++){
       if(Atom_Records[i].is_Atom_ff_int_type){
          if(Atom_Records[i].Atom_ff_int_type==Atom_ff_int_type){
             indx = i;
             break;
          }
       }
   }

   return indx;
}

int ForceField::Atom_Record_Index(std::string Atom_ff_type){
/****************************************************************
   This function searches for index of Atom_Records vector in which
   data about atom of force field type "Atom_ff_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Atom_Records.size();

   for(int i=0;i<sz;i++){
       if(Atom_Records[i].is_Atom_ff_type){
          if(Atom_Records[i].Atom_ff_type==Atom_ff_type){
             indx = i;
             break;
          }
       }
   }

   return indx;
}


int ForceField::Atom_Record_Index_by_Element(std::string Atom_element){
/****************************************************************
   This function searches for index of Atom_Records vector in which
   data about atom of element name "Atom_element" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Atom_Records.size();

   for(int i=0;i<sz;i++){
       if(Atom_Records[i].is_Atom_element){
          if(Atom_Records[i].Atom_element==Atom_element){
             indx = i;
             break;
          }
       }
   }

   return indx;
}


int ForceField::Atom_Record_Index_by_Element(std::string Atom_element,vector<int>& res){
/****************************************************************
   This function searches for index of Atom_Records vector in which
   data about atom of element name "Atom_element" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Atom_Records.size();
   if(res.size()>0) { res.clear(); }

   for(int i=0;i<sz;i++){
       if(Atom_Records[i].is_Atom_element){
          if(Atom_Records[i].Atom_element==Atom_element){
             indx = i;
             res.push_back(i);
          }
       }
   }

   return indx;
}

int ForceField::Atom_Record_Index_by_Element(int Atom_atomic_number,vector<int>& res){
/****************************************************************
   This function searches for index of Atom_Records vector in which
   data about atom of element number "Atom_atomic_number" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Atom_Records.size();
   if(res.size()>0) { res.clear(); }

   for(int i=0;i<sz;i++){
       if(Atom_Records[i].is_Atom_atomic_number){
          if(Atom_Records[i].Atom_atomic_number==Atom_atomic_number){
             indx = i;
             res.push_back(i);          
          }
       }
   }

   return indx;
}

int ForceField::Add_Atom_Record(Atom_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Atom_Record rec in array Atom_Records. If this atom type (atom record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the atom record rec is not valid - will not be added.
*********************************************************************/
   int res = 1;
   int sz = Atom_Records.size();

   int int_types = (rec.is_Atom_ff_int_type);
   int sym_types = (rec.is_Atom_ff_type);
   int elt_types = (rec.is_Atom_atomic_number);

   vector<int> at_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_type
   last_indx = Atom_Record_Index(rec.Atom_ff_int_type);
   at_indxs.push_back(last_indx);
   }
   else if(sym_types){
   // Search on the basis of symbolic type
   last_indx = Atom_Record_Index(rec.Atom_ff_type);
   at_indxs.push_back(last_indx);
   }
   else if(elt_types){
   // Search on the basis of atomic numbers
   last_indx = Atom_Record_Index_by_Element(rec.Atom_atomic_number,at_indxs);
   }
   else{
   // Record will not be added to the table of Angle_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2;
   }

   if(last_indx==-2){ // In this case the record was not added to Atom_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The atom record is valid but does not exist in Atom_Records array
                           // so it will be added completely
       Atom_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar atom record is already defined - try to add new properties

       for(int i=0;i<at_indxs.size();i++){
           Atom_Records[at_indxs[i]].merge(rec);
       }
       res = 0;
   }

   return res;

}

/*
int ForceField::Add_Atom_Record(Atom_Record rec){
*********************************************************************
   This is user-interface function. It checks for existance of the 
   Atom_Record rec in array Atom_Records. If this atom type (atom record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array (we did not add it this time). Returns -1
   if there is an error.
*********************************************************************
   int res = 1;   
   int sz = Atom_Records.size();
   int indx = sz;
   Atom_Record record; 
   record = rec;  

   if(rec.is_Atom_ff_type){
   
      for(int i=0;i<sz;i++){
       
         if(Atom_Records[i].is_Atom_ff_type){ 
            // This type already exist
            if(Atom_Records[i].Atom_ff_type==rec.Atom_ff_type){
                //res = 0;  - this is old version: means do not modify existing data
                res = 1; // new - add new data to existing record
                record = Atom_Records[i];
                indx = i; break; 
            }            
         }        
      }// for i

   }// is_Atom_ff_type == 1 (defined)
   else{
    std::cout<<"To add atom type record to the force field its Atom_ff_type should be defined"<<std::endl;
    res = -1;
   }

   if(res==1){
      // We will add this atom type to the force field 
      if(!rec.is_Atom_ff_int_type){
          // If Atom_ff_int_type not defined we define it
          // to be an index of this record in array
          rec.is_Atom_ff_int_type = 1;
          rec.Atom_ff_int_type = indx;
      }
      Atom_Records.push_back(record);
   }
   
   return res;

}
*/



int ForceField::Bond_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type){
/****************************************************************
   This function searches for index of Bond_Records vector in which
   data about bond formed by 2 atoms of of force field types
   "Atom1_ff_int_type" and "Atom2_ff_int_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Bond_Records.size();

   for(int i=0;i<sz;i++){
       if(Bond_Records[i].is_Atom1_ff_int_type && Bond_Records[i].is_Atom2_ff_int_type){
          if(((Bond_Records[i].Atom1_ff_int_type==Atom1_ff_int_type)&&(Bond_Records[i].Atom2_ff_int_type==Atom2_ff_int_type) )|| 
             ((Bond_Records[i].Atom1_ff_int_type==Atom2_ff_int_type)&&(Bond_Records[i].Atom2_ff_int_type==Atom1_ff_int_type) )
            )
          {
             indx = i;
             break;
          }
       }
   }

   return indx;
}


int ForceField::Bond_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type,int order){
/****************************************************************
   This function searches for index of Bond_Records vector in which
   data about bond formed by 2 atoms of of force field types
   "Atom1_ff_int_type" and "Atom2_ff_int_type" are stored
   Returns -1 if such index has not been found
   order parameter distinguish bonds 1-2 and 2-1
*****************************************************************/

   int indx = -1;
   int sz   = Bond_Records.size();

   int cmpr11,cmpr12,cmpr21,cmpr22,res;
   for(int i=0;i<sz;i++){
       if(Bond_Records[i].is_Atom1_ff_int_type && Bond_Records[i].is_Atom2_ff_int_type){
          cmpr11 = (Bond_Records[i].Atom1_ff_int_type==Atom1_ff_int_type);
          cmpr22 = (Bond_Records[i].Atom2_ff_int_type==Atom2_ff_int_type);
          cmpr12 = (Bond_Records[i].Atom1_ff_int_type==Atom2_ff_int_type);
          cmpr21 = (Bond_Records[i].Atom2_ff_int_type==Atom1_ff_int_type);

          res = 0; 
          if(order==1){  res = (cmpr11 && cmpr22);       }
          else{          res = ( (cmpr11 && cmpr22) || (cmpr12 && cmpr21) );   }


          if(res){ indx = i;  break;   }

       }// if is_Atom1_ff_int_type && is_Atom2_ff_int_type

   }// for i

   return indx;
}


int ForceField::Bond_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Bond_type_index, vector<int>& res){
/****************************************************************
   This function searches for indexes of Bond_Records vector in which
   data about bond formed by 2 atoms of of force field types
   "Atom1_ff_int_type" and "Atom2_ff_int_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of bond_records which match the search parameters
   The bond record will be chosen on the basis of additional comparison
   the Bond_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Bond_Records.size();
   int cmpr11,cmpr12,cmpr21,cmpr22,cmpr;

   for(int i=0;i<sz;i++){

       cmpr11 = cmpr12 = cmpr21 = cmpr22 = 0;
       if(Bond_Records[i].is_Atom1_ff_int_type && Bond_Records[i].is_Atom2_ff_int_type){

          cmpr11 = (Bond_Records[i].Atom1_ff_int_type == Atom1_ff_int_type);
          cmpr22 = (Bond_Records[i].Atom2_ff_int_type == Atom2_ff_int_type);
          cmpr12 = (Bond_Records[i].Atom1_ff_int_type == Atom2_ff_int_type);
          cmpr21 = (Bond_Records[i].Atom2_ff_int_type == Atom1_ff_int_type);

       }// if both Atom1_ff_int_type and Atom2_ff_int_type are defined    
       

       if(Bond_Records[i].is_Bond_type_index){
          
          cmpr = (Bond_Records[i].Bond_type_index == Bond_type_index);
          
       }// if Bond_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if((cmpr11&&cmpr22&&cmpr)||(cmpr12&&cmpr21&&cmpr)){
          indx = i; 
          res.push_back(i);
       }
   }// for i
   return indx;
}


int ForceField::Bond_Record_Index_by_Element(int Atom1_atomic_number, int Atom2_atomic_number, int Bond_type_index, vector<int>& res){
/****************************************************************
   This function searches for indexes of Bond_Records vector in which
   data about bond formed by 2 atoms of of atomic numbers
   "Atom1_atomic_number" and "Atom2_atomic_number" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of bond_records which match the search parameters
   The bond record will be chosen on the basis of additional comparison
   the Bond_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Bond_Records.size();
   int cmpr11,cmpr12,cmpr21,cmpr22,cmpr;

   for(int i=0;i<sz;i++){

       cmpr11 = cmpr12 = cmpr21 = cmpr22 = 0;
       if(Bond_Records[i].is_Atom1_atomic_number && Bond_Records[i].is_Atom2_atomic_number){

          cmpr11 = (Bond_Records[i].Atom1_atomic_number == Atom1_atomic_number);
          cmpr22 = (Bond_Records[i].Atom2_atomic_number == Atom2_atomic_number);
          cmpr12 = (Bond_Records[i].Atom1_atomic_number == Atom2_atomic_number);
          cmpr21 = (Bond_Records[i].Atom2_atomic_number == Atom1_atomic_number);

       }// if both Atom1_ff_int_type and Atom2_ff_int_type are defined


       if(Bond_Records[i].is_Bond_type_index){

          cmpr = (Bond_Records[i].Bond_type_index == Bond_type_index);

       }// if Bond_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if((cmpr11&&cmpr22&&cmpr)||(cmpr12&&cmpr21&&cmpr)){
          indx = i;
          res.push_back(i);
       }
   }// for i
   return indx;
}

int ForceField::Bond_Record_Index_by_Element(std::string Atom1_element, std::string Atom2_element, int Bond_type_index, vector<int>& res){
/****************************************************************
   This function searches for indexes of Bond_Records vector in which
   data about bond formed by 2 atoms of elements
   "Atom1_element" and "Atom2_element" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of bond_records which match the search parameters
   The bond record will be chosen on the basis of additional comparison
   the Bond_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Bond_Records.size();
   int cmpr11,cmpr12,cmpr21,cmpr22,cmpr;

   for(int i=0;i<sz;i++){

       cmpr11 = cmpr12 = cmpr21 = cmpr22 = 0;
       if(Bond_Records[i].is_Atom1_element && Bond_Records[i].is_Atom2_element){

          cmpr11 = (Bond_Records[i].Atom1_element == Atom1_element);
          cmpr22 = (Bond_Records[i].Atom2_element == Atom2_element);
          cmpr12 = (Bond_Records[i].Atom1_element == Atom2_element);
          cmpr21 = (Bond_Records[i].Atom2_element == Atom1_element);

       }// if both Atom1_ff_int_type and Atom2_ff_int_type are defined


       if(Bond_Records[i].is_Bond_type_index){

          cmpr = (Bond_Records[i].Bond_type_index == Bond_type_index);

       }// if Bond_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if((cmpr11&&cmpr22&&cmpr)||(cmpr12&&cmpr21&&cmpr)){
          indx = i;
          res.push_back(i);
       }
   }// for i
   return indx;
}

int ForceField::Bond_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type){
/****************************************************************
   This function searches for indexes of Bond_Records vector in which
   data about bond formed by 2 atoms of of force field types
   "Atom1_ff_type" and "Atom2_ff_type" are stored
   Returns -1 if such index has not been found, otherwise returns the index
   of last of bond_records which match the search parameters
   The bond record will be chosen on the basis of additional comparison
   the Bond_type_index properties if they are available
*****************************************************************/
  int indx = -1;
  int sz   = Bond_Records.size();
  int cmpr11,cmpr12,cmpr21,cmpr22,cmpr;

  for(int i=0;i<sz;i++){
    cmpr11 = cmpr12 = cmpr21 = cmpr22 = 0;
    if(Bond_Records[i].is_Atom1_ff_type && Bond_Records[i].is_Atom2_ff_type){
      cmpr11 = (Bond_Records[i].Atom1_ff_type == Atom1_ff_type);
      cmpr22 = (Bond_Records[i].Atom2_ff_type == Atom2_ff_type);
      cmpr12 = (Bond_Records[i].Atom1_ff_type == Atom2_ff_type);
      cmpr21 = (Bond_Records[i].Atom2_ff_type == Atom1_ff_type);
    }// if both Atom1_ff_int_type and Atom2_ff_int_type are defined
    // --------- Conclusions -------
    if((cmpr11&&cmpr22&&cmpr)||(cmpr12&&cmpr21&&cmpr)){ indx = i;  }
  }// for i
  return indx;
}


int ForceField::Bond_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type, int Bond_type_index, vector<int>& res){
/****************************************************************
   This function searches for indexes of Bond_Records vector in which
   data about bond formed by 2 atoms of of force field types
   "Atom1_ff_type" and "Atom2_ff_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of bond_records which match the search parameters
   The bond record will be chosen on the basis of additional comparison
   the Bond_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Bond_Records.size();
   int cmpr11,cmpr12,cmpr21,cmpr22,cmpr;

   for(int i=0;i<sz;i++){

       cmpr11 = cmpr12 = cmpr21 = cmpr22 = 0;
       if(Bond_Records[i].is_Atom1_ff_type && Bond_Records[i].is_Atom2_ff_type){

          cmpr11 = (Bond_Records[i].Atom1_ff_type == Atom1_ff_type);
          cmpr22 = (Bond_Records[i].Atom2_ff_type == Atom2_ff_type);
          cmpr12 = (Bond_Records[i].Atom1_ff_type == Atom2_ff_type);
          cmpr21 = (Bond_Records[i].Atom2_ff_type == Atom1_ff_type);

       }// if both Atom1_ff_int_type and Atom2_ff_int_type are defined


       if(Bond_Records[i].is_Bond_type_index){

          cmpr = (Bond_Records[i].Bond_type_index == Bond_type_index);

       }// if Bond_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if((cmpr11&&cmpr22&&cmpr)||(cmpr12&&cmpr21&&cmpr)){
          indx = i;
          res.push_back(i);
       }
   }// for i
   return indx;
}





int ForceField::Add_Bond_Record(Bond_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Bond_Record rec in array Bond_Records. If this bond type (bond record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the bond record rec is not valid - will not be added.
   If the Bond_type_index property is defined then it is taken into account
   during check on bond existence
*********************************************************************/
   int res = 1;
   int sz = Bond_Records.size();

   int int_types = (rec.is_Atom1_ff_int_type && rec.is_Atom2_ff_int_type );
   int sym_types = (rec.is_Atom1_ff_type && rec.is_Atom2_ff_type);
   int elt_types = (rec.is_Atom1_atomic_number && rec.is_Atom2_atomic_number);
   int is_bt_indx = rec.is_Bond_type_index;
   int bt_indx = 0;
   if(is_bt_indx){ bt_indx = rec.Bond_type_index; }

   vector<int> bnd_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_types
   last_indx = Bond_Record_Index(rec.Atom1_ff_int_type,rec.Atom2_ff_int_type,bt_indx,bnd_indxs);
   }
   else if(sym_types){
   // Search on the basis of symbolic types
   last_indx = Bond_Record_Index(rec.Atom1_ff_type,rec.Atom2_ff_type,bt_indx,bnd_indxs);
   }
   else if(elt_types){
   // Search on the basis of atomic numbers
   last_indx = Bond_Record_Index_by_Element(rec.Atom1_atomic_number,rec.Atom2_atomic_number,bt_indx,bnd_indxs);
   }
   else{ 
   // Record will not be added to the table of Bond_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2; 
   }

  
   if(last_indx==-2){ // In this case the record was not added to Bond_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The bond record is valid but does not exist in Bond_Records array
                           // so it will be added completely
       Bond_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar bond record is already defined - try to add new properties

       for(int i=0;i<bnd_indxs.size();i++){
           Bond_Records[bnd_indxs[i]].merge(rec);
       }
       res = 0;
   }


/* Old function body

   if(rec.is_Atom1_ff_type&&rec.is_Atom2_ff_type){

      for(int i=0;i<sz;i++){

         if(Bond_Records[i].is_Atom1_ff_type&&Bond_Records[i].is_Atom2_ff_type){
            // This type already exist
            if(((Bond_Records[i].Atom1_ff_type==rec.Atom1_ff_type)&&(Bond_Records[i].Atom2_ff_type==rec.Atom2_ff_type))||
               ((Bond_Records[i].Atom1_ff_type==rec.Atom2_ff_type)&&(Bond_Records[i].Atom2_ff_type==rec.Atom1_ff_type))
              )
              {
               if(Bond_Records[i].is_Bond_type_index){
               res = 0; break; 
               }// is_Bond_type_index
              }// Atom1_ff_type and Atom2_ff_type
         }// is_Atom1_ff_type and is_Atom2_ff_type
      }// for i

   }
   else{
    std::cout<<"To add bond type record to the force field its Atom1_ff_type and Atom2_ff_type should be defined"<<std::endl;
    res = -1;
   }

   if(res==1){
      // Before adding bond record to array - update Atom1(2)_ff_int_type
      sz = Atom_Records.size();
      for(int i=0;i<sz;i++){
          if(Atom_Records[i].is_Atom_ff_int_type){

             // Atom1
             if(Atom_Records[i].Atom_ff_type == rec.Atom1_ff_type){               
                   rec.Atom1_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom1_ff_int_type = 1;                
             }// if Atom1

             // Atom2  
             if(Atom_Records[i].Atom_ff_type == rec.Atom2_ff_type){               
                   rec.Atom2_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom2_ff_int_type = 1;                
             }// if Atom2

          }
      }// for i     
      Bond_Records.push_back(rec);
   }
*/

   return res;
}


int ForceField::Angle_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Atom3_ff_int_type){
/****************************************************************
   This function searches for index of Angle_Records vector in which
   data about angle formed by 3 atoms of of force field types
   "Atom1_ff_int_type", "Atom2_ff_int_type" and "Atom3_ff_int_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Angle_Records.size();

   for(int i=0;i<sz;i++){
       if(Angle_Records[i].is_Atom2_ff_int_type){
          if(Angle_Records[i].Atom2_ff_int_type==Atom2_ff_int_type){
             if(Angle_Records[i].is_Atom1_ff_int_type && Angle_Records[i].is_Atom3_ff_int_type){
                if(((Angle_Records[i].Atom1_ff_int_type==Atom1_ff_int_type)&&(Angle_Records[i].Atom3_ff_int_type==Atom3_ff_int_type) )||
                   ((Angle_Records[i].Atom1_ff_int_type==Atom3_ff_int_type)&&(Angle_Records[i].Atom3_ff_int_type==Atom1_ff_int_type) )
                  )
                  {
                    indx = i;
                    break;
                  } // if side atom types matches
             }// if side atom int_types are defined
          }// if center atom matches
       } // if Atom2_ff_int_type is defined
   }// for i

   return indx;
}

int ForceField::Angle_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Atom3_ff_int_type, int Angle_type_index, int order, vector<int>& res){
/****************************************************************
   This function searches for indexes of Atom_Records vector in which
   data about angle formed by 3 atoms of of force field types
   "Atom1_ff_int_type" "Atom2_ff_int_type" and "Atom3_ff_int_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of angle_records which match the search parameters
   The angle record will be chosen on the basis of additional comparison
   the Angle_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Angle_Records.size();
   int cmpr11,cmpr13,cmpr31,cmpr33,cmpr,cmpr2;

   for(int i=0;i<sz;i++){

       cmpr2 = 0;
       if(Angle_Records[i].is_Atom2_ff_int_type){
           cmpr2 = (Angle_Records[i].Atom2_ff_int_type == Atom2_ff_int_type);
       }

       if(cmpr2){

       cmpr11 = cmpr13 = cmpr31 = cmpr33 = 0;
       if(Angle_Records[i].is_Atom1_ff_int_type && Angle_Records[i].is_Atom3_ff_int_type){

          cmpr11 = (Angle_Records[i].Atom1_ff_int_type == Atom1_ff_int_type);
          cmpr33 = (Angle_Records[i].Atom3_ff_int_type == Atom3_ff_int_type);
          cmpr13 = (Angle_Records[i].Atom1_ff_int_type == Atom3_ff_int_type);
          cmpr31 = (Angle_Records[i].Atom3_ff_int_type == Atom1_ff_int_type);

       }// if both Atom1_ff_int_type and Atom3_ff_int_type are defined


       if(Angle_Records[i].is_Angle_type_index){

          cmpr = (Angle_Records[i].Angle_type_index == Angle_type_index);

       }// if Angle_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if(order==1){
           if(cmpr11&&cmpr33&&cmpr){
               indx = i;
               res.push_back(i);
           }
       }else{
           if((cmpr11&&cmpr33&&cmpr)||(cmpr13&&cmpr31&&cmpr)){
               indx = i;
               res.push_back(i);
           }
       }// else order==0

       }// if cmpr2==1

   }// for i
   return indx;
}


int ForceField::Angle_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type, std::string Atom3_ff_type, int Angle_type_index, int order, vector<int>& res){
/****************************************************************
   order = 0 - no direction 
   order = 1 - should be directed
   This function searches for indexes of Atom_Records vector in which
   data about angle formed by 3 atoms of of force field types
   "Atom1_ff_type" "Atom2_ff_type" and "Atom3_ff_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of angle_records which match the search parameters
   The angle record will be chosen on the basis of additional comparison
   the Angle_type_index properties if they are available
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Angle_Records.size();
   int cmpr11,cmpr13,cmpr31,cmpr33,cmpr,cmpr2;

   for(int i=0;i<sz;i++){

       cmpr2 = 0;
       if(Angle_Records[i].is_Atom2_ff_type){
           cmpr2 = (Angle_Records[i].Atom2_ff_type == Atom2_ff_type);
       }

       if(cmpr2){

       cmpr11 = cmpr13 = cmpr31 = cmpr33 = 0;
       if(Angle_Records[i].is_Atom1_ff_type && Angle_Records[i].is_Atom3_ff_type){

          cmpr11 = (Angle_Records[i].Atom1_ff_type == Atom1_ff_type);
          cmpr33 = (Angle_Records[i].Atom3_ff_type == Atom3_ff_type);
          cmpr13 = (Angle_Records[i].Atom1_ff_type == Atom3_ff_type);
          cmpr31 = (Angle_Records[i].Atom3_ff_type == Atom1_ff_type);

       }// if both Atom1_ff_int_type and Atom3_ff_int_type are defined


       if(Angle_Records[i].is_Angle_type_index){

          cmpr = (Angle_Records[i].Angle_type_index == Angle_type_index);

       }// if Angle_type_index is defined
       else{ cmpr = 1; }
       // --------- Conclusions -------
       if(order==1){
           if(cmpr11&&cmpr33&&cmpr){
               indx = i;
               res.push_back(i);
           }
       }else{
           if((cmpr11&&cmpr33&&cmpr)||(cmpr13&&cmpr31&&cmpr)){
               indx = i;
               res.push_back(i);
           }
       }// else

       }// if cmpr2==1

   }// for i
   return indx;
}


int ForceField::Angle_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type, std::string Atom3_ff_type){
/****************************************************************
   This function searches for indexes of Atom_Records vector in which
   data about angle formed by 3 atoms of of force field types
   "Atom1_ff_type" "Atom2_ff_type" and "Atom3_ff_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of angle_records which match the search parameters
   The angle record will be chosen on the basis of additional comparison
   the Angle_type_index properties if they are available
*****************************************************************/
  int indx = -1;
  int sz   = Angle_Records.size();
  int cmpr11,cmpr13,cmpr31,cmpr33,cmpr,cmpr2;
  for(int i=0;i<sz;i++){
    cmpr2 = 0;
    if(Angle_Records[i].is_Atom2_ff_type){
      cmpr2 = (Angle_Records[i].Atom2_ff_type == Atom2_ff_type);
    }
    if(cmpr2){
      cmpr11 = cmpr13 = cmpr31 = cmpr33 = 0;
      if(Angle_Records[i].is_Atom1_ff_type && Angle_Records[i].is_Atom3_ff_type){
        cmpr11 = (Angle_Records[i].Atom1_ff_type == Atom1_ff_type);
        cmpr33 = (Angle_Records[i].Atom3_ff_type == Atom3_ff_type);
        cmpr13 = (Angle_Records[i].Atom1_ff_type == Atom3_ff_type);
        cmpr31 = (Angle_Records[i].Atom3_ff_type == Atom1_ff_type);
      }// if both Atom1_ff_int_type and Atom3_ff_int_type are defined

      // --------- Conclusions -------
      if((cmpr11&&cmpr33&&cmpr)||(cmpr13&&cmpr31&&cmpr)){ indx = i;  }
    }// if cmpr2==1
  }// for i
  return indx;
}



int ForceField::Add_Angle_Record(Angle_Record rec,int order){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Angle_Record rec in array Angle_Records. If this angle type (angle record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the angle record rec is not valid - will not be added.
   If the Angle_type_index property is defined then it is taken into account
   during check on angle existence
   order parameter controls if we should consider order of indices
   (that is ijk is not the same as kji) or to treat them similarly
*********************************************************************/
   int res = 1;
   int sz = Angle_Records.size();

   int int_types = (rec.is_Atom1_ff_int_type && rec.is_Atom2_ff_int_type && rec.is_Atom3_ff_int_type);
   int sym_types = (rec.is_Atom1_ff_type && rec.is_Atom2_ff_type && rec.is_Atom3_ff_type);
//   int elt_types = (rec.is_Atom1_atomic_number && rec.is_Atom2_atomic_number && rec.is_Atom3_atomic_number);
   int is_at_indx = rec.is_Angle_type_index;
   int at_indx = 0;
   if(is_at_indx){ at_indx = rec.Angle_type_index; }

   vector<int> ang_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_types
   last_indx = Angle_Record_Index(rec.Atom1_ff_int_type,rec.Atom2_ff_int_type,rec.Atom3_ff_int_type,at_indx, order, ang_indxs);
   }
   else if(sym_types){
   // Search on the basis of symbolic types
   last_indx = Angle_Record_Index(rec.Atom1_ff_type,rec.Atom2_ff_type, rec.Atom3_ff_type, at_indx, order,ang_indxs);
   }
   else{
   // Record will not be added to the table of Angle_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2;
   }

   if(last_indx==-2){ // In this case the record was not added to Angle_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The angle record is valid but does not exist in Angle_Records array
                           // so it will be added completely
       Angle_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar bond record is already defined - try to add new properties

       for(int i=0;i<ang_indxs.size();i++){
           Angle_Records[ang_indxs[i]].merge(rec);
       }
       res = 0;
   }


   return res;
}



int ForceField::Add_Angle_Record(Angle_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Angle_Record rec in array Angle_Records. If this angle type (angle record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array (we did not add it this time). Returns -1
   if there is an error.
*********************************************************************/
   int res = 1;
   int sz = Angle_Records.size();

   if(rec.is_Atom1_ff_type&&rec.is_Atom2_ff_type&&rec.is_Atom3_ff_type){

      for(int i=0;i<sz;i++){

         if(Angle_Records[i].is_Atom1_ff_type&&Angle_Records[i].is_Atom2_ff_type&&Angle_Records[i].is_Atom3_ff_type){
            // This type already exist
	   if(Angle_Records[i].Atom2_ff_type==rec.Atom2_ff_type){
             if(((Angle_Records[i].Atom1_ff_type==rec.Atom1_ff_type)&&(Angle_Records[i].Atom2_ff_type==rec.Atom3_ff_type))||
                ((Angle_Records[i].Atom1_ff_type==rec.Atom3_ff_type)&&(Angle_Records[i].Atom2_ff_type==rec.Atom1_ff_type))
               ){ res = 0; break; }
	   }
         }
      }// for i

   }
   else{
    std::cout<<"To add angle type record to the force field its Atom1_ff_type, Atom2_ff_type  and Atom3_ff_type should be defined"<<std::endl;
    res = -1;
   }

   if(res==1){
      // Before adding angle record to array - update Atom1(2,3)_ff_int_type
      sz = Atom_Records.size();
      for(int i=0;i<sz;i++){
          if(Atom_Records[i].is_Atom_ff_int_type){

             // Atom1
             if(Atom_Records[i].Atom_ff_type == rec.Atom1_ff_type){               
                   rec.Atom1_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom1_ff_int_type = 1;                
             }// if Atom1

             // Atom2  
             if(Atom_Records[i].Atom_ff_type == rec.Atom2_ff_type){               
                   rec.Atom2_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom2_ff_int_type = 1;                
             }// if Atom2

             // Atom3  
             if(Atom_Records[i].Atom_ff_type == rec.Atom3_ff_type){               
                   rec.Atom3_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom3_ff_int_type = 1;                
             }// if Atom3

          }
      }// for i     
      Angle_Records.push_back(rec);
   }

   return res;
}


int ForceField::Dihedral_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Atom3_ff_int_type, int Atom4_ff_int_type){
/****************************************************************
   This function searches for index of Dihedral_Records vector in which
   data about dihedral formed by 4 atoms of of force field types
   "Atom1_ff_int_type", "Atom2_ff_int_type","Atom3_ff_int_type" and 
   "Atom4_ff_int_type" are stored
    Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Dihedral_Records.size();

   for(int i=0;i<sz;i++){
       if(Dihedral_Records[i].is_Atom2_ff_int_type&&Dihedral_Records[i].is_Atom3_ff_int_type){

          if((Dihedral_Records[i].Atom2_ff_int_type==Atom2_ff_int_type)&&(Dihedral_Records[i].Atom3_ff_int_type==Atom3_ff_int_type)){
             if(Dihedral_Records[i].is_Atom1_ff_int_type && Dihedral_Records[i].is_Atom4_ff_int_type){
                if((Dihedral_Records[i].Atom1_ff_int_type==Atom1_ff_int_type)&&(Dihedral_Records[i].Atom4_ff_int_type==Atom4_ff_int_type) )                  
                  {
                    indx = i;
                    break;
                  } // if side atom types matches: 1=1 and 4=4
             }// if side atoms int_types are defined
          }// if center 2 atoms matches : 2=2 and 3=3

          if((Dihedral_Records[i].Atom2_ff_int_type==Atom3_ff_int_type)&&(Dihedral_Records[i].Atom3_ff_int_type==Atom2_ff_int_type)){
             if(Dihedral_Records[i].is_Atom1_ff_int_type && Dihedral_Records[i].is_Atom4_ff_int_type){
                if((Dihedral_Records[i].Atom1_ff_int_type==Atom4_ff_int_type)&&(Dihedral_Records[i].Atom4_ff_int_type==Atom1_ff_int_type) )
                  {
                    indx = i;
                    break;
                  } // if side atom types matches: 1=4 and 4=1
             }// if side atoms int_types are defined
          }// if center 2 atoms matches : 2=3 and 3=2


       } // if Atom2_ff_int_type and Atom3_ff_int_type are defined
   }// for i

   return indx;
}

int ForceField::Improper_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Atom3_ff_int_type, int Atom4_ff_int_type){
/****************************************************************
   This function searches for index of Improper_Records vector in which
   data about improper formed by 4 atoms of of force field types
   "Atom1_ff_int_type", "Atom2_ff_int_type","Atom3_ff_int_type" and
   "Atom4_ff_int_type" are stored
    Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int res = 0;
   int sz   = Improper_Records.size();

   for(int i=0;i<sz;i++){
       if(Improper_Records[i].is_Atom1_ff_int_type &&
          Improper_Records[i].is_Atom2_ff_int_type &&
          Improper_Records[i].is_Atom3_ff_int_type &&
          Improper_Records[i].is_Atom4_ff_int_type)
       {

          res = 0;
          if(Improper_Records[i].Atom2_ff_int_type==Atom2_ff_int_type){

          // Now search match of any of 6 permutations of 3 periferal atoms
              if(Improper_Records[i].Atom1_ff_int_type==Atom1_ff_int_type){
                  if(Improper_Records[i].Atom3_ff_int_type==Atom3_ff_int_type && Improper_Records[i].Atom4_ff_int_type==Atom4_ff_int_type){
                      res = 1;
                  }
                  else if(Improper_Records[i].Atom3_ff_int_type==Atom4_ff_int_type && Improper_Records[i].Atom4_ff_int_type==Atom3_ff_int_type){
                      res = 1;
                  }
              }

              else if(Improper_Records[i].Atom3_ff_int_type==Atom3_ff_int_type){
                  if(Improper_Records[i].Atom1_ff_int_type==Atom1_ff_int_type && Improper_Records[i].Atom4_ff_int_type==Atom4_ff_int_type){
                      res = 1;
                  }
                  else if(Improper_Records[i].Atom1_ff_int_type==Atom4_ff_int_type && Improper_Records[i].Atom4_ff_int_type==Atom1_ff_int_type){
                      res = 1;
                  }
              }

              else if(Improper_Records[i].Atom4_ff_int_type==Atom4_ff_int_type){
                  if(Improper_Records[i].Atom3_ff_int_type==Atom3_ff_int_type && Improper_Records[i].Atom1_ff_int_type==Atom1_ff_int_type){
                      res = 1;
                  }
                  else if(Improper_Records[i].Atom3_ff_int_type==Atom1_ff_int_type && Improper_Records[i].Atom1_ff_int_type==Atom3_ff_int_type){
                      res = 1;
                  }
              }

          }// if 2 == 2 (central atoms are the same)

          if(res){ indx = i; break; }

       }// if all in types are defined

   }// for i

   return indx;
}



int ForceField::Dihedral_Record_Index(int Atom2_ff_int_type, int Atom3_ff_int_type){
/****************************************************************
   This function searches for index of Dihedral_Records vector based on
   only middle 2 atoms of force field types "Atom2_ff_int_type" and 
   "Atom3_ff_int_type". This is necessary because of in many cases
   the dihedral parmeter are not sensitive to 2 terminal atoms.
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Dihedral_Records.size();

   for(int i=0;i<sz;i++){
       if(Dihedral_Records[i].is_Atom2_ff_int_type&&Dihedral_Records[i].is_Atom3_ff_int_type){

          if(((Dihedral_Records[i].Atom2_ff_int_type==Atom2_ff_int_type)&&(Dihedral_Records[i].Atom3_ff_int_type==Atom3_ff_int_type))|| 
             ((Dihedral_Records[i].Atom2_ff_int_type==Atom3_ff_int_type)&&(Dihedral_Records[i].Atom3_ff_int_type==Atom2_ff_int_type))
            ){
                    indx = i;
                    break;              

             }// if center 2 atoms matches : 2=2, 3=3 or 2=3, 3=2


       } // if Atom2_ff_int_type and Atom3_ff_int_type are defined
   }// for i

   return indx;
}

int ForceField::Dihedral_Record_Index(int Atom1_ff_int_type, int Atom2_ff_int_type, int Atom3_ff_int_type, int Atom4_ff_int_type, int Dihedral_type_index, int order, vector<int>& res){
/****************************************************************
   This function searches for indexes of Dihedral_Records vector in which
   data about dihedral formed by 4 atoms of of force field types
   "Atom1_ff_int_type" "Atom2_ff_int_type" "Atom3_ff_int_type" and "
   Atom4_ff_int_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of dihedral_records which match the search parameters
   The dihedral record will be chosen on the basis of additional comparison
   the Dihedral_type_index properties if they are available
   order - is a parameter defining importance of indices order
   if order = 1 => ijkl is not the same as lkji
*****************************************************************/

   if(res.size()>0) { res.clear(); }
   int indx = -1;
   int sz   = Dihedral_Records.size();
   int cmpr11,cmpr14,cmpr41,cmpr23,cmpr32,cmpr22,cmpr33,cmpr44,cmpr;

   for(int i=0;i<sz;i++){

       if(Dihedral_Records[i].is_Dihedral_type_index){

          cmpr = (Dihedral_Records[i].Dihedral_type_index == Dihedral_type_index);

       }// if Dihedral_type_index is defined
       else{ cmpr = 1; }


       cmpr22 = cmpr33 = 0;
       if(Dihedral_Records[i].is_Atom1_ff_int_type &&
          Dihedral_Records[i].is_Atom2_ff_int_type &&
          Dihedral_Records[i].is_Atom3_ff_int_type &&
          Dihedral_Records[i].is_Atom4_ff_int_type
       ){
           cmpr11 = (Dihedral_Records[i].Atom1_ff_int_type == Atom1_ff_int_type);
           cmpr22 = (Dihedral_Records[i].Atom2_ff_int_type == Atom2_ff_int_type);
           cmpr33 = (Dihedral_Records[i].Atom3_ff_int_type == Atom3_ff_int_type);
           cmpr23 = (Dihedral_Records[i].Atom2_ff_int_type == Atom3_ff_int_type);
           cmpr32 = (Dihedral_Records[i].Atom3_ff_int_type == Atom2_ff_int_type);
           cmpr44 = (Dihedral_Records[i].Atom4_ff_int_type == Atom4_ff_int_type);
           cmpr14 = (Dihedral_Records[i].Atom1_ff_int_type == Atom4_ff_int_type);
           cmpr41 = (Dihedral_Records[i].Atom4_ff_int_type == Atom1_ff_int_type);

       }


       if(order==1){     

           if(cmpr&&cmpr11&&cmpr22&&cmpr33&&cmpr44){
               indx = i;
               res.push_back(i);
           }

       }else{

           if( (cmpr&&cmpr11&&cmpr22&&cmpr33&&cmpr44) ||
               (cmpr&&cmpr14&&cmpr23&&cmpr32&&cmpr41)
             ){
               indx = i;
               res.push_back(i);
           }

       }// order != 1


   }// for i
   return indx;
}

int ForceField::Dihedral_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type, std::string Atom3_ff_type, std::string Atom4_ff_type){
/****************************************************************
   This function searches for indexes of Dihedral_Records vector in which
   data about dihedral formed by 4 atoms of of force field types
   "Atom1_ff_int_type" "Atom2_ff_int_type" "Atom3_ff_int_type" and "
   Atom4_ff_int_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of dihedral_records which match the search parameters
   The dihedral record will be chosen on the basis of additional comparison
   the Dihedral_type_index properties if they are available
   order - is a parameter defining importance of indices order
   if order = 1 => ijkl is not the same as lkji
*****************************************************************/
  vector<int> tmp;
  return Dihedral_Record_Index(Atom1_ff_type,Atom2_ff_type,Atom3_ff_type,Atom4_ff_type,-1,0,tmp);
}

int ForceField::Dihedral_Record_Index(std::string Atom1_ff_type, std::string Atom2_ff_type, std::string Atom3_ff_type, std::string Atom4_ff_type, int Dihedral_type_index, int order, vector<int>& res){
/****************************************************************
   This function searches for indexes of Dihedral_Records vector in which
   data about dihedral formed by 4 atoms of of force field types
   "Atom1_ff_int_type" "Atom2_ff_int_type" "Atom3_ff_int_type" and "
   Atom4_ff_int_type" are stored
   Returns -1 if such index has not been found otherwise, returns the index
   of last of dihedral_records which match the search parameters
   The dihedral record will be chosen on the basis of additional comparison
   the Dihedral_type_index properties if they are available
   order - is a parameter defining importance of indices order
   if order = 1 => ijkl is not the same as lkji
   If Dihedral_type_index == -1 it will not be taken into account
*****************************************************************/
  if(res.size()>0) { res.clear(); }
  int indx = -1;
  int sz   = Dihedral_Records.size();
  int cmpr11,cmpr14,cmpr41,cmpr23,cmpr32,cmpr22,cmpr33,cmpr44,cmpr;
  for(int i=0;i<sz;i++){
    //--------------------------------------------------
    if(Dihedral_type_index==-1){ cmpr = 1; }
    else {
      if(Dihedral_Records[i].is_Dihedral_type_index){
        cmpr = (Dihedral_Records[i].Dihedral_type_index == Dihedral_type_index);
      }// if Dihedral_type_index is defined
      else{ cmpr = 1; }
    }

    //----------------------------------------------------
    cmpr22 = cmpr33 = 0;
       if(Dihedral_Records[i].is_Atom1_ff_type &&
          Dihedral_Records[i].is_Atom2_ff_type &&
          Dihedral_Records[i].is_Atom3_ff_type &&
          Dihedral_Records[i].is_Atom4_ff_type
       ){
           cmpr11 = (Dihedral_Records[i].Atom1_ff_type == Atom1_ff_type);
           cmpr22 = (Dihedral_Records[i].Atom2_ff_type == Atom2_ff_type);
           cmpr33 = (Dihedral_Records[i].Atom3_ff_type == Atom3_ff_type);
           cmpr23 = (Dihedral_Records[i].Atom2_ff_type == Atom3_ff_type);
           cmpr32 = (Dihedral_Records[i].Atom3_ff_type == Atom2_ff_type);
           cmpr44 = (Dihedral_Records[i].Atom4_ff_type == Atom4_ff_type);
           cmpr14 = (Dihedral_Records[i].Atom1_ff_type == Atom4_ff_type);
           cmpr41 = (Dihedral_Records[i].Atom4_ff_type == Atom1_ff_type);
       }


       if(order==1){

           if(cmpr&&cmpr11&&cmpr22&&cmpr33&&cmpr44){
               indx = i;
               res.push_back(i);
           }

       }else{

           if( (cmpr&&cmpr11&&cmpr22&&cmpr33&&cmpr44) ||
               (cmpr&&cmpr14&&cmpr23&&cmpr32&&cmpr41)
             ){
               indx = i;
               res.push_back(i);
           }

       }// order != 1


   }// for i
   return indx;
}

int ForceField::Add_Dihedral_Record(Dihedral_Record rec,int order){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Dihedral_Record rec in array Dihedral_Records. If this dihedral type 
   (dihedral record) exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the dihedral record rec is not valid - will not be added.
   If the Dihedral_type_index property is defined then it is taken into account
   during check on dihedral existence
   order parameter controls if we should consider order of indices
   (that is ijkl is not the same as lkji) or to treat them similarly
*********************************************************************/
   int res = 1;
   int sz = Dihedral_Records.size();

   int int_types = (rec.is_Atom1_ff_int_type && rec.is_Atom2_ff_int_type && rec.is_Atom3_ff_int_type && rec.is_Atom4_ff_int_type);
   int sym_types = (rec.is_Atom1_ff_type && rec.is_Atom2_ff_type && rec.is_Atom3_ff_type && rec.is_Atom4_ff_type);
//   int elt_types = (rec.is_Atom1_atomic_number && rec.is_Atom2_atomic_number && rec.is_Atom3_atomic_number);
   int is_dt_indx = rec.is_Dihedral_type_index;
   int dt_indx = 0;
   if(is_dt_indx){ dt_indx = rec.Dihedral_type_index; }

   vector<int> dih_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_types
   last_indx = Dihedral_Record_Index(rec.Atom1_ff_int_type,rec.Atom2_ff_int_type,rec.Atom3_ff_int_type, rec.Atom4_ff_int_type,dt_indx, order, dih_indxs);
   }
   else if(sym_types){
   // Search on the basis of symbolic types
   last_indx = Dihedral_Record_Index(rec.Atom1_ff_type,rec.Atom2_ff_type, rec.Atom3_ff_type, rec.Atom4_ff_type,dt_indx, order,dih_indxs);
   }
   else{
   // Record will not be added to the table of Dihedral_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2;
   }

   if(last_indx==-2){ // In this case the record was not added to Dihedral_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The dihedral record is valid but does not exist in Dihedral_Records array
                           // so it will be added completely
       Dihedral_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar bond record is already defined - try to add new properties

       for(int i=0;i<dih_indxs.size();i++){
           Dihedral_Records[dih_indxs[i]].merge(rec);
       }
       res = 0;
   }


   return res;
}





int ForceField::Add_Dihedral_Record(Dihedral_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Dihedral_Record rec in array Dihedral_Records. If this dihedral type (dihedral record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array (we did not add it this time). Returns -1
   if there is an error.
*********************************************************************/
   int res = 1;
   int sz = Dihedral_Records.size();

   if(rec.is_Atom1_ff_type&&rec.is_Atom2_ff_type&&rec.is_Atom3_ff_type&&rec.is_Atom4_ff_type){

      for(int i=0;i<sz;i++){

         if(Dihedral_Records[i].is_Atom1_ff_type&&Dihedral_Records[i].is_Atom2_ff_type
          &&Dihedral_Records[i].is_Atom3_ff_type&&Dihedral_Records[i].is_Atom4_ff_type){
           // This type already exist
	   if((Dihedral_Records[i].Atom2_ff_type==rec.Atom2_ff_type)&&(Dihedral_Records[i].Atom3_ff_type==rec.Atom3_ff_type)){
             if((Dihedral_Records[i].Atom1_ff_type==rec.Atom1_ff_type)&&(Dihedral_Records[i].Atom4_ff_type==rec.Atom4_ff_type))               
               { res = 0; break; }
	   }
           if((Dihedral_Records[i].Atom2_ff_type==rec.Atom3_ff_type)&&(Dihedral_Records[i].Atom3_ff_type==rec.Atom2_ff_type)){
             if((Dihedral_Records[i].Atom1_ff_type==rec.Atom4_ff_type)&&(Dihedral_Records[i].Atom4_ff_type==rec.Atom1_ff_type))               
               { res = 0; break; }
	   }

         }
      }// for i

   }
   else{
    std::cout<<"To add dihedral type record to the force field its Atom1_ff_type, Atom2_ff_type, Atom3_ff_type and Atom4_ff_type should be defined"<<std::endl;
    res = -1;
   }

   if(res==1){
      // Before adding angle record to array - update Atom1(2,3,4)_ff_int_type
      sz = Atom_Records.size();
      for(int i=0;i<sz;i++){
          if(Atom_Records[i].is_Atom_ff_int_type){

             // Atom1
             if(Atom_Records[i].Atom_ff_type == rec.Atom1_ff_type){               
                   rec.Atom1_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom1_ff_int_type = 1;                
             }// if Atom1

             // Atom2  
             if(Atom_Records[i].Atom_ff_type == rec.Atom2_ff_type){               
                   rec.Atom2_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom2_ff_int_type = 1;                
             }// if Atom2

             // Atom3  
             if(Atom_Records[i].Atom_ff_type == rec.Atom3_ff_type){               
                   rec.Atom3_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom3_ff_int_type = 1;                
             }// if Atom3

             // Atom4  
             if(Atom_Records[i].Atom_ff_type == rec.Atom4_ff_type){               
                   rec.Atom4_ff_int_type = Atom_Records[i].Atom_ff_int_type;
                   rec.is_Atom4_ff_int_type = 1;                
             }// if Atom4

          }
      }// for i     
      Dihedral_Records.push_back(rec);
   }

   return res;
}

int ForceField::Add_Improper_Record(Dihedral_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Dihedral_Record rec in array Improper_Records. If this improper type
   (improper record) exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the improper record rec is not valid - will not be added.

*********************************************************************/
   int res = 1;
   int sz = Improper_Records.size();

   int int_types = (rec.is_Atom1_ff_int_type && rec.is_Atom2_ff_int_type && rec.is_Atom3_ff_int_type && rec.is_Atom4_ff_int_type);
   int sym_types = (rec.is_Atom1_ff_type && rec.is_Atom2_ff_type && rec.is_Atom3_ff_type && rec.is_Atom4_ff_type);

   int is_dt_indx = rec.is_Dihedral_type_index;
   int dt_indx = 0;
   if(is_dt_indx){ dt_indx = rec.Dihedral_type_index; }

   vector<int> dih_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_types
   last_indx = Improper_Record_Index(rec.Atom1_ff_int_type,rec.Atom2_ff_int_type,rec.Atom3_ff_int_type, rec.Atom4_ff_int_type);
   dih_indxs.push_back(last_indx);
   }
/*
   else if(sym_types){
   // Search on the basis of symbolic types
   last_indx = Dihedral_Record_Index(rec.Atom1_ff_type,rec.Atom2_ff_type, rec.Atom3_ff_type, rec.Atom4_ff_type,dt_indx, order,dih_indxs);
   }
*/
   else{
   // Record will not be added to the table of Improper_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2;
   }


   if(last_indx==-2){ // In this case the record was not added to Improper_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The dihedral record is valid but does not exist in Improper_Records array
                           // so it will be added completely
       Improper_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar improper record is already defined - try to add new properties

       for(int i=0;i<dih_indxs.size();i++){
           Improper_Records[dih_indxs[i]].merge(rec);
       }
       res = 0;
   }


   return res;
}


int ForceField::Fragment_Record_Index(int Fragment_ff_int_type){
/****************************************************************
   This function searches for index of Fragment_Records vector in which
   data about fragment of force field type "Fragment_ff_int_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Fragment_Records.size();

   for(int i=0;i<sz;i++){
       if(Fragment_Records[i].is_Fragment_ff_int_type){
          if(Fragment_Records[i].Fragment_ff_int_type==Fragment_ff_int_type){
             indx = i;
             break;
          }
       }
   }

   return indx;
}

int ForceField::Fragment_Record_Index(std::string Fragment_ff_type){
/****************************************************************
   This function searches for index of Fragment_Records vector in which
   data about fragment of force field type "Fragment_ff_type" are stored
   Returns -1 if such index has not been found
*****************************************************************/

   int indx = -1;
   int sz   = Fragment_Records.size();

   for(int i=0;i<sz;i++){
       if(Fragment_Records[i].is_Fragment_ff_type){
          if(Fragment_Records[i].Fragment_ff_type==Fragment_ff_type){
             indx = i;
             break;
          }
       }
   }

   return indx;
}

int ForceField::Add_Fragment_Record(Fragment_Record rec){
/*********************************************************************
   This is user-interface function. It checks for existance of the
   Fragment_Record rec in array Fragment_Records. If this fragment type (fragment record)
   exists - do nothing if it does not - add it to array
   Returns 1 if this is a new record (it just has been added) and 0 if
   record existed in array but we added new properties to it. Returns -1
   if the atom record rec is not valid - will not be added.
*********************************************************************/
   int res = 1;
   int sz = Fragment_Records.size();

   int int_types = (rec.is_Fragment_ff_int_type);
   int sym_types = (rec.is_Fragment_ff_type);

   vector<int> at_indxs;
   int last_indx;

   if(int_types){
   // Search on the basis of int_type
   last_indx = Fragment_Record_Index(rec.Fragment_ff_int_type);
   at_indxs.push_back(last_indx);
   }
   else if(sym_types){
   // Search on the basis of symbolic type
   last_indx = Fragment_Record_Index(rec.Fragment_ff_type);
   at_indxs.push_back(last_indx);
   }
   else{
   // Record will not be added to the table of Fragment_Records since it is
   // useless (no labels defined neither integer nor symbolic)
   last_indx = -2;
   }

   if(last_indx==-2){ // In this case the record was not added to Fragment_Records array
       res = -1;
   }
   else if(last_indx==-1){ // The atom record is valid but does not exist in Fragment_Records array
                           // so it will be added completely
       Fragment_Records.push_back(rec);
       res = 1;
   }
   else{  // Similar atom record is already defined - try to add new properties

       for(int i=0;i<at_indxs.size();i++){
           Fragment_Records[at_indxs[i]].merge(rec);
       }
       res = 0;
   }

   return res;

}



int ForceField::show_atom_records(){

    for(int i=0;i<Atom_Records.size();i++){

        Atom_Records[i].show_info();
    }

    return 1;
}

int ForceField::show_bond_records(){

    for(int i=0;i<Bond_Records.size();i++){

        Bond_Records[i].show_info();
    }

    return 1;
}

int ForceField::show_angle_records(){

    for(int i=0;i<Angle_Records.size();i++){

        Angle_Records[i].show_info();
    }

    return 1;
}

int ForceField::show_dihedral_records(){

    for(int i=0;i<Dihedral_Records.size();i++){

        Dihedral_Records[i].show_info();
    }

    return 1;
}

int ForceField::show_improper_records(){

    for(int i=0;i<Improper_Records.size();i++){

        Improper_Records[i].show_info();
    }

    return 1;
}

int ForceField::show_fragment_records(){

    for(int i=0;i<Fragment_Records.size();i++){

        Fragment_Records[i].show_info();
    }

    return 1;
}


}// namespace libforcefield
}// namespace libhamiltonian_mm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian



