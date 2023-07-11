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

#include "Atom_Record.h"

/// liblibra namespace
namespace liblibra{

using namespace libio;

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



int Atom_Record::show_info() const{

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

  if(is_Atom_ff_type){  libio::save(pt,path+".Atom_ff_type",Atom_ff_type);    }
  if(is_Atom_ff_type_H){  libio::save(pt,path+".Atom_ff_type_H",Atom_ff_type_H);    }
  if(is_Atom_ff_int_type){  libio::save(pt,path+".Atom_ff_int_type",Atom_ff_int_type);    }
  if(is_Atom_element){  libio::save(pt,path+".Atom_element",Atom_element);    }
  if(is_Atom_atomic_number){  libio::save(pt,path+".Atom_atomic_number",Atom_atomic_number);    }
  if(is_Atom_electronegativity){  libio::save(pt,path+".Atom_electronegativity",Atom_electronegativity);    }
  if(is_Atom_partial_charge){  libio::save(pt,path+".Atom_partial_charge",Atom_partial_charge);    }
  if(is_Atom_ff_eq_int_type2){  libio::save(pt,path+".Atom_ff_eq_int_type2",Atom_ff_eq_int_type2);    }
  if(is_Atom_ff_eq_int_type3){  libio::save(pt,path+".Atom_ff_eq_int_type3",Atom_ff_eq_int_type3);    }
  if(is_Atom_ff_eq_int_type4){  libio::save(pt,path+".Atom_ff_eq_int_type4",Atom_ff_eq_int_type4);    }
  if(is_Atom_ff_eq_int_type5){  libio::save(pt,path+".Atom_ff_eq_int_type5",Atom_ff_eq_int_type5);    }
  if(is_Atom_radius){  libio::save(pt,path+".Atom_radius",Atom_radius);    }
  if(is_Atom_Z_star){  libio::save(pt,path+".Atom_Z_star",Atom_Z_star);    }
  if(is_Atom_theta){  libio::save(pt,path+".Atom_theta",Atom_theta);    }
  if(is_Atom_sigma){  libio::save(pt,path+".Atom_sigma",Atom_sigma);    }
  if(is_Atom_epsilon){  libio::save(pt,path+".Atom_epsilon",Atom_epsilon);    }
  if(is_Atom_GMP){  libio::save(pt,path+".Atom_GMP",Atom_GMP);    }
  if(is_Atom_crd){  libio::save(pt,path+".Atom_crd",Atom_crd);    }
  if(is_Atom_val){  libio::save(pt,path+".Atom_val",Atom_val);    }
  if(is_Atom_pilp){  libio::save(pt,path+".Atom_pilp",Atom_pilp);    }
  if(is_Atom_mltb){  libio::save(pt,path+".Atom_mltb",Atom_mltb);    }
  if(is_Atom_arom){  libio::save(pt,path+".Atom_arom",Atom_arom);    }
  if(is_Atom_lin){  libio::save(pt,path+".Atom_lin",Atom_lin);    }
  if(is_Atom_sbmb){  libio::save(pt,path+".Atom_sbmb",Atom_sbmb);    }
  if(is_Atom_alpha){  libio::save(pt,path+".Atom_alpha",Atom_alpha);    }
  if(is_Atom_N_eff){  libio::save(pt,path+".Atom_N_eff",Atom_N_eff);    }
  if(is_Atom_A_scale){  libio::save(pt,path+".Atom_A_scale",Atom_A_scale);    }
  if(is_Atom_G_scale){  libio::save(pt,path+".Atom_G_scale",Atom_G_scale);    }
  if(is_Atom_DAN){  libio::save(pt,path+".Atom_DAN",Atom_DAN);    }
  if(is_Atom_pbci){  libio::save(pt,path+".Atom_pbci",Atom_pbci);    }
  if(is_Atom_fcadj){  libio::save(pt,path+".Atom_fcadj",Atom_fcadj);    }
  if(is_Atom_dative){  libio::save(pt,path+".Atom_dative",Atom_dative);    }
  if(is_Atom_brdr1){  libio::save(pt,path+".Atom_brdr1",Atom_brdr1);    }
  if(is_Atom_brdr2){  libio::save(pt,path+".Atom_brdr2",Atom_brdr2);    }
  if(is_Atom_brdr3){  libio::save(pt,path+".Atom_brdr3",Atom_brdr3);    }
  if(is_Atom_such_n){  libio::save(pt,path+".Atom_such_n",Atom_such_n);    }
  if(is_Atom_such_m){  libio::save(pt,path+".Atom_such_m",Atom_such_m);    }
  if(is_Atom_such_a){  libio::save(pt,path+".Atom_such_a",Atom_such_a);    }
  if(is_Atom_such_D){  libio::save(pt,path+".Atom_such_D",Atom_such_D);    }
  if(is_Atom_such_c){  libio::save(pt,path+".Atom_such_c",Atom_such_c);    }


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

  libio::load(pt,path+".Atom_ff_type",Atom_ff_type,is_Atom_ff_type); if(is_Atom_ff_type==1) { status=1;}
  libio::load(pt,path+".Atom_ff_type_H",Atom_ff_type_H,is_Atom_ff_type_H); if(is_Atom_ff_type_H==1) { status=1;}
  libio::load(pt,path+".Atom_ff_int_type",Atom_ff_int_type,is_Atom_ff_int_type); if(is_Atom_ff_int_type==1) { status=1;}
  libio::load(pt,path+".Atom_element",Atom_element,is_Atom_element); if(is_Atom_element==1) { status=1;}
  libio::load(pt,path+".Atom_atomic_number",Atom_atomic_number,is_Atom_atomic_number); if(is_Atom_atomic_number==1) { status=1;}

  libio::load(pt,path+".Atom_electronegativity",Atom_electronegativity,is_Atom_electronegativity); if(is_Atom_electronegativity==1) { status=1;}
  libio::load(pt,path+".Atom_partial_charge",Atom_partial_charge,is_Atom_partial_charge); if(is_Atom_partial_charge==1) { status=1;}
  libio::load(pt,path+".Atom_ff_eq_int_type2",Atom_ff_eq_int_type2,is_Atom_ff_eq_int_type2); if(is_Atom_ff_eq_int_type2==1) { status=1;}
  libio::load(pt,path+".Atom_ff_eq_int_type3",Atom_ff_eq_int_type3,is_Atom_ff_eq_int_type3); if(is_Atom_ff_eq_int_type3==1) { status=1;}
  libio::load(pt,path+".Atom_ff_eq_int_type4",Atom_ff_eq_int_type4,is_Atom_ff_eq_int_type4); if(is_Atom_ff_eq_int_type4==1) { status=1;}
  libio::load(pt,path+".Atom_ff_eq_int_type5",Atom_ff_eq_int_type5,is_Atom_ff_eq_int_type5); if(is_Atom_ff_eq_int_type5==1) { status=1;}

  libio::load(pt,path+".Atom_radius",Atom_radius,is_Atom_radius); if(is_Atom_radius==1) { status=1;}
  libio::load(pt,path+".Atom_Z_star",Atom_Z_star,is_Atom_Z_star); if(is_Atom_Z_star==1) { status=1;}
  libio::load(pt,path+".Atom_theta",Atom_theta,is_Atom_theta); if(is_Atom_theta==1) { status=1;}
  libio::load(pt,path+".Atom_sigma",Atom_sigma,is_Atom_sigma); if(is_Atom_sigma==1) { status=1;}
  libio::load(pt,path+".Atom_epsilon",Atom_epsilon,is_Atom_epsilon); if(is_Atom_epsilon==1) { status=1;}
  libio::load(pt,path+".Atom_GMP",Atom_GMP,is_Atom_GMP); if(is_Atom_GMP==1) { status=1;}

  libio::load(pt,path+".Atom_crd",Atom_crd,is_Atom_crd); if(is_Atom_crd==1) { status=1;}
  libio::load(pt,path+".Atom_val",Atom_val,is_Atom_val); if(is_Atom_val==1) { status=1;}
  libio::load(pt,path+".Atom_pilp",Atom_pilp,is_Atom_pilp); if(is_Atom_pilp==1) { status=1;}
  libio::load(pt,path+".Atom_mltb",Atom_mltb,is_Atom_mltb); if(is_Atom_mltb==1) { status=1;}
  libio::load(pt,path+".Atom_arom",Atom_arom,is_Atom_arom); if(is_Atom_arom==1) { status=1;}
  libio::load(pt,path+".Atom_lin",Atom_lin,is_Atom_lin); if(is_Atom_lin==1) { status=1;}
  libio::load(pt,path+".Atom_sbmb",Atom_sbmb,is_Atom_sbmb); if(is_Atom_sbmb==1) { status=1;}
  libio::load(pt,path+".Atom_alpha",Atom_alpha,is_Atom_alpha); if(is_Atom_alpha==1) { status=1;}
  libio::load(pt,path+".Atom_N_eff",Atom_N_eff,is_Atom_N_eff); if(is_Atom_N_eff==1) { status=1;}
  libio::load(pt,path+".Atom_A_scale",Atom_A_scale,is_Atom_A_scale); if(is_Atom_A_scale==1) { status=1;}
  libio::load(pt,path+".Atom_G_scale",Atom_G_scale,is_Atom_G_scale); if(is_Atom_G_scale==1) { status=1;}
  libio::load(pt,path+".Atom_DAN",Atom_DAN,is_Atom_DAN); if(is_Atom_DAN==1) { status=1;}
  libio::load(pt,path+".Atom_pbci",Atom_pbci,is_Atom_pbci); if(is_Atom_pbci==1) { status=1;}
  libio::load(pt,path+".Atom_fcadj",Atom_fcadj,is_Atom_fcadj); if(is_Atom_fcadj==1) { status=1;}

  libio::load(pt,path+".Atom_dative",Atom_dative,is_Atom_dative); if(is_Atom_dative==1) { status=1;}
  libio::load(pt,path+".Atom_brdr1",Atom_brdr1,is_Atom_brdr1); if(is_Atom_brdr1==1) { status=1;}
  libio::load(pt,path+".Atom_brdr2",Atom_brdr2,is_Atom_brdr2); if(is_Atom_brdr2==1) { status=1;}
  libio::load(pt,path+".Atom_brdr3",Atom_brdr3,is_Atom_brdr3); if(is_Atom_brdr3==1) { status=1;}

  libio::load(pt,path+".Atom_such_n",Atom_such_n,is_Atom_such_n); if(is_Atom_such_n==1) { status=1;}
  libio::load(pt,path+".Atom_such_m",Atom_such_m,is_Atom_such_m); if(is_Atom_such_m==1) { status=1;}
  libio::load(pt,path+".Atom_such_a",Atom_such_a,is_Atom_such_a); if(is_Atom_such_a==1) { status=1;}
  libio::load(pt,path+".Atom_such_D",Atom_such_D,is_Atom_such_D); if(is_Atom_such_D==1) { status=1;}
  libio::load(pt,path+".Atom_such_c",Atom_such_c,is_Atom_such_c); if(is_Atom_such_c==1) { status=1;}


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





}// namespace libforcefield
}// namespace liblibra


