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

#include "Angle_Record.h"

/// liblibra namespace
namespace liblibra{

using namespace libio;


namespace libforcefield{

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
   set_value(is_Angle_coordination, Angle_coordination, at, "Angle_coordination");

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
   is_Angle_coordination = 0;


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
   if(at.is_Angle_coordination){ Angle_coordination = at.Angle_coordination; is_Angle_coordination = 1;}

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
   if(!is_Angle_coordination && at.is_Angle_coordination) { Angle_coordination = at.Angle_coordination; is_Angle_coordination = 1;}

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
   if(is_Angle_coordination)   {std::cout<<"Angle_coordination = "<<Angle_coordination<<" some units"<<std::endl; }


   std::cout<<std::endl;

   return 1;

}

void Angle_Record::save(boost::property_tree::ptree& pt,std::string path){

  if(is_Atom1_ff_type){  libio::save(pt,path+".Atom1_ff_type",Atom1_ff_type);    }
  if(is_Atom2_ff_type){  libio::save(pt,path+".Atom2_ff_type",Atom2_ff_type);    }
  if(is_Atom3_ff_type){  libio::save(pt,path+".Atom3_ff_type",Atom3_ff_type);    }
  if(is_Atom1_ff_int_type){  libio::save(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type);    }
  if(is_Atom2_ff_int_type){  libio::save(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type);    }
  if(is_Atom3_ff_int_type){  libio::save(pt,path+".Atom3_ff_int_type",Atom3_ff_int_type);    }

  if(is_Angle_type_index){  libio::save(pt,path+".Angle_type_index",Angle_type_index);    }
  if(is_Angle_theta_eq){  libio::save(pt,path+".Angle_theta_eq",Angle_theta_eq);    }
  if(is_Angle_k_angle){  libio::save(pt,path+".Angle_k_angle",Angle_k_angle);    }
  if(is_Angle_r_eq){  libio::save(pt,path+".Angle_r_eq",Angle_r_eq);    }
  if(is_Angle_k_ub){  libio::save(pt,path+".Angle_k_ub",Angle_k_ub);    }
  if(is_Angle_kijk_sb){  libio::save(pt,path+".Angle_kijk_sb",Angle_kijk_sb);    }
  if(is_Angle_kkji_sb){  libio::save(pt,path+".Angle_kkji_sb",Angle_kkji_sb);    }
  if(is_Angle_coordination){  libio::save(pt,path+".Angle_coordination",Angle_coordination);    }

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

  libio::load(pt,path+".Atom1_ff_type",Atom1_ff_type,is_Atom1_ff_type); if(is_Atom1_ff_type==1) { status=1;}
  libio::load(pt,path+".Atom2_ff_type",Atom2_ff_type,is_Atom2_ff_type); if(is_Atom2_ff_type==1) { status=1;}
  libio::load(pt,path+".Atom3_ff_type",Atom3_ff_type,is_Atom3_ff_type); if(is_Atom3_ff_type==1) { status=1;}
  libio::load(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type,is_Atom1_ff_int_type); if(is_Atom1_ff_int_type==1) { status=1;}
  libio::load(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type,is_Atom2_ff_int_type); if(is_Atom2_ff_int_type==1) { status=1;}
  libio::load(pt,path+".Atom3_ff_int_type",Atom3_ff_int_type,is_Atom3_ff_int_type); if(is_Atom3_ff_int_type==1) { status=1;}

  libio::load(pt,path+".Angle_type_index",Angle_type_index,is_Angle_type_index); if(is_Angle_type_index==1) { status=1;}
  libio::load(pt,path+".Angle_theta_eq",Angle_theta_eq,is_Angle_theta_eq); if(is_Angle_theta_eq==1) { status=1;}
  libio::load(pt,path+".Angle_k_angle",Angle_k_angle,is_Angle_k_angle); if(is_Angle_k_angle==1) { status=1;}
  libio::load(pt,path+".Angle_r_eq",Angle_r_eq,is_Angle_r_eq); if(is_Angle_r_eq==1) { status=1;}
  libio::load(pt,path+".Angle_k_ub",Angle_k_ub,is_Angle_k_ub); if(is_Angle_k_ub==1) { status=1;}
  libio::load(pt,path+".Angle_kijk_sb",Angle_kijk_sb,is_Angle_kijk_sb); if(is_Angle_kijk_sb==1) { status=1;}
  libio::load(pt,path+".Angle_kkji_sb",Angle_kkji_sb,is_Angle_kkji_sb); if(is_Angle_kkji_sb==1) { status=1;}
  libio::load(pt,path+".Angle_coordination",Angle_coordination,is_Angle_coordination); if(is_Angle_coordination==1) { status=1;}

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




}// namespace libforcefield
}// liblibra

