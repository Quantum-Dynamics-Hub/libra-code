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

#include "Dihedral_Record.h"

/// liblibra namespace
namespace liblibra{

using namespace libio;

namespace libforcefield{

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

  if(is_Atom1_ff_type){  libio::save(pt,path+".Atom1_ff_type",Atom1_ff_type);    }
  if(is_Atom2_ff_type){  libio::save(pt,path+".Atom2_ff_type",Atom2_ff_type);    }
  if(is_Atom3_ff_type){  libio::save(pt,path+".Atom3_ff_type",Atom3_ff_type);    }
  if(is_Atom4_ff_type){  libio::save(pt,path+".Atom4_ff_type",Atom4_ff_type);    }
  if(is_Atom1_ff_int_type){  libio::save(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type);    }
  if(is_Atom2_ff_int_type){  libio::save(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type);    }
  if(is_Atom3_ff_int_type){  libio::save(pt,path+".Atom3_ff_int_type",Atom3_ff_int_type);    }
  if(is_Atom4_ff_int_type){  libio::save(pt,path+".Atom4_ff_int_type",Atom4_ff_int_type);    }

  if(is_Dihedral_type_index){  libio::save(pt,path+".Dihedral_type_index",Dihedral_type_index);    }
  if(is_Dihedral_vphi){  libio::save(pt,path+".Dihedral_vphi",Dihedral_vphi);    }
  if(is_Dihedral_vphi1){  libio::save(pt,path+".Dihedral_vphi1",Dihedral_vphi1);    }
  if(is_Dihedral_vphi2){  libio::save(pt,path+".Dihedral_vphi2",Dihedral_vphi2);    }
  if(is_Dihedral_vphi3){  libio::save(pt,path+".Dihedral_vphi3",Dihedral_vphi3);    }
  if(is_Dihedral_phase){  libio::save(pt,path+".Dihedral_phase",Dihedral_phase);    }
  if(is_Dihedral_mult){  libio::save(pt,path+".Dihedral_mult",Dihedral_mult);    }

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

  libio::load(pt,path+".Atom1_ff_type",Atom1_ff_type,is_Atom1_ff_type); if(is_Atom1_ff_type==1) { status=1;}
  libio::load(pt,path+".Atom2_ff_type",Atom2_ff_type,is_Atom2_ff_type); if(is_Atom2_ff_type==1) { status=1;}
  libio::load(pt,path+".Atom3_ff_type",Atom3_ff_type,is_Atom3_ff_type); if(is_Atom3_ff_type==1) { status=1;}
  libio::load(pt,path+".Atom4_ff_type",Atom4_ff_type,is_Atom4_ff_type); if(is_Atom4_ff_type==1) { status=1;}
  libio::load(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type,is_Atom1_ff_int_type); if(is_Atom1_ff_int_type==1) { status=1;}
  libio::load(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type,is_Atom2_ff_int_type); if(is_Atom2_ff_int_type==1) { status=1;}
  libio::load(pt,path+".Atom3_ff_int_type",Atom3_ff_int_type,is_Atom3_ff_int_type); if(is_Atom3_ff_int_type==1) { status=1;}
  libio::load(pt,path+".Atom4_ff_int_type",Atom4_ff_int_type,is_Atom4_ff_int_type); if(is_Atom4_ff_int_type==1) { status=1;}

  libio::load(pt,path+".Dihedral_type_index",Dihedral_type_index,is_Dihedral_type_index); if(is_Dihedral_type_index==1) { status=1;}
  libio::load(pt,path+".Dihedral_vphi",Dihedral_vphi,is_Dihedral_vphi); if(is_Dihedral_vphi==1) { status=1;}
  libio::load(pt,path+".Dihedral_vphi1",Dihedral_vphi1,is_Dihedral_vphi1); if(is_Dihedral_vphi1==1) { status=1;}
  libio::load(pt,path+".Dihedral_vphi2",Dihedral_vphi2,is_Dihedral_vphi2); if(is_Dihedral_vphi2==1) { status=1;}
  libio::load(pt,path+".Dihedral_vphi3",Dihedral_vphi3,is_Dihedral_vphi3); if(is_Dihedral_vphi3==1) { status=1;}
  libio::load(pt,path+".Dihedral_phase",Dihedral_phase,is_Dihedral_phase); if(is_Dihedral_phase==1) { status=1;}
  libio::load(pt,path+".Dihedral_mult",Dihedral_mult,is_Dihedral_mult); if(is_Dihedral_mult==1) { status=1;}

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




}// namespace libforcefield
}// namespace liblibra


