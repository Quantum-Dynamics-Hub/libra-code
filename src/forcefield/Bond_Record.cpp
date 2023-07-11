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

#include "Bond_Record.h"

/// liblibra namespace
namespace liblibra{

using namespace libio;

namespace libforcefield{


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

  if(is_Atom1_ff_type){  libio::save(pt,path+".Atom1_ff_type",Atom1_ff_type);    }
  if(is_Atom2_ff_type){  libio::save(pt,path+".Atom2_ff_type",Atom2_ff_type);    }
  if(is_Atom1_ff_int_type){  libio::save(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type);    }
  if(is_Atom2_ff_int_type){  libio::save(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type);    }
  if(is_Atom1_atomic_number){  libio::save(pt,path+".Atom1_atomic_number",Atom1_atomic_number);    }
  if(is_Atom2_atomic_number){  libio::save(pt,path+".Atom2_atomic_number",Atom2_atomic_number);    }
  if(is_Atom1_element){  libio::save(pt,path+".Atom1_element",Atom1_element);    }
  if(is_Atom2_element){  libio::save(pt,path+".Atom2_element",Atom2_element);    }

  if(is_Bond_type_index){  libio::save(pt,path+".Bond_type_index",Bond_type_index);    }
  if(is_Bond_r_eq){  libio::save(pt,path+".Bond_r_eq",Bond_r_eq);    }
  if(is_Bond_k_bond){  libio::save(pt,path+".Bond_k_bond",Bond_k_bond);    }
  if(is_Bond_D_bond){  libio::save(pt,path+".Bond_D_bond",Bond_D_bond);    }
  if(is_Bond_alpha){  libio::save(pt,path+".Bond_alpha",Bond_alpha);    }
  if(is_Bond_r_eq_ref){  libio::save(pt,path+".Bond_r_eq_ref",Bond_r_eq_ref);    }
  if(is_Bond_k_bond_ref){  libio::save(pt,path+".Bond_k_bond_ref",Bond_k_bond_ref);    }
  if(is_Bond_shift_elec){  libio::save(pt,path+".Bond_shift_elec",Bond_shift_elec);    }
  if(is_Bond_bci){  libio::save(pt,path+".Bond_bci",Bond_bci);    }

  if(is_Bond_wij){  libio::save(pt,path+".Bond_wij",Bond_wij);    }
  if(is_Bond_wij_1){  libio::save(pt,path+".Bond_wij_1",Bond_wij_1);    }
  if(is_Bond_wij_2){  libio::save(pt,path+".Bond_wij_2",Bond_wij_2);    }
  if(is_Bond_alpij){  libio::save(pt,path+".Bond_alpij",Bond_alpij);    }
  if(is_Bond_alpij_1){  libio::save(pt,path+".Bond_alpij_1",Bond_alpij_1);    }
  if(is_Bond_alpij_2){  libio::save(pt,path+".Bond_alpij_2",Bond_alpij_2);    }
  if(is_Bond_rij_1){  libio::save(pt,path+".Bond_rij_1",Bond_rij_1);    }
  if(is_Bond_rij_2){  libio::save(pt,path+".Bond_rij_2",Bond_rij_2);    }

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

  libio::load(pt,path+".Atom1_ff_type",Atom1_ff_type,is_Atom1_ff_type); if(is_Atom1_ff_type==1) { status=1;}
  libio::load(pt,path+".Atom2_ff_type",Atom2_ff_type,is_Atom2_ff_type); if(is_Atom2_ff_type==1) { status=1;}
  libio::load(pt,path+".Atom1_ff_int_type",Atom1_ff_int_type,is_Atom1_ff_int_type); if(is_Atom1_ff_int_type==1) { status=1;}
  libio::load(pt,path+".Atom2_ff_int_type",Atom2_ff_int_type,is_Atom2_ff_int_type); if(is_Atom2_ff_int_type==1) { status=1;}
  libio::load(pt,path+".Atom1_atomic_number",Atom1_atomic_number,is_Atom1_atomic_number); if(is_Atom1_atomic_number==1) { status=1;}
  libio::load(pt,path+".Atom2_atomic_number",Atom2_atomic_number,is_Atom2_atomic_number); if(is_Atom2_atomic_number==1) { status=1;}
  libio::load(pt,path+".Atom1_element",Atom1_element,is_Atom1_element); if(is_Atom1_element==1) { status=1;}
  libio::load(pt,path+".Atom2_element",Atom2_element,is_Atom2_element); if(is_Atom2_element==1) { status=1;}

  libio::load(pt,path+".Bond_type_index",Bond_type_index,is_Bond_type_index); if(is_Bond_type_index==1) { status=1;}
  libio::load(pt,path+".Bond_r_eq",Bond_r_eq,is_Bond_r_eq); if(is_Bond_r_eq==1) { status=1;}
  libio::load(pt,path+".Bond_k_bond",Bond_k_bond,is_Bond_k_bond); if(is_Bond_k_bond==1) { status=1;}
  libio::load(pt,path+".Bond_D_bond",Bond_D_bond,is_Bond_D_bond); if(is_Bond_D_bond==1) { status=1;}
  libio::load(pt,path+".Bond_alpha",Bond_alpha,is_Bond_alpha); if(is_Bond_alpha==1) { status=1;}
  libio::load(pt,path+".Bond_r_eq_ref",Bond_r_eq_ref,is_Bond_r_eq_ref); if(is_Bond_r_eq_ref==1) { status=1;}
  libio::load(pt,path+".Bond_k_bond_ref",Bond_k_bond_ref,is_Bond_k_bond_ref); if(is_Bond_k_bond_ref==1) { status=1;}
  libio::load(pt,path+".Bond_shift_elec",Bond_shift_elec,is_Bond_shift_elec); if(is_Bond_shift_elec==1) { status=1;}
  libio::load(pt,path+".Bond_bci",Bond_bci,is_Bond_bci); if(is_Bond_bci==1) { status=1;}

  libio::load(pt,path+".Bond_wij",Bond_wij,is_Bond_wij); if(is_Bond_wij==1) { status=1;}
  libio::load(pt,path+".Bond_wij_1",Bond_wij_1,is_Bond_wij_1); if(is_Bond_wij_1==1) { status=1;}
  libio::load(pt,path+".Bond_wij_2",Bond_wij_2,is_Bond_wij_2); if(is_Bond_wij_2==1) { status=1;}
  libio::load(pt,path+".Bond_alpij",Bond_alpij,is_Bond_alpij); if(is_Bond_alpij==1) { status=1;}
  libio::load(pt,path+".Bond_alpij_1",Bond_alpij_1,is_Bond_alpij_1); if(is_Bond_alpij_1==1) { status=1;}
  libio::load(pt,path+".Bond_alpij_2",Bond_alpij_2,is_Bond_alpij_2); if(is_Bond_alpij_2==1) { status=1;}
  libio::load(pt,path+".Bond_rij_1",Bond_rij_1,is_Bond_rij_1); if(is_Bond_rij_1==1) { status=1;}
  libio::load(pt,path+".Bond_rij_2",Bond_rij_2,is_Bond_rij_2); if(is_Bond_rij_2==1) { status=1;}

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



}// namespace libforcefield
}// namespace liblibra


