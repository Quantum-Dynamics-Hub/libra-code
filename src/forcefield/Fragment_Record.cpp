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

#include "Fragment_Record.h"

/// liblibra namespace
namespace liblibra{

//using namespace liblinalg;
using namespace libio;

namespace libforcefield{

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

  if(is_Fragment_ff_type){  libio::save(pt,path+".Fragment_ff_type",Fragment_ff_type);    }
  if(is_Fragment_ff_int_type){  libio::save(pt,path+".Fragment_ff_int_type",Fragment_ff_int_type);    }

  if(is_Fragment_di){  libio::save(pt,path+".Fragment_di",Fragment_di);    }
  if(is_Fragment_li){  libio::save(pt,path+".Fragment_li",Fragment_li);    }
  if(is_Fragment_e0){  libio::save(pt,path+".Fragment_e0",Fragment_e0);    }
  if(is_Fragment_rat){  libio::save(pt,path+".Fragment_rat",Fragment_rat);    }
  if(is_Fragment_dw){  libio::save(pt,path+".Fragment_dw",Fragment_dw);    }
  if(is_Fragment_mu){  libio::save(pt,path+".Fragment_mu",Fragment_mu);    }
  if(is_Fragment_nu){  libio::save(pt,path+".Fragment_nu",Fragment_nu);    }

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

  libio::load(pt,path+".Fragment_ff_type",Fragment_ff_type,is_Fragment_ff_type); if(is_Fragment_ff_type==1) { status=1;}
  libio::load(pt,path+".Fragment_ff_int_type",Fragment_ff_int_type,is_Fragment_ff_int_type); if(is_Fragment_ff_int_type==1) { status=1;}

  libio::load(pt,path+".Fragment_di",Fragment_di,is_Fragment_di); if(is_Fragment_di==1) { status=1;}
  libio::load(pt,path+".Fragment_li",Fragment_li,is_Fragment_li); if(is_Fragment_li==1) { status=1;}
  libio::load(pt,path+".Fragment_e0",Fragment_e0,is_Fragment_e0); if(is_Fragment_e0==1) { status=1;}
  libio::load(pt,path+".Fragment_rat",Fragment_rat,is_Fragment_rat); if(is_Fragment_rat==1) { status=1;}
  libio::load(pt,path+".Fragment_dw",Fragment_dw,is_Fragment_dw); if(is_Fragment_dw==1) { status=1;}
  libio::load(pt,path+".Fragment_mu",Fragment_mu,is_Fragment_mu); if(is_Fragment_mu==1) { status=1;}
  libio::load(pt,path+".Fragment_nu",Fragment_nu,is_Fragment_nu); if(is_Fragment_nu==1) { status=1;}

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




}// namespace libforcefield
}// namespace liblibra


