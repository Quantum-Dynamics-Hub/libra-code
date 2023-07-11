/*********************************************************************************
* Copyright (C) 2017-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Interactions_2_Body.cpp
  \brief The file implements functions for molecular-mechanical Hamiltonian calculations as well
  as the classes for organizing such computations in an object-oriented way: basic functionality of the class
*/

#include "Interactions_2_Body.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_mm namespace
namespace libhamiltonian_mm{


using namespace libpot;


void Interaction_2_Body::set_coords(VECTOR* r1_, VECTOR* r2_){  r[0] = r1_; r[1] = r2_; }
void Interaction_2_Body::set_coords(VECTOR& r1_, VECTOR& r2_){  r[0] = &r1_; r[1] = &r2_; }

void Interaction_2_Body::set_transl(VECTOR* t1_, VECTOR* t2_){  t[0] = t1_; t[1] = t2_; }
void Interaction_2_Body::set_transl(VECTOR& t1_, VECTOR& t2_){  t[0] = &t1_; t[1] = &t2_; }

void Interaction_2_Body::set_forces(VECTOR* f1_, VECTOR* f2_){  f[0] = f1_; f[1] = f2_; }
void Interaction_2_Body::set_forces(VECTOR& f1_, VECTOR& f2_){  f[0] = &f1_; f[1] = &f2_; }

void Interaction_2_Body::set_charges(double* q1_, double* q2_){  q[0] = q1_; q[1] = q2_; }
void Interaction_2_Body::set_charges(double& q1_, double& q2_){  q[0] = &q1_; q[1] = &q2_; }







void Bond_Interaction::set_functional(std::string f){
/**
                                int_type           functional

      Bond_Interaction            20               Harmonic     (0)
                                                    Quartic     (1)
                                                     Morse      (2)
*/
    if(f=="Harmonic")    { functional = 0;  }
    else if(f=="Quartic"){ functional = 1;  }
    else if(f=="Morse")  { functional = 2;  }
    else{ std::cout<<"Warning: Bond functional "<<f<<" is not implemented\n";  }

}

void Bond_Interaction::set_params(map<std::string,double>& params){

    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="K"){ K = it->second; }
      else if(it->first=="D"){ D = it->second; }
      else if(it->first=="r0"){ r0 = it->second; }
      else if(it->first=="alpha"){ alpha = it->second; }
    }
}

void Bond_Interaction::set_params(boost::python::dict params){

    map<std::string, double> int_params = dict2map(params);
    set_params(int_params);
}


void Bond_Interaction::compute(){
/** 
*/
    if(is_active){

        VECTOR r1, r2;
        if(r[0]==NULL || r[1]==NULL){
            cout<<"Error in Bond_Interaction::compute() : coordinates are not set\nExiting...";
            exit(0);
        }


        if(t[0]==NULL || t[1]==NULL){
            // if the translations are not set, they are not used (this effectively means they are all zeroes)
            cout<<"Warning in Bond_Interaction::compute() : translations are not set. Not using them\n";
            r1 = *r[0];
            r2 = *r[1];
        }
        else{
            r1 = *r[0]; r1 += *t[0];
            r2 = *r[1]; r2 += *t[1];
        }

        VECTOR f1, f2; 
        MATRIX hess(6,6);

        if(functional==0){ energy += Bond_Harmonic(r1,r2,f1,f2,hess,K,r0,2); }
        else if(functional==1){ energy += Bond_Quartic(r1,r2,f1,f2,K,r0); }
        else if(functional==2){ energy += Bond_Morse(r1,r2,f1,f2,D,r0,alpha); }


        if(f[0]==NULL || f[1]==NULL){
            cout<<"Error in Bond_Interaction::compute() : forces are not set\nWill not be assigned...";
        }
        else{

            *f[0] += f1;
            *f[1] += f2;

            MATRIX3x3 tmp;
            tmp.tensor_product(r1 - r2 , f1); 
            stress_at += tmp;
        }

        if(Hess==NULL){
        }
        else{
            add_submatrix(Hess, &hess, Hess_stenc, 1.0);
        }

    }// if is_active
}





void VdW_Interaction::set_functional(std::string f){
/**
                                int_type           functional

      VdW_Interaction             21                   LJ       (0)
                                                  Buffered14_7  (1)
                                                     Morse      (2)


*/
    if(f=="LJ")    { functional = 0;  }
    else if(f=="Buffered14_7"){ functional = 1;  }
    else if(f=="Morse")  { functional = 2;  }
    else{ std::cout<<"Warning: VdW functional "<<f<<" is not implemented\n";  }

}

void VdW_Interaction::set_params(map<std::string,double>& params){

    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="sigma"){ sigma = it->second; }
      else if(it->first=="epsilon"){ epsilon = it->second; }
      else if(it->first=="D"){ D = it->second; }
      else if(it->first=="r0"){ r0 = it->second; }
      else if(it->first=="alpha"){ alpha = it->second; }
      else if(it->first=="scale"){ scale = it->second; }

    }
}
void VdW_Interaction::set_params(boost::python::dict params){

    map<std::string, double> int_params = dict2map(params);
    set_params(int_params);
}


void VdW_Interaction::compute(double R_on, double R_off){
/** 

    R_on       Distance when the switching function starts (and is = 1)
    R_off      Distance when the switching function stops (and is = 0)

*/
    if(is_active){

        //======== Get the positions and translations ============
        VECTOR r1, r2;
        if(r[0]==NULL || r[1]==NULL){
            cout<<"Error in VdW_Interaction::compute() : coordinates are not set\nExiting...";
            exit(0);
        }


        if(t[0]==NULL || t[1]==NULL){
            // if the translations are not set, they are not used (this effectively means they are all zeroes)
            cout<<"Warning in VdW_Interaction::compute() : translations are not set. Not using them\n";
            r1 = *r[0];
            r2 = *r[1];
        }
        else{
            r1 = *r[0]; r1 += *t[0];
            r2 = *r[1]; r2 += *t[1];
        }

        //=============== Switching ==================
        double R_on2 = R_on*R_on;
        double R_off2 = R_off*R_off;
        double SW = 1.0; 
        VECTOR dSW; dSW = 0.0;
        double en = 0.0;

        double d12 = (r1-r2).length2();
        if(d12<=R_off2){
            if(d12>=R_on2){
                SWITCH(r1,r2,R_on,R_off,SW,dSW);
            }
        }else{ SW = 0.0; dSW = 0.0; }


        //========= Run actual calculations ============== 
        VECTOR f1, f2, f12; 

        if(SW>0.0){
            if(functional==0)     { en = Vdw_LJ(r1,r2,f1,f2,sigma,scale*epsilon);     }      
            else if(functional==1){ en = Vdw_Buffered14_7(r1,r2,f1,f2,sigma,scale*epsilon); }
            else if(functional==2){ en = Vdw_Morse(r1,r2,f1,f2,scale*D,r0,alpha); }

            energy += SW*en;

            if(f[0]==NULL || f[1]==NULL){
                cout<<"Warning in VdW_Interaction::compute() : forces are not set\nWill not be assigned...";
            }
            else{
            
                VECTOR f12 = (SW*f1 - en*dSW);
                *f[0] += f12;
                *f[1] -= f12;

                MATRIX3x3 tmp;
                tmp.tensor_product(r1 -r2 , f12);   
                stress_at += tmp;
            }

        }// SW>0.0

    }// if is_active
}

void VdW_Interaction::compute(){   compute(1e+10, 1.1e+10);  }




void Elec_Interaction::set_functional(std::string f){
/**
                                int_type           functional

      Elec_Interaction            22                 Coulomb    (0)



*/
    if(f=="Coulomb")    { functional = 0;  }
    else{ std::cout<<"Warning: Elec functional "<<f<<" is not implemented\n";  }

}


void Elec_Interaction::set_params(map<std::string,double>& params){

    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="xi1"){ xi1 = it->second; }
      else if(it->first=="xi2"){ xi2 = it->second; }
      else if(it->first=="eps"){ eps = it->second; }
      else if(it->first=="J"){ J = it->second; }
      else if(it->first=="delta"){ delta = it->second; }
      else if(it->first=="scale"){ scale = it->second; }
    }
}

void Elec_Interaction::set_params(boost::python::dict params){

    map<std::string, double> int_params = dict2map(params);
    set_params(int_params);
}

void Elec_Interaction::compute(double R_on, double R_off){
/** 

    R_on       Distance when the switching function starts (and is = 1)
    R_off      Distance when the switching function stops (and is = 0)

*/
    if(is_active){

        //======== Get the positions and translations ============
        VECTOR r1, r2;
        if(r[0]==NULL || r[1]==NULL){
            cout<<"Error in Elec_Interaction::compute() : coordinates are not set\nExiting...";
            exit(0);
        }


        if(t[0]==NULL || t[1]==NULL){
            // if the translations are not set, they are not used (this effectively means they are all zeroes)
            cout<<"Warning in Elec_Interaction::compute() : translations are not set. Not using them\n";
            r1 = *r[0];
            r2 = *r[1];
        }
        else{
            r1 = *r[0]; r1 += *t[0];
            r2 = *r[1]; r2 += *t[1];
        }

        double q1, q2;
        if(q[0]==NULL || q[1]==NULL){
            cout<<"Warning in Elec_Interaction::compute() : charges are not set\nUsing zero charges...";
        }
        else{ 
            q1 = *q[0]; q2 = *q[1]; 
        }

        //=============== Switching ==================
        double R_on2 = R_on*R_on;
        double R_off2 = R_off*R_off;
        double SW = 1.0; 
        VECTOR dSW; dSW = 0.0;
        double en = 0.0;

        double d12 = (r1-r2).length2();
        if(d12<=R_off2){
            if(d12>=R_on2){
                SWITCH(r1,r2,R_on,R_off,SW,dSW);
            }
        }else{ SW = 0.0; dSW = 0.0; }



        //========= Run actual calculations ============== 
        VECTOR f1, f2, f12; 

        if(SW>0.0){
            if(functional==0){ en = Elec_Coulomb(r1,r2,f1,f2,q1,q2,eps,delta); }

            energy += SW*en;

            if(f[0]==NULL || f[1]==NULL){
                cout<<"Error in Elec_Interaction::compute() : forces are not set\nWill not be assigned...";
            }
            else{
            
                VECTOR f12 = (SW*f1 - en*dSW);
                *f[0] += f12;
                *f[1] -= f12;

                MATRIX3x3 tmp;
                tmp.tensor_product(r1 -r2 , f12);   
                stress_at += tmp;
            }

        }// SW>0.0

    }// if is_active
}

void Elec_Interaction::compute(){   compute(1e+10, 1.1e+10);  }




}// namespace libhamiltonian_mm
}// namespace libatomistic
}// liblibra
