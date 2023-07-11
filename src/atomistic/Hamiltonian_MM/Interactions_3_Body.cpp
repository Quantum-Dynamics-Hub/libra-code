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
  \file Interactions_3_Body.cpp
  \brief The file implements functions for molecular-mechanical Hamiltonian calculations as well
  as the classes for organizing such computations in an object-oriented way: basic functionality of the class
*/

#include "Interactions_3_Body.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_mm namespace
namespace libhamiltonian_mm{




void Interaction_3_Body::set_coords(VECTOR* r1_, VECTOR* r2_, VECTOR* r3_){  r[0] = r1_; r[1] = r2_; r[2] = r3_;}
void Interaction_3_Body::set_coords(VECTOR& r1_, VECTOR& r2_, VECTOR& r3_){  r[0] = &r1_; r[1] = &r2_; r[2] = &r3_;}

void Interaction_3_Body::set_transl(VECTOR* t1_, VECTOR* t2_, VECTOR* t3_){  t[0] = t1_; t[1] = t2_; t[2] = t3_;}
void Interaction_3_Body::set_transl(VECTOR& t1_, VECTOR& t2_, VECTOR& t3_){  t[0] = &t1_; t[1] = &t2_; t[2] = &t3_;}

void Interaction_3_Body::set_forces(VECTOR* f1_, VECTOR* f2_, VECTOR* f3_){  f[0] = f1_; f[1] = f2_; f[2] = f3_;}
void Interaction_3_Body::set_forces(VECTOR& f1_, VECTOR& f2_, VECTOR& f3_){  f[0] = &f1_; f[2] = &f2_; f[2] = &f3_;}

void Interaction_3_Body::set_charges(double* q1_, double* q2_, double* q3_){  q[0] = q1_; q[1] = q2_; q[2] = q3_;}
void Interaction_3_Body::set_charges(double& q1_, double& q2_, double& q3_){  q[0] = &q1_; q[1] = &q2_; q[2] = &q3_;}


void Angle_Interaction::set_functional(std::string f){
/**
                                int_type           functional

      Angle_Interaction           30                Harmonic    (0)
                                                    Fourier     (1)
                                                Fourier_General (2)
                                                Fourier_Special (3)
                                                  Harmonic_Cos  (4)
                                            Harmonic_Cos_General(5)
                                                      Cubic     (6)

*/

    if(f=="Harmonic")                 { functional = 0;  }
    else if(f=="Fourier")             { functional = 1;  }
    else if(f=="Fourier_General")     { functional = 2;  }
    else if(f=="Fourier_Special")     { functional = 3;  }
    else if(f=="Harmonic_Cos")        { functional = 4;  }
    else if(f=="Harmonic_Cos_General"){ functional = 5;  }
    else if(f=="Cubic")               { functional = 6;  }
    else{ std::cout<<"Warning: Angle functional "<<f<<" is not implemented\n"; }

}


void Angle_Interaction::set_params(map<std::string,double>& params){

    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="k_theta"){ k_theta = it->second; }
      else if(it->first=="theta_0"){ theta_0 = it->second; }
      else if(it->first=="cos_theta_0"){ cos_theta_0 = it->second; }
      else if(it->first=="C0"){ C0 = it->second; }
      else if(it->first=="C1"){ C1 = it->second; }
      else if(it->first=="C2"){ C2 = it->second; }
      else if(it->first=="coordination"){ coordination = (int)it->second; }
    }

}


void Angle_Interaction::compute(){
/** 
*/
    if(is_active){

        VECTOR r1, r2, r3;
        if(r[0]==NULL || r[1]==NULL || r[2]==NULL){
            cout<<"Error in Angle_Interaction::compute() : coordinates are not set\nExiting...";
            exit(0);
        }


        if(t[0]==NULL || t[1]==NULL || t[2]==NULL){
            // if the translations are not set, they are not used (this effectively means they are all zeroes)
            cout<<"Warning in Angle_Interaction::compute() : translations are not set. Not using them\n";
            r1 = *r[0];
            r2 = *r[1];
            r3 = *r[2];
        }
        else{
            r1 = *r[0]; r1 += *t[0];
            r2 = *r[1]; r2 += *t[1];
            r3 = *r[1]; r3 += *t[2];
        }

        VECTOR f1, f2, f3; 

        if(functional==0){ energy += Angle_Harmonic(r1,r2,r3,f1,f2,f3,k_theta,theta_0); }
        else if(functional==1){ energy += Angle_Fourier(r1,r2,r3,f1,f2,f3,k_theta,C0, C1,C2,coordination); }
        else if(functional==2){ energy += Angle_Fourier_General(r1,r2,r3,f1,f2,f3,k_theta,C0,C1,C2); }
        else if(functional==3){ energy += Angle_Fourier_Special(r1,r2,r3,f1,f2,f3,k_theta,coordination); }
        else if(functional==4){ energy += Angle_Harmonic_Cos(r1,r2,r3,f1,f2,f3,k_theta,cos_theta_0,coordination); }
        else if(functional==5){ energy += Angle_Harmonic_Cos_General(r1,r2,r3,f1,f2,f3,k_theta,cos_theta_0); }
        else if(functional==6){ energy += Angle_Cubic(r1,r2,r3,f1,f2,f3,k_theta,theta_0); }


        if(f[0]==NULL || f[1]==NULL || f[2]==NULL){
            cout<<"Warning in Angle_Interaction::compute() : forces are not set\nWill not be assigned...";
        }
        else{

            *f[0] += f1;
            *f[1] += f2;
            *f[2] += f3;

            MATRIX3x3 tmp;
            tmp.tensor_product(r1 , f1);   stress_at += tmp;
            tmp.tensor_product(r2 , f2);   stress_at += tmp;
            tmp.tensor_product(r3 , f3);   stress_at += tmp;
        }

    }// if is_active
}




}// namespace libhamiltonian_mm
}// namespace libatomistic
}// liblibra
