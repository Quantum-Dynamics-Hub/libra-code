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
  \file Interactions_4_Body.cpp
  \brief The file implements functions for molecular-mechanical Hamiltonian calculations as well
  as the classes for organizing such computations in an object-oriented way: basic functionality of the class
*/

#include "Interactions_4_Body.h"

/// liblibra namespace
namespace liblibra{


namespace libatomistic{

/// libhamiltonian_mm namespace
namespace libhamiltonian_mm{



void Interaction_4_Body::set_coords(VECTOR* r1_, VECTOR* r2_, VECTOR* r3_, VECTOR* r4_){  r[0] = r1_; r[1] = r2_; r[2] = r3_; r[3] = r4_;}
void Interaction_4_Body::set_coords(VECTOR& r1_, VECTOR& r2_, VECTOR& r3_, VECTOR& r4_){  r[0] = &r1_; r[1] = &r2_; r[2] = &r3_; r[3] = &r4_;}

void Interaction_4_Body::set_transl(VECTOR* t1_, VECTOR* t2_, VECTOR* t3_, VECTOR* t4_){  t[0] = t1_; t[1] = t2_; t[2] = t3_; t[3] = t4_;}
void Interaction_4_Body::set_transl(VECTOR& t1_, VECTOR& t2_, VECTOR& t3_, VECTOR& t4_){  t[0] = &t1_; t[1] = &t2_; t[2] = &t3_; t[3] = &t4_;}

void Interaction_4_Body::set_forces(VECTOR* f1_, VECTOR* f2_, VECTOR* f3_, VECTOR* f4_){  f[0] = f1_; f[1] = f2_; f[2] = f3_; f[3] = f4_;}
void Interaction_4_Body::set_forces(VECTOR& f1_, VECTOR& f2_, VECTOR& f3_, VECTOR& f4_){  f[0] = &f1_; f[2] = &f2_; f[2] = &f3_; f[3] = &f4_;}

void Interaction_4_Body::set_charges(double* q1_, double* q2_, double* q3_, double* q4_){  q[0] = q1_; q[1] = q2_; q[2] = q3_; q[3] = q4_;}
void Interaction_4_Body::set_charges(double& q1_, double& q2_, double& q3_, double& q4_){  q[0] = &q1_; q[1] = &q2_; q[2] = &q3_; q[3] = &q4_;}


void Dihedral_Interaction::set_functional(std::string f){
/**
                                int_type           functional

      Dihedral_Interaction        40               General0     (0)
                                                   General1     (1)
                                                   General2     (2)
                                                   General3     (3)
                                                   Fourier0     (4)
                                                   Fourier1     (5)

*/

    if(f=="General0")       { functional = 0;  }
    else if(f=="General1")  { functional = 1;  }
    else if(f=="General2")  { functional = 2;  }
    else if(f=="General3")  { functional = 3;  }
    else if(f=="Fourier0")  { functional = 4;  }
    else if(f=="Fourier1")  { functional = 5;  }
    else{ std::cout<<"Warning: Dihedral functional "<<f<<" is not implemented\n"; }

}


void Dihedral_Interaction::set_params(map<std::string,double>& params){

    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="Vphi"){ Vphi = it->second; }
      else if(it->first=="phi0"){ phi0 = it->second; }
      else if(it->first=="Vphi1"){ Vphi1 = it->second; }
      else if(it->first=="Vphi2"){ Vphi2 = it->second; }
      else if(it->first=="Vphi3"){ Vphi3 = it->second; }
      else if(it->first=="opt"){ opt = (int)it->second; }
      else if(it->first=="n"){ n = (int)it->second; }
    }

}


void Dihedral_Interaction::compute(){
/** 
*/
    if(is_active){

        VECTOR r1, r2, r3, r4;
        if(r[0]==NULL || r[1]==NULL || r[2]==NULL || r[3]==NULL ){
            cout<<"Error in Dihedral_Interaction::compute() : coordinates are not set\nExiting...";
            exit(0);
        }


        if(t[0]==NULL || t[1]==NULL || t[2]==NULL || t[3]==NULL ){
            // if the translations are not set, they are not used (this effectively means they are all zeroes)
            cout<<"Warning in Dihedral_Interaction::compute() : translations are not set. Not using them\n";
            r1 = *r[0];
            r2 = *r[1];
            r3 = *r[2];
            r4 = *r[3];
        }
        else{
            r1 = *r[0]; r1 += *t[0];
            r2 = *r[1]; r2 += *t[1];
            r3 = *r[1]; r3 += *t[2];
            r4 = *r[1]; r4 += *t[3];
        }

        VECTOR f1, f2, f3, f4; 

        if(functional==0){ energy += Dihedral_General(r1,r2,r3,r4,f1,f2,f3,f4,Vphi,phi0,n,opt); }
        else if(functional==1){ energy += Dihedral_Fourier(r1,r2,r3,r4,f1,f2,f3,f4,Vphi1,Vphi2,Vphi3,opt); }


        if(f[0]==NULL || f[1]==NULL || f[2]==NULL || f[3]==NULL ){
            cout<<"Warning in Dihedral_Interaction::compute() : forces are not set\nWill not be assigned...";
        }
        else{

            *f[0] += f1;
            *f[1] += f2;
            *f[2] += f3;
            *f[3] += f4;

            MATRIX3x3 tmp;
            tmp.tensor_product(r1 , f1);   stress_at += tmp;
            tmp.tensor_product(r2 , f2);   stress_at += tmp;
            tmp.tensor_product(r3 , f3);   stress_at += tmp;
            tmp.tensor_product(r4 , f4);   stress_at += tmp;
        }

    }// if is_active
}



void OOP_Interaction::set_functional(std::string f){
/**
                                int_type           functional

      OOP_Interaction             41               Fourier      (0)
                                                   Wilson       (1)
                                                   Harmonic     (2)

*/

    if(f=="Fourier")        { functional = 0;  }
    else if(f=="Wilson")    { functional = 1;  }
    else if(f=="Harmonic")  { functional = 2;  }
    else{ std::cout<<"Warning: OOP functional "<<f<<" is not implemented\n"; }

}

void OOP_Interaction::set_params(map<std::string,double>& params){

    for(map<std::string,double>::iterator it=params.begin();it!=params.end();it++){
      if(it->first=="K"){ K = it->second; }
      else if(it->first=="xi_0"){ xi_0 = it->second; }
      else if(it->first=="C0"){ C0 = it->second; }
      else if(it->first=="C1"){ C1 = it->second; }
      else if(it->first=="C2"){ C2 = it->second; }
      else if(it->first=="opt"){ opt = (int)it->second; }
    }


}


void OOP_Interaction::compute(){
/** 
*/
    if(is_active){

        VECTOR r1, r2, r3, r4;
        if(r[0]==NULL || r[1]==NULL || r[2]==NULL || r[3]==NULL ){
            cout<<"Error in OOP_Interaction::compute() : coordinates are not set\nExiting...";
            exit(0);
        }


        if(t[0]==NULL || t[1]==NULL || t[2]==NULL || t[3]==NULL ){
            // if the translations are not set, they are not used (this effectively means they are all zeroes)
            cout<<"Warning in OOP_Interaction::compute() : translations are not set. Not using them\n";
            r1 = *r[0];
            r2 = *r[1];
            r3 = *r[2];
            r4 = *r[3];
        }
        else{
            r1 = *r[0]; r1 += *t[0];
            r2 = *r[1]; r2 += *t[1];
            r3 = *r[1]; r3 += *t[2];
            r4 = *r[1]; r4 += *t[3];
        }

        VECTOR f1, f2, f3, f4; 

        if(functional==0){ energy += OOP_Fourier(r1,r2,r3,r4,f1,f2,f3,f4,K,C0,C1,C2,opt); }
        else if(functional==1){ energy += OOP_Wilson(r1,r2,r3,r4,f1,f2,f3,f4,K,xi_0); }
        else if(functional==2){ energy += OOP_Harmonic(r1,r2,r3,r4,f1,f2,f3,f4,K); }



        if(f[0]==NULL || f[1]==NULL || f[2]==NULL || f[3]==NULL ){
            cout<<"Warning in Dihedral_Interaction::compute() : forces are not set\nWill not be assigned...";
        }
        else{

            *f[0] += f1;
            *f[1] += f2;
            *f[2] += f3;
            *f[3] += f4;

            MATRIX3x3 tmp;
            tmp.tensor_product(r1 , f1);   stress_at += tmp;
            tmp.tensor_product(r2 , f2);   stress_at += tmp;
            tmp.tensor_product(r3 , f3);   stress_at += tmp;
            tmp.tensor_product(r4 , f4);   stress_at += tmp;
        }

    }// if is_active
}






}// namespace libhamiltonian_mm
}// namespace libatomistic
}// liblibra
