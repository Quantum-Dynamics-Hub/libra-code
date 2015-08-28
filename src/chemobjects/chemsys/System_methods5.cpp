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

#include "System.h"

namespace libchemobjects{
namespace libchemsys{

void System::extract_atomic_q(vector<double>& q){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(q.size() != 3.0 * Number_of_atoms){
    cout<<"Error: in System::extract_atomic_q: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_atoms;i++){

    q[3*i+0] = Atoms[i].Atom_RB.rb_cm.x; 
    q[3*i+1] = Atoms[i].Atom_RB.rb_cm.y; 
    q[3*i+2] = Atoms[i].Atom_RB.rb_cm.z; 

  }// for i


}

void System::set_atomic_q(vector<double>& q){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(q.size() != 3.0 * Number_of_atoms){
    cout<<"Error: in System::set_atomic_q: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_atoms;i++){

    VECTOR q_(q[3*i+0], q[3*i+1], q[3*i+2]);
    Atoms[i].Atom_RB.set_position(q_);

  }// for i


}


void System::extract_atomic_p(vector<double>& p){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(p.size() != 3.0 * Number_of_atoms){
    cout<<"Error: in System::extract_atomic_p: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_atoms;i++){

    p[3*i+0] = Atoms[i].Atom_RB.rb_p.x; 
    p[3*i+1] = Atoms[i].Atom_RB.rb_p.y; 
    p[3*i+2] = Atoms[i].Atom_RB.rb_p.z; 

  }// for i



}

void System::set_atomic_p(vector<double>& p){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(p.size() != 3.0 * Number_of_atoms){
    cout<<"Error: in System::set_atomic_p: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_atoms;i++){

    VECTOR p_(p[3*i+0], p[3*i+1], p[3*i+2]);
    Atoms[i].Atom_RB.set_momentum(p_);

  }// for i

}

void System::extract_atomic_v(vector<double>& v){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(v.size() != 3.0 * Number_of_atoms){
    cout<<"Error: in System::extract_atomic_v: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_atoms;i++){

    v[3*i+0] = Atoms[i].Atom_RB.rb_v.x; 
    v[3*i+1] = Atoms[i].Atom_RB.rb_v.y; 
    v[3*i+2] = Atoms[i].Atom_RB.rb_v.z; 

  }// for i

}

void System::set_atomic_v(vector<double>& v){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(v.size() != 3.0 * Number_of_atoms){
    cout<<"Error: in System::set_atomic_v: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_atoms;i++){

    VECTOR v_(v[3*i+0], v[3*i+1], v[3*i+2]);
    Atoms[i].Atom_RB.set_velocity(v_);

  }// for i

}



void System::extract_atomic_f(vector<double>& f){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(f.size() != 3.0 * Number_of_atoms){
    cout<<"Error: in System::extract_atomic_f: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_atoms;i++){

    f[3*i+0] = Atoms[i].Atom_RB.rb_force.x; 
    f[3*i+1] = Atoms[i].Atom_RB.rb_force.y; 
    f[3*i+2] = Atoms[i].Atom_RB.rb_force.z; 

  }// for i

}

void System::set_atomic_f(vector<double>& f){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(f.size() != 3.0 * Number_of_atoms){
    cout<<"Error: in System::set_atomic_f: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_atoms;i++){

    VECTOR f_(f[3*i+0], f[3*i+1], f[3*i+2]);
    Atoms[i].Atom_RB.set_force(f_);

  }// for i

}


void System::extract_atomic_mass(vector<double>& mass){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(mass.size() != 3.0 * Number_of_atoms){
    cout<<"Error: in System::extract_atomic_mass: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_atoms;i++){

    mass[3*i+0] = Atoms[i].Atom_RB.rb_mass; 
    mass[3*i+1] = Atoms[i].Atom_RB.rb_mass; 
    mass[3*i+2] = Atoms[i].Atom_RB.rb_mass; 

  }// for i

}

void System::set_atomic_mass(vector<double>& mass){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(mass.size() != 3.0 * Number_of_atoms){
    cout<<"Error: in System::set_atomic_mass: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_atoms;i++){

    Atoms[i].Atom_RB.set_mass(mass[3*i+0]); // assume mass[3*i+0] = mass[3*i+1] = mass[3*i+2]

  }// for i

}



void System::extract_fragment_q(vector<double>& q){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(q.size() != 3.0 * Number_of_fragments){
    cout<<"Error: in System::extract_fragment_q: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_fragments;i++){

    q[3*i+0] = Fragments[i].Group_RB.rb_cm.x; 
    q[3*i+1] = Fragments[i].Group_RB.rb_cm.y; 
    q[3*i+2] = Fragments[i].Group_RB.rb_cm.z; 

  }// for i


}

void System::set_fragment_q(vector<double>& q){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(q.size() != 3.0 * Number_of_fragments){
    cout<<"Error: in System::set_fragment_q: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_fragments;i++){

    VECTOR q_(q[3*i+0], q[3*i+1], q[3*i+2]);
    Fragments[i].Group_RB.set_position(q_);

  }// for i


}



void System::extract_fragment_p(vector<double>& p){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(p.size() != 3.0 * Number_of_fragments){
    cout<<"Error: in System::extract_fragment_p: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_fragments;i++){

    p[3*i+0] = Fragments[i].Group_RB.rb_p.x; 
    p[3*i+1] = Fragments[i].Group_RB.rb_p.y; 
    p[3*i+2] = Fragments[i].Group_RB.rb_p.z; 

  }// for i


}

void System::set_fragment_p(vector<double>& p){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(p.size() != 3.0 * Number_of_fragments){
    cout<<"Error: in System::set_fragment_p: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_fragments;i++){

    VECTOR p_(p[3*i+0], p[3*i+1], p[3*i+2]);
    Fragments[i].Group_RB.set_momentum(p_);

  }// for i


}

void System::extract_fragment_v(vector<double>& v){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(v.size() != 3.0 * Number_of_fragments){
    cout<<"Error: in System::extract_fragment_v: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_fragments;i++){

    v[3*i+0] = Fragments[i].Group_RB.rb_v.x; 
    v[3*i+1] = Fragments[i].Group_RB.rb_v.y; 
    v[3*i+2] = Fragments[i].Group_RB.rb_v.z; 

  }// for i

}

void System::set_fragment_v(vector<double>& v){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(v.size() != 3.0 * Number_of_fragments){
    cout<<"Error: in System::set_fragments_v: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_fragments;i++){

    VECTOR v_(v[3*i+0], v[3*i+1], v[3*i+2]);
    Fragments[i].Group_RB.set_velocity(v_);

  }// for i

}




void System::extract_fragment_f(vector<double>& f){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(f.size() != 3.0 * Number_of_fragments){
    cout<<"Error: in System::extract_fragment_f: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_fragments;i++){

    f[3*i+0] = Fragments[i].Group_RB.rb_force.x; 
    f[3*i+1] = Fragments[i].Group_RB.rb_force.y; 
    f[3*i+2] = Fragments[i].Group_RB.rb_force.z; 

  }// for i


}

void System::set_fragment_f(vector<double>& f){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(f.size() != 3.0 * Number_of_fragments){
    cout<<"Error: in System::set_fragment_f: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_fragments;i++){

    VECTOR f_(f[3*i+0], f[3*i+1], f[3*i+2]);
    Fragments[i].Group_RB.set_force(f_);

  }// for i


}


void System::extract_fragment_mass(vector<double>& mass){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(mass.size() != 3.0 * Number_of_fragments){
    cout<<"Error: in System::extract_fragment_q: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_fragments;i++){

    mass[3*i+0] = Fragments[i].Group_RB.rb_mass; 
    mass[3*i+1] = Fragments[i].Group_RB.rb_mass; 
    mass[3*i+2] = Fragments[i].Group_RB.rb_mass; 

  }// for i


}


void System::set_fragment_mass(vector<double>& mass){

  // Simple data converter - no object creation (nor reduction of dimension)

  if(mass.size() != 3.0 * Number_of_fragments){
    cout<<"Error: in System::extract_fragment_q: dimensions of System and vector<double> objects are not compatible\n";  
    exit(0);
  }
  
  for(int i=0;i<Number_of_fragments;i++){

    Fragments[i].Group_RB.set_mass(mass[3*i+0]); // assume mass[3*i+0] = mass[3*i+1] = mass[3*i+2]

  }// for i


}




}// namespace libchemsys
}// namespace libchemobjects
