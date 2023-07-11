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
/**
  \file Nuclear.cpp
  \brief The file implements Nuclear class and functions for propagation of nuclear DOF
    
*/

#include "Nuclear.h"

/// liblibra namespace
namespace liblibra{


/// libdyn namespace
namespace libdyn{

/// libnuclear namespace
namespace libnuclear{


Nuclear::Nuclear(){  
/** \brief Default constructor

  Does not do anything special
*/
 nnucl = 0; 
}

Nuclear::Nuclear(int _nnucl){  
/** 
  \brief Constructor with parameters

  Allocates memory and initializes variables to the default values
  Default masses = 2000 (a.u., so approximately the mass of hydrogen atom)
  Default coordinates, momenta, forces = 0.0

  \param[in] _nnucl The number of nuclear DOF
*/

    nnucl = _nnucl;
    mass = vector<double>(nnucl,2000.0); 
    q = vector<double>(nnucl,0.0);
    p = vector<double>(nnucl,0.0);
    f = vector<double>(nnucl,0.0);
    ctyp = vector<int>(nnucl,0);  
}

Nuclear::Nuclear(const Nuclear& ob){
/** 
  \brief Copy constructor

  Allocates memory and initializes variables using existing Nuclear object

  \param[in] ob The object used to initialize the present object
*/

  nnucl = ob.nnucl;
  mass = ob.mass;
  q = ob.q; 
  p = ob.p;
  f = ob.f;
  ctyp = ob.ctyp;
}

Nuclear::~Nuclear(){
/** 
  \brief Destructor

*/

  if(mass.size()>0){  mass.clear(); }
  if(q.size()>0){  q.clear(); }
  if(p.size()>0){  p.clear(); }
  if(f.size()>0){  f.clear(); }
  if(ctyp.size()>0){  ctyp.clear(); }
}


void Nuclear::propagate_p(int i,double dt){ 
/** 
  \brief Propagate momentum for a selected DOF

  p[i] --> p[i] + dt*f[i]

  \param[in] i The index of the selected DOF
  \param[in] dt The time step
*/

  p[i] += dt*f[i]; 
}

void Nuclear::propagate_p(double dt){
/** 
  \brief Propagate momenta of all DOF

  \param[in] dt The time step
*/

  for(int i=0;i<nnucl;i++){  p[i] += dt*f[i];   }  
}

void Nuclear::propagate_p(double dt,vector<int>& active){  
/** 
  \brief Propagate only selected (active) DOF

  \param[in] dt The time step
  \param[in] active The vector of N integers: active[i] is just the i-th active (participating in propagation) DOF
*/

  int sz = active.size();  for(int i=0;i<sz;i++){ int I = active[i]; p[I] += dt*f[I];   }
}

void Nuclear::scale_p(int i,double scl){
/** 
  \brief Scale momentum for selected DOF

  p[i] -->  p[i] * slc

  \param[in] i The index of the selected DOF
  \param[in] sc The scaling parameter
*/

  p[i] *= scl;  
}

void Nuclear::scale_p(double scl){
/** 
  \brief Scale momenta for all DOF

  \param[in] sc The scaling parameter
*/

  for(int i=0;i<nnucl;i++){  p[i] *= scl;   }  
}

void Nuclear::scale_p(double scl,vector<int>& active){  
/** 
  \brief Scale momenta for active DOF

  \param[in] sc The scaling parameter
  \param[in] active The vector of N integers: active[i] is just the i-th active (participating in propagation) DOF
*/

  int sz = active.size(); for(int i=0;i<sz;i++){ int I = active[i]; p[I] *= scl;   } 
}


void Nuclear::propagate_q(int i,double dt){
/** 
  \brief Propagate coordinate for selected DOF as:

  q[i] -->  q[i] + dt*p[i]/mass[i]

  \param[in] i The index of the selected DOF
  \param[in] dt The integration time step
*/

  q[i] += dt*p[i]/mass[i];  
}

void Nuclear::propagate_q(double dt){
/** 
  \brief Propagate coordinate for all DOF

  \param[in] dt The integration time step
*/

  for(int i=0;i<nnucl;i++){  q[i] += dt*p[i]/mass[i];   } 
}

void Nuclear::propagate_q(double dt,vector<int>& active){  
/** 
  \brief Propagate coordinate for all active DOF

  \param[in] dt The integration time step
  \param[in] active The vector of N integers: active[i] is just the i-th active (participating in propagation) DOF
*/

  int sz = active.size();  for(int i=0;i<sz;i++){ int I = active[i]; q[I] += dt*p[I]/mass[I];   }
}


void Nuclear::scale_q(int i,double scl){
/** 
  \brief Scale coordinate for selected DOF

  q[i] -->  q[i] * slc

  \param[in] i The index of the selected DOF
  \param[in] sc The scaling parameter
*/

  q[i] *= scl;  
}

void Nuclear::scale_q(double scl){
/** 
  \brief Scale coordinates for all DOF

  \param[in] sc The scaling parameter
*/

  for(int i=0;i<nnucl;i++){  q[i] *= scl;   }  
}

void Nuclear::scale_q(double scl,vector<int>& active){  
/** 
  \brief Scale coordinates for active DOF

  \param[in] sc The scaling parameter
  \param[in] active The vector of N integers: active[i] is just the i-th active (participating in propagation) DOF
*/

  int sz = active.size(); for(int i=0;i<sz;i++){ int I = active[i]; q[I] *= scl;   } 
}




}// namespace libnuclear
}// namespace libdyn
}// liblibra

