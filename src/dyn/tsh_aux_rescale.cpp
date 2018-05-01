/*********************************************************************************
* Copyright (C) 2015-2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file tsh_aux_rescale.cpp
  \brief The file implements the momenta rescaling procedures used in the TSH algorithms
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;
using namespace libhamiltonian;
namespace bp = boost::python;


/// libdyn namespace
namespace libdyn{



int rescale_velocities_adiabatic
(MATRIX& p, MATRIX& invM, CMATRIX* ham_adi, vector<CMATRIX*>& dc1_adi, int new_st,int old_st, int do_reverse){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another
  \param[in,out] p [ndof x 1] matrix of nuclear DOF momenta
  \param[in] invM [ndox x 1] matrix of inverse masses for all nuclear DOFs
  \param[in,out] ham_adi [nadi x nadi] Hamiltonian matrix in the adiabatic representation
  \param[in,out] dc1_adi ntraj x [nadi x nadi] List of pointer to derivative coupling matrices for all nuclear DOFs
  \param[in] new_st The proposed new state
  \param[in] old_st The original state
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  Return: the function return the final state: the proposed if the velocity rescaling is successfull, or the initial, if not.

  This version implies that the adiabatic representation is used

*/
  int ndof = p.n_rows;

  int st;

  if(new_st!=old_st){

//    cout<<"old_st = "<<old_st<<"  new_st = "<<new_st<<endl;

    // According to Fabiano
    double a_ij = 0.0;
    double b_ij = 0.0;

    for(int k=0;k<ndof;k++){

      double D = dc1_adi[k]->get(old_st,new_st).real(); // derivative coupling w.r.t. nuclear DOF k

      a_ij += 0.5*(D*D*invM.get(k,0)); 
      b_ij += (D*p.get(k,0))*invM.get(k,0);
    }
    double det = b_ij*b_ij + 4.0*a_ij*(ham_adi->get(old_st,old_st).real() - ham_adi->get(new_st,new_st).real());

    // Calculate the scaling factor and new state
    double gamma_ij = 0.0;

    if(det<0.0){

      if(fabs(a_ij)>1e-100){  // only consider reversals, if the couplings are sizable
        if(do_reverse){     gamma_ij = b_ij / a_ij;}
         else{ gamma_ij = 0.0;  }
      }
      else{  gamma_ij = 0.0;    } // don't consider reversal, if the couplings are too small

      st = old_st; // # hop does not occur - frustrated hop

    }
    else{
      if(fabs(a_ij)>1e-100){ // only compute the rescaling factor and do the hop, if the couplings are sizable

        if(b_ij<0){ gamma_ij = 0.5*(b_ij + sqrt(det))/a_ij; }
        else{       gamma_ij = 0.5*(b_ij - sqrt(det))/a_ij; }

        st = new_st;
      } 
      else{   // otherwise - don't 
        gamma_ij = 0.0; 
        st = old_st;
      } 

    }
//    cout<<"Velocity rescaling: factor = "<<gamma_ij<<endl;
//    cout<<"old_st = "<<old_st<<"  new_st = "<<new_st<<endl;

    //Rescale velocities and do the hop
    for(int k=0;k<ndof;k++){ 
      double D = dc1_adi[k]->get(old_st,new_st).real(); 
      p.add(k, 0, - gamma_ij * D); 
    }

  }
  else{ st = old_st; }

  return st;

} // rescale velocities adiabatic


int rescale_velocities_adiabatic
(MATRIX& p, MATRIX& invM, CMATRIX& ham_adi, vector<CMATRIX>& dc1_adi, int new_st,int old_st, int do_reverse){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another
  \param[in,out] p [ndof x 1] matrix of nuclear DOF momenta
  \param[in] invM [ndox x 1] matrix of inverse masses for all nuclear DOFs
  \param[in,out] ham_adi [nadi x nadi] Hamiltonian matrix in the adiabatic representation
  \param[in,out] dc1_adi ntraj x [nadi x nadi] List of pointer to derivative coupling matrices for all nuclear DOFs
  \param[in] new_st The proposed new state
  \param[in] old_st The original state
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  Return: the function return the final state: the proposed if the velocity rescaling is successfull, or the initial, if not.

  This version implies that the adiabatic representation is used

*/
  int ndof = p.n_rows;

  int st;

  if(new_st!=old_st){

    // According to Fabiano
    double a_ij = 0.0;
    double b_ij = 0.0;

    for(int k=0;k<ndof;k++){

      double D = dc1_adi[k].get(old_st,new_st).real(); // derivative coupling w.r.t. nuclear DOF k

      a_ij += 0.5*(D*D*invM.get(k,0)); 
      b_ij += (D*p.get(k,0))*invM.get(k,0);
    }
    double det = b_ij*b_ij + 4.0*a_ij*(ham_adi.get(old_st,old_st).real() - ham_adi.get(new_st,new_st).real());

    // Calculate the scaling factor and new state
    double gamma_ij = 0.0;

    if(det<0.0){

      if(fabs(a_ij)>1e-100){  // only consider reversals, if the couplings are sizable
        if(do_reverse){     gamma_ij = b_ij / a_ij;}
         else{ gamma_ij = 0.0;  }
      }
      else{  gamma_ij = 0.0;    } // don't consider reversal, if the couplings are too small

      st = old_st; // # hop does not occur - frustrated hop

    }
    else{
      if(fabs(a_ij)>1e-100){ // only compute the rescaling factor and do the hop, if the couplings are sizable

        if(b_ij<0){ gamma_ij = 0.5*(b_ij + sqrt(det))/a_ij; }
        else{       gamma_ij = 0.5*(b_ij - sqrt(det))/a_ij; }

        st = new_st;
      } 
      else{   // otherwise - don't 
        gamma_ij = 0.0; 
        st = old_st;
      } 

    }

    //Rescale velocities and do the hop
    for(int k=0;k<ndof;k++){ 
      double D = dc1_adi[k].get(old_st,new_st).real(); 
      p.add(k, 0, - gamma_ij * D); 
    }

  }
  else{ st = old_st; }

  return st;

} // rescale velocities adiabatic




int rescale_velocities_adiabatic
(MATRIX& p, MATRIX& invM, nHamiltonian* ham, int new_st,int old_st, int do_reverse){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another and do the
  rescaling if the hop is allowed

  \param[in,out] p [ndof x 1] matrix of the nuclear DOF momenta
  \param[in] invM [ndof x 1] matrix with the inverse masses of nuclear DOFs
  \param[in] ham The Hamiltonian class object that takes care of all energy-related calculations
  \param[in] new_st The proposed new state
  \param[in] old_st The original state
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  Return: the function return the final state: the proposed if the velocity rescaling is successfull, or the initial, if not.

  This version implies that the adiabatic representation is used
*/


  return rescale_velocities_adiabatic(p, invM, ham->ham_adi, ham->dc1_adi, new_st, old_st, do_reverse);

}

int rescale_velocities_adiabatic
(MATRIX& p, MATRIX& invM, nHamiltonian& ham, int new_st,int old_st, int do_reverse){
/**
  \brief See the description of the
  int rescale_velocities_adiabatic
  (MATRIX& p, MATRIX& invM, nHamiltonian* ham, int new_st,int old_st, int do_reverse) function
*/


  return rescale_velocities_adiabatic(p, invM, ham.ham_adi, ham.dc1_adi, new_st, old_st, do_reverse);

}




int rescale_velocities_adiabatic(vector<double>& p, vector<double>& masses, 
                                 CMATRIX& ham_adi, vector<CMATRIX>& dc1_adi,
                                 int new_st,int old_st, int do_reverse){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another and do the
  rescaling if the hop is allowed

  \param[in,out] p Momenta of nuclear DOFs
  \param[in,out] masses Masses of nuclear DOFs
  \param[in,out] ham_adi Energies of all adiabatic states
  \param[in,out] dc1_adi Derivative coupling matrices

  \param[in] new_st The proposed new state
  \param[in] old_st The original state
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  Return: the function return the final state: the proposed if the velocity rescaling is successfull, or the initial, if not.

  This version implies that the adiabatic representation is used
*/
  int nnucl = p.size();

  int st;

  if(new_st!=old_st){

    // According to Fabiano
    double a_ij = 0.0;
    double b_ij = 0.0;

    for(int k=0;k<nnucl;k++){

      double D = dc1_adi[k].get(old_st,new_st).real(); // derivative coupling w.r.t. nuclear DOF k

      a_ij += 0.5*(D*D / masses[k]); 
      b_ij += (D*p[k])/masses[k];
    }
    double det = b_ij*b_ij + 4.0*a_ij*(ham_adi.get(old_st,old_st).real() - ham_adi.get(new_st,new_st).real());

    // Calculate the scaling factor and new state
    double gamma_ij = 0.0;

    if(det<0.0){

      if(fabs(a_ij)>1e-100){  // only consider reversals, if the couplings are sizable
        if(do_reverse){     gamma_ij = b_ij / a_ij;}
         else{ gamma_ij = 0.0;  }
      }
      else{  gamma_ij = 0.0;    } // don't consider reversal, if the couplings are too small

      st = old_st; // # hop does not occur - frustrated hop

    }
    else{
      if(fabs(a_ij)>1e-100){ // only compute the rescaling factor and do the hop, if the couplings are sizable

        if(b_ij<0){ gamma_ij = 0.5*(b_ij + sqrt(det))/a_ij; }
        else{       gamma_ij = 0.5*(b_ij - sqrt(det))/a_ij; }

        st = new_st;
      } 
      else{   // otherwise - don't 
        gamma_ij = 0.0; 
        st = old_st;
      } 

    }

    //Rescale velocities and do the hop
    for(int k=0;k<nnucl;k++){ 
      double D = dc1_adi[k].get(old_st,new_st).real(); 
      p[k] = p[k] - gamma_ij * D; 
    }

  }
  else{ st = old_st; }

  return st;

} // rescale velocities adiabatic




void rescale_velocities_adiabatic(Nuclear* mol, Hamiltonian* ham, int& new_st,int& old_st, int do_reverse){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in,out] new_st The index of the new state: this is a truly new state if the attempted hop was successfull, or just can be 
             an old state, otherwise
  \param[in] old_st The index of the old state (from which we try to hop)
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  This verions implies that the adiabatic representation is used
*/

  int st;

  //cout<<"in rescale_velocities_adiabatic\n";

  if(new_st!=old_st){

    int final_st = old_st;  // no hopping by default

    // According to Fabiano
    double a_ij = 0.0;
    double b_ij = 0.0;

    for(int k=0;k<mol->nnucl;k++){

      double D = ham->D(old_st,new_st,k).real(); // derivative coupling w.r.t. nuclear DOF k
      cout<<"D= "<<D<<endl;

      a_ij += 0.5*(D*D / mol->mass[k]); 
      b_ij += (D*mol->p[k])/mol->mass[k];
    }
    double det = b_ij*b_ij + 4.0*a_ij*(ham->Hvib(old_st,old_st).real() - ham->Hvib(new_st,new_st).real());

    cout<<"det= "<<det<<" a_ij= "<<a_ij<<"  b_ij= "<<b_ij<<endl;

    // Calculate the scaling factor and new state
    double gamma_ij = 0.0;

    if(det<0.0){

      if(do_reverse){     gamma_ij = b_ij / a_ij;}
      else{ gamma_ij = 0.0;  }

      final_st = old_st; // # hop does not occur - frustrated hop

    }
    else{
      cout<<"gamma_ij= "<<gamma_ij<<endl;

      if(b_ij<0){ gamma_ij = 0.5*(b_ij + sqrt(det))/a_ij; }
      else{       gamma_ij = 0.5*(b_ij - sqrt(det))/a_ij; }
      final_st = new_st;
    }

    //Rescale velocities and do the hop
    for(int k=0;k<mol->nnucl;k++){ 
      double D = ham->D(old_st,new_st,k).real(); 
      cout<<"D= "<<D<<endl;
      mol->p[k] = mol->p[k] - gamma_ij * D; 
    }

    st = final_st;

  }
  else{ st = old_st; }

  new_st = st;


} // rescale velocities adiabatic

int rescale_velocities_adiabatic(Nuclear& mol, Hamiltonian& ham, int old_st, int do_reverse){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another - Python-friendly
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in] old_st The index of the old state (from which we try to hop)
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  The function returns the index of the new state: this is a truly new state if the attempted hop was successfull, or just can be 
  an old state, otherwise

  This verions implies that the adiabatic representation is used
*/


  int new_st = old_st;
  rescale_velocities_adiabatic(&mol, &ham, new_st, old_st, do_reverse);
  return new_st;

}




int rescale_velocities_diabatic(MATRIX& p, MATRIX& invM, CMATRIX* ham_adi, int new_st,int old_st){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another and do the
  momenta rescaling if the transition is allowed.

  \param[in,out] p [ndof x 1] nucler momenta. Can be updated (velocity rescaling or reversal)
  \param[in] invM [ndof x 1] matrix of inverse masses for all DOFs
  \param[in] ham_adi [nadi x nadi]  matrix of adiabatic Hamiltonian 
  \param[in] new_st The index of the proposed state
  \param[in] old_st The index of the old state (from which we try to hop)

  This versions implies that the diabatic representation is used

  In this case, derivative couplings are zero by definition of diabatic (position-independent) states
  so one can not rescale velocities along directions of derivative couplings, since there is no such directions
  We just scale velocities uniformly - based on energy conservation principle
  In principle, this rescaling procedure can be applied to surface hopping scheme in adiabatic basis, too

  Return: the function returns the final state: the proposed, if the velocity rescaling is successfull, or the initial, if not.

*/

  int st = old_st;

  double T_i = compute_kinetic_energy(p, invM); // initial kinetic energy

  double E_i = ham_adi->get(old_st,old_st).real();// initial potential energy
  double E_f = ham_adi->get(new_st,new_st).real();// final potential energy
  
  double T_f = T_i + E_i - E_f;             // predicted final kinetic energy

  if(T_f>0.0){  // hop is possible - accept it
 
    double scl = 0.0; // rescaling factor
     
    if(T_i>0.0){   scl = std::sqrt(T_f/T_i);  }
    else{ scl = 0.0; }

    // Rescale velocities
    p *= scl;

    // state stays the same
    st = new_st;

  }
  else{  // hop is not possible - reject it
    st = old_st;
  }

  return st;

} // rescale velocities diabatic



int rescale_velocities_diabatic(MATRIX& p, MATRIX& invM, CMATRIX& ham_adi, int new_st,int old_st){
/**
  \brief See the documentation of the 
  int rescale_velocities_diabatic(MATRIX& p, MATRIX& invM, CMATRIX* ham_adi, int new_st,int old_st)
  function
*/

  return rescale_velocities_diabatic(p, invM, &ham_adi, new_st, old_st);
}


int rescale_velocities_diabatic(MATRIX& p, MATRIX& invM, nHamiltonian* ham, int new_st,int old_st){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another and do the
  momenta rescaling if the transition is allowed.

  \param[in,out] p [ndof x 1] nucler momenta. Can be updated (velocity rescaling or reversal)
  \param[in] invM [ndof x 1] matrix of inverse masses for all DOFs
  \param[in] ham Is a nHamiltonian object that takes care of all energy-related transformations
  \param[in] new_st The index of the proposed state
  \param[in] old_st The index of the old state (from which we try to hop)

  This version implies that the diabatic representation is used

  In this case, derivative couplings are zero by definition of diabatic (position-independent) states
  so one can not rescale velocities along directions of derivative couplings, since there is no such directions
  We just scale velocities uniformly - based on energy conservation principle
  In principle, this rescaling procedure can be applied to surface hopping scheme in adiabatic basis, too

  Return: the function returns the final state: the proposed, if the velocity rescaling is successfull, or the initial, if not.

*/

  return rescale_velocities_diabatic(p, invM, ham->ham_adi, new_st, old_st);
}

int rescale_velocities_diabatic(MATRIX& p, MATRIX& invM, nHamiltonian& ham, int new_st,int old_st){
/**
  \brief See the documentation of the 
  int rescale_velocities_diabatic(MATRIX& p, MATRIX& invM, nHamiltonian* ham, int new_st,int old_st)
  function
*/

  return rescale_velocities_diabatic(p, invM, ham.ham_adi, new_st, old_st);

}



void rescale_velocities_diabatic(Nuclear* mol, Hamiltonian* ham, int& new_st,int& old_st){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in,out] new_st The index of the new state: this is a truly new state if the attempted hop was successfull, or just can be 
             an old state, otherwise
  \param[in] old_st The index of the old state (from which we try to hop)

  This verions implies that the diabatic representation is used

  In this case, derivative couplings are zero by definition of diabatic (position-independent) states
  so one can not rescale velocities along directions of derivative couplings, since there is no such directions
  We just scale velocities uniformly - based on energy conservation principle
  In principle, this rescaling procedure can be applied to surface hopping scheme in adiabatic basis, too
*/


  double T_i = compute_kinetic_energy(mol); // initial kinetic energy
  double E_i = ham->Hvib(old_st,old_st).real();// initial potential energy
  double E_f = ham->Hvib(new_st,new_st).real();// final potential energy
  
  double T_f = T_i + E_i - E_f;             // predicted final kinetic energy

  if(T_f>0.0){  // hop is possible - accept it
 
    double scl = 0.0; // rescaling factor
     
    if(T_i>0.0){   scl = std::sqrt(T_f/T_i);  }
    else{ scl = 0.0; }

    // Rescale velocities
    for(int i=0;i<mol->nnucl;i++){  mol->p[i] *= scl; }

    // state stays the same

  }
  else{  // hop is not possible - reject it

    new_st = old_st;

  }

} // rescale velocities diabatic

int rescale_velocities_diabatic(Nuclear& mol, Hamiltonian& ham, int old_st){
/**
  \brief Determine whether we need to do velocity rescaling/reversal when going from one state to another - Python-friendly
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in] old_st The index of the old state (from which we try to hop)

  The function returns the index of the new state: this is a truly new state if the attempted hop was successfull, or just can be 
  an old state, otherwise

  This verions implies that the diabatic representation is used

  In this case, derivative couplings are zero by definition of diabatic (position-independent) states
  so one can not rescale velocities along directions of derivative couplings, since there is no such directions
  We just scale velocities uniformly - based on energy conservation principle
  In principle, this rescaling procedure can be applied to surface hopping scheme in adiabatic basis, too
*/


  int new_st = old_st;
  rescale_velocities_diabatic(&mol, &ham, new_st, old_st);
  return new_st;

}





void rescale_velocities_adiabatic(int ntraj, vector<Nuclear*>& mol, vector<Hamiltonian*>& ham,
      vector<int>& new_st, vector<int>& old_st, int do_reverse){

  // Calculate auxiliary variables to determine the case
  // Here i - old state, j - new state
  vector<int> st;

  int status = 1; // new_st == old_st
  for(int traj=0; traj<ntraj; traj++){
    if(new_st[traj]!=old_st[traj]){ status=0; }
  }

  if(status==0){  // different states

    vector<int> final_st; final_st = old_st;  // no hopping by default

    // According to Fabiano
    double a_ij = 0.0;
    double b_ij = 0.0;
    double E_old = 0.0;
    double E_new = 0.0;
    
    for(int traj=0;traj<ntraj;traj++){
      for(int k=0;k<mol[traj]->nnucl;k++){

        double D = ham[traj]->D(old_st[traj],new_st[traj],k).real(); // derivative coupling w.r.t. nuclear DOF k

        a_ij += 0.5*(D*D / mol[traj]->mass[k]); 
        b_ij += (D*mol[traj]->p[k])/mol[traj]->mass[k];

      }// for k

      E_old += ham[traj]->Hvib(old_st[traj],old_st[traj]).real();
      E_new += ham[traj]->Hvib(new_st[traj],new_st[traj]).real();

    }


    double det = b_ij*b_ij + 4.0*a_ij*(E_old - E_new);


    // Calculate the scaling factor and new state
    double gamma_ij = 0.0;

    if(det<0.0){

      if(do_reverse){     gamma_ij = b_ij / a_ij;}
      else{ gamma_ij = 0.0;  }

      final_st = old_st; // # multi-particle hop does not occur - frustrated hop

    }
    else{
      if(b_ij<0){ gamma_ij = 0.5*(b_ij + sqrt(det))/a_ij; }
      else{       gamma_ij = 0.5*(b_ij - sqrt(det))/a_ij; }
      final_st = new_st;
    }

    //Rescale velocities and do the hop
    for(int traj=0; traj<ntraj; traj++){
      for(int k=0;k<mol[traj]->nnucl;k++){ 

        double D = ham[traj]->D(old_st[traj],new_st[traj],k).real(); 
        mol[traj]->p[k] = mol[traj]->p[k] - gamma_ij * D; 

      }// for k
    }// for traj

    st = final_st;

  }
  else{ st = old_st; }

  new_st = st;


} // rescale velocities adiabatic



void rescale_velocities_adiabatic(int ntraj, vector<Nuclear>& mol, vector<Hamiltonian>& ham,
      vector<int>& new_st, vector<int>& old_st, int do_reverse){

  // Calculate auxiliary variables to determine the case
  // Here i - old state, j - new state
  vector<int> st;

  int status = 1; // new_st == old_st
  for(int traj=0; traj<ntraj; traj++){
    if(new_st[traj]!=old_st[traj]){ status=0; }
  }

  if(status==0){  // different states

    vector<int> final_st; final_st = old_st;  // no hopping by default

    // According to Fabiano
    double a_ij = 0.0;
    double b_ij = 0.0;
    double E_old = 0.0;
    double E_new = 0.0;
    
    for(int traj=0;traj<ntraj;traj++){
      for(int k=0;k<mol[traj].nnucl;k++){

        double D = ham[traj].D(old_st[traj],new_st[traj],k).real(); // derivative coupling w.r.t. nuclear DOF k

        a_ij += 0.5*(D*D / mol[traj].mass[k]); 
        b_ij += (D*mol[traj].p[k])/mol[traj].mass[k];

      }// for k

      E_old += ham[traj].Hvib(old_st[traj],old_st[traj]).real();
      E_new += ham[traj].Hvib(new_st[traj],new_st[traj]).real();

    }


    double det = b_ij*b_ij + 4.0*a_ij*(E_old - E_new);


    // Calculate the scaling factor and new state
    double gamma_ij = 0.0;

    if(det<0.0){

      if(do_reverse){     gamma_ij = b_ij / a_ij;}
      else{ gamma_ij = 0.0;  }

      final_st = old_st; // # multi-particle hop does not occur - frustrated hop

    }
    else{
      if(b_ij<0){ gamma_ij = 0.5*(b_ij + sqrt(det))/a_ij; }
      else{       gamma_ij = 0.5*(b_ij - sqrt(det))/a_ij; }
      final_st = new_st;
    }

    //Rescale velocities and do the hop
    for(int traj=0; traj<ntraj; traj++){
      for(int k=0;k<mol[traj].nnucl;k++){ 

        double D = ham[traj].D(old_st[traj],new_st[traj],k).real(); 
        mol[traj].p[k] = mol[traj].p[k] - gamma_ij * D; 

      }// for k
    }// for traj

    st = final_st;

  }
  else{ st = old_st; }

  new_st = st;


} // rescale velocities adiabatic



/*
int rescale_velocities_adiabatic(Nuclear& mol, Hamiltonian& ham, int old_st, int do_reverse){

  int new_st = old_st;
  rescale_velocities_adiabatic(&mol, &ham, new_st, old_st, do_reverse);
  return new_st;

}
*/






}// namespace libdyn
}// liblibra

