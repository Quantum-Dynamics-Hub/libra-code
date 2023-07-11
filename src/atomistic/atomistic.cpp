/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Hamiltonian_Atomistic.cpp
  \brief The file implements the high-level computations in the atomistic Hamiltonians - allowing to
  mix different resolutions (e.g. classical + quantum) and approaches (different force fields or quantum methods), fragmentation.
*/


#if defined(USING_PCH)
#include "../pch.h"
#else
#include <complex>
#include <cmath>
#endif 

#include "../math_linalg/liblinalg.h"
#include "../math_meigen/libmeigen.h"
#include "../calculators/libcalculators.h"
#include "atomistic.h"


/// liblibra namespace
namespace liblibra{

namespace libatomistic{


using namespace libcalculators;
using namespace libchemobjects;
using namespace libchemobjects::libchemsys;
using namespace libhamiltonian_mm;
using namespace libhamiltonian_qm;
using namespace liblinalg;
using namespace libmeigen;
using std::complex;
using std::sin;
using std::cos;
using std::exp;
using std::sqrt;


Hamiltonian_Atomistic::Hamiltonian_Atomistic(int _nelec, int _nnucl){
/**
  \param[in] _nelec The number of electronic basis states 
  \param[in] _nnucl The number of nuclear DOF

  This Hamiltonian_Atomistic constructor allocates internal memory for ham_dia, ham_adi, d1ham_dia, d1ham_adi, 
  and d2ham_dia and populates them with zeroes. The default representation is selected, status flags are 
  set to zero (Hamiltonian is not up to date) and the hamiltonian types are initialized to 0 (meaning actual Hamiltonians
  are not initialized)
*/

  _syst = NULL;
  
  int i;

  nelec = _nelec;
  nnucl = _nnucl;

  ham_dia = new MATRIX(nelec,nelec); *ham_dia = 0.0;
  ham_adi = new MATRIX(nelec,nelec); *ham_adi = 0.0;

  d1ham_dia = vector<MATRIX*>(nnucl);
  d1ham_adi = vector<MATRIX*>(nnucl);

  for(i=0;i<nnucl;i++){  
    d1ham_dia[i] = new MATRIX(nelec,nelec); *d1ham_dia[i] = 0.0; 
    d1ham_adi[i] = new MATRIX(nelec,nelec); *d1ham_adi[i] = 0.0; 
  }

  d2ham_dia = vector<MATRIX*>(nnucl*nnucl);
  for(i=0;i<nnucl*nnucl;i++){  d2ham_dia[i] = new MATRIX(nelec,nelec); *d2ham_dia[i] = 0.0;  }


  rep = 1; // default representation is adiabatic  
  status_dia = 0;
  status_adi = 0;


  // Setup Hamiltonian types:
  ham_types = vector<int>(5,0);  // which Hamiltonians are initialized ham_types[0] = MM, ham_types[1] = QM

}


void Hamiltonian_Atomistic::set_Hamiltonian_type(std::string ham_type){ // libchemobjects::libchemsys::System& syst){
/**
  \param[in] ham_type The type of Hamiltonian to initialize. Options: "MM" and "QM"

  Initialize the selected Hamiltonian type, if not done yet. The adiabatic representation will be chosen.
*/

  if(ham_type=="MM"){
    rep = 1;  // diabatic is same as adiabatic

    if(ham_types[0]==0){    
      mm_ham = new listHamiltonian_MM();
      ham_types[0] = 1;
    }

  }
  else if(ham_type=="QM"){
    rep = 1;  // adiabatic

    if(ham_types[1]==0){    
      qm_ham = new listHamiltonian_QM();
      ham_types[1] = 1;
    }

  }

  else{
    cout<<"Error: Unrecognized Hamiltonian type. Supported types are: \"MM\" and \"QM\" \n";
    cout<<"Exiting...\n"; exit(0);
  }

}


void Hamiltonian_Atomistic::show_interactions_statistics(){
  mm_ham->show_interactions_statistics();
}

void Hamiltonian_Atomistic::show_interactions(){
  mm_ham->show_interactions();
}


void Hamiltonian_Atomistic::set_atom_types(System& syst, vector<int>& lst, ForceField& ff){
  mm_ham->set_atom_types(syst, lst, ff);
}

void Hamiltonian_Atomistic::set_fragment_types(System& syst, vector<int>& lst, ForceField& ff){
  mm_ham->set_fragment_types(syst, lst, ff);  
}

bool Hamiltonian_Atomistic::is_active(Atom& a1, Atom& a2){
  return mm_ham->is_active(a1,a2);
}

bool Hamiltonian_Atomistic::is_active(Atom& a1, Atom& a2, Atom& a3){
  return mm_ham->is_active(a1,a2,a3);
}

bool Hamiltonian_Atomistic::is_active(Atom& a1, Atom& a2, Atom& a3, Atom& a4){
  return mm_ham->is_active(a1,a2,a3,a4);
}


void Hamiltonian_Atomistic::set_atom_interactions_for_atoms
(System& syst,string int_type,vector<Atom>& top_elt,vector<int>& lst1,vector<int>& lst2,ForceField& ff,int verb){

  mm_ham->set_atom_interactions_for_atoms(syst,int_type,top_elt,lst1,lst2,ff,verb);

}

void Hamiltonian_Atomistic::set_group_interactions_for_atoms
(System& syst,string int_type,vector<Group>& top_elt,vector<int>& lst1,vector<int>& lst2,ForceField& ff){

  mm_ham->set_group_interactions_for_atoms(syst,int_type,top_elt,lst1,lst2,ff);
}

void Hamiltonian_Atomistic::set_interactions_for_atoms
(System& syst, boost::python::list lst1, boost::python::list lst2, ForceField& ff, int verb, int assign_rings){

  mm_ham->set_interactions_for_atoms(syst,lst1,lst2,ff,verb,assign_rings);
}

void Hamiltonian_Atomistic::set_interactions_for_fragments
(System& syst, boost::python::list lst1, boost::python::list lst2, ForceField& ff){

  mm_ham->set_interactions_for_fragments(syst,lst1,lst2,ff);

}

void Hamiltonian_Atomistic::apply_pbc_to_interactions(System& syst, int int_type,int nx,int ny,int nz){

  mm_ham->apply_pbc_to_interactions(syst, int_type, nx, ny, nz);
}

void Hamiltonian_Atomistic::set_respa_types(std::string inter_type,std::string respa_type){

  mm_ham->set_respa_types(inter_type, respa_type);
}




Hamiltonian_Atomistic::~Hamiltonian_Atomistic(){
/**
  Destructor: Free memory occupied by ham_dia, ham_adi, d1ham_dia, d1ham_adi, d2ham_dia
*/

  int i;

  delete ham_dia;
  delete ham_adi;

  for(i=0;i<nnucl;i++){  
    delete d1ham_dia[i];
    delete d1ham_adi[i];
  }
  for(i=0;i<nnucl*nnucl;i++){  delete d2ham_dia[i]; }

  if(d1ham_dia.size()>0) { d1ham_dia.clear(); }
  if(d2ham_dia.size()>0) { d2ham_dia.clear(); }
  if(d1ham_adi.size()>0) { d1ham_adi.clear(); }


}


void Hamiltonian_Atomistic::set_q(vector<double>& q_){
/**
  Update Hamiltonian coordinates (all are the real-valued scalars)

  \param[in] q The vector of real-valued coordinates to be used for Hamiltonian calculations.
  Note: this also sets the status_dia and status_adi variables to zero, impliying the Hamiltonian is not 
  up to date - which is what we want: since the coordinates are changed, the Hamiltonian must be recomputed
  From the prafmatic point of view, if you call this function - expect slower performance.

  Keep in mind that this will also change the coordinates of the external system object bound to this Hamiltonian
*/


  if(q_.size()!=nnucl){
    cout<<"Error in Hamiltonian_Atomistic::set_q - the size of input array ("<<q_.size()<<") does not match";
    cout<<"the number of nuclear degrees of freedom ("<<nnucl<<")\nExiting...\n"; exit(0);
  }

  q = q_;
  status_dia = 0;
  status_adi = 0;

  _syst->set_atomic_q(q);
  //_syst->set_fragment_q(q);

}

void Hamiltonian_Atomistic::set_q(boost::python::list q_){
/**
  Update Hamiltonian coordinates (all are the real-valued scalars) - Python-friendly

  \param[in] q The list of real-valued coordinates to be used for Hamiltonian calculations.
  Note: this also sets the status_dia and status_adi variables to zero, impliying the Hamiltonian is not 
  up to date - which is what we want: since the coordinates are changed, the Hamiltonian must be recomputed
  From the prafmatic point of view, if you call this function - expect slower performance.

  Keep in mind that this will also change the coordinates of the external system object bound to this Hamiltonian
*/

  int sz = len(q_);
  vector<double> q__(sz,0.0);
  for(int i=0;i<sz;i++){ q__[i] = extract<double>(q_[i]); }

  set_q(q__);

}

void Hamiltonian_Atomistic::set_v(vector<double>& v_){
/**
  Update Hamiltonian velocities (all are real-valued scalars)

  \param[in] v The vector of real-valued velocities to be used for Hamiltonian calculations.

  The velocities are only needed for vibronic Hamiltonian (adiabatic representation) calculations. 
  Otherwise, they are not used.
  Only status_adi is set to 0, so only adiabatic Hamiltonian is recomputed.
  For future: in fact, we only need to update the vibronic Hamiltonian, so we still may save a lot, when adiabatic
  calculations imply electronic structure calculations

  Keep in mind that this will also change the velocities/momenta of the external system object bound to this Hamiltonian
*/


  if(v_.size()!=nnucl){
    cout<<"Error in Hamiltonian_Atomistic::set_v - the size of input array ("<<v_.size()<<") does not match";
    cout<<"the number of nuclear degrees of freedom ("<<nnucl<<")\nExiting...\n"; exit(0);
  }

  v = v_;
  status_adi = 0;  // only affects adiabatic computations

  _syst->set_atomic_v(v);
  //_syst->set_fragment_v(v);
  

}

void Hamiltonian_Atomistic::set_v(boost::python::list v_){
/**
  Update Hamiltonian velocities (all are real-valued scalars) - Python-friendly

  \param[in] v The Python list of real-valued velocities to be used for Hamiltonian calculations.

  The velocities are only needed for vibronic Hamiltonian (adiabatic representation) calculations. 
  Otherwise, they are not used.
  Only status_adi is set to 0, so only adiabatic Hamiltonian is recomputed.
  For future: in fact, we only need to update the vibronic Hamiltonian, so we still may save a lot, when adiabatic
  calculations imply electronic structure calculations

  Keep in mind that this will also change the velocities/momenta of the external system object bound to this Hamiltonian
*/


  int sz = len(v_);
  vector<double> v__(sz,0.0);
  for(int i=0;i<sz;i++){ v__[i] = extract<double>(v_[i]); }

  set_v(v__);

}




/*
void Hamiltonian_Atomistic::set_rep(int rep_){
  rep = rep_;
}

//  ham->set_q(mol->q);

void Hamiltonian_Atomistic::set_params(vector<double>& params_){

  
  params = vector<double>(params_.size(), 0.0);

  // Now copy input params:
  for(int i=0;i<params_.size(); i++){
    params[i] = params_[i];
  }

  // Since the parameters have changed - we need to recompute everything
  status_dia = 0;
  status_adi = 0;

}

void Hamiltonian_Atomistic::set_params(boost::python::list params_){

  int sz = boost::python::len(params_);
  vector<double> tmp_params(sz, 0.0);

  // Now copy input params:
  for(int i=0;i<sz; i++){
    tmp_params[i] = boost::python::extract<double>(params_[i]);
  }

  set_params(tmp_params);

}



void Hamiltonian_Atomistic::set_v(vector<double>& v_){
  v = v_;
  status_adi = 0;  // only affects adiabatic computations
}

void Hamiltonian_Atomistic::set_v(boost::python::list v_){

  int sz = boost::python::len(v_);
  vector<double> tmp_v(sz,0.0);

  for(int i=0;i<sz; i++){
    tmp_v[i] = boost::python::extract<double>(v_[i]);
  }

  set_v(tmp_v);

}


void Hamiltonian_Atomistic::compute(){

  if(rep==0){  compute_diabatic();  }
  else if(rep==1){  compute_adiabatic(); }

}

*/
void Hamiltonian_Atomistic::compute_diabatic(){
/**
  This function computes diabatic PESs and derivative couplings

  For MM - actual computations are performed here and then are used in compute_adiabatic (if needed)
  For QM - we would call the SCF solver for orbital picture description, but we actually use the 
  true excitonic description, so the nothing is actually done here for such type of calculations
  
  In other words, we use method_option = 1 
  method_option = 0 is here, but is never reached, since this parameter is 
  set inside of the code - you'd need to recompile the code to get to tht part). That part is left in the 
  code for future developments.
*/

  int i;

  //cout<<"in Hamiltonian_Atomistic::compute_diabatic()\n";

  if(status_dia == 0){ // only compute this if the result is not up to date

    // Zero forces   
    _syst->zero_forces_and_torques();
    *ham_dia = 0.0;
    for(int n=0;n<nnucl;n++){   *d1ham_dia[n] = 0.0;   }
     
  
    if(ham_types[0]==1){

      // Do actual computations
      int sz = mm_ham->interactions.size();
      double res = 0.0;
      int tmp;
      for(int i=0;i<sz;i++){
        res += mm_ham->interactions[i].calculate(tmp);
        //cout<<"interactions #"<<i<<", energy = "<<mm_ham->interactions[i].calculate(tmp)<<endl;
      }

      for(int st=0;st<nelec;st++){
        // Energies
        ham_dia->M[st*nelec+st] += res;
      
        // First derivatives     
        // But take only atomistic forces at this time
        for(int i=0;i<_syst->Number_of_atoms;i++){         
          d1ham_dia[3*i  ]->M[st*nelec+st] -= _syst->Atoms[i].Atom_RB.rb_force.x;
          d1ham_dia[3*i+1]->M[st*nelec+st] -= _syst->Atoms[i].Atom_RB.rb_force.y;
          d1ham_dia[3*i+2]->M[st*nelec+st] -= _syst->Atoms[i].Atom_RB.rb_force.z;
        }

        // Second derivatives
        for(i=0;i<nnucl*nnucl;i++){    d2ham_dia[i]->M[st*nelec+st] = 0.0;      }

      } // for i


    }// MM

    if(ham_types[1]==1){

      int method_option = 1;  // 
      
      if(method_option==0){ // Orbital picture

        _syst->zero_forces_and_torques();
        qm_ham->compute_scf(*_syst);

        // Energies
        MATRIX* tmp; tmp = new MATRIX(nelec,nelec);
        *tmp = 0.0;

        if(nelec < qm_ham->el->Fao_alp->n_cols){
      
          vector<int> subset(nelec);
          for(int i=0;i<nelec;i++){ subset[i] = i; }  
          pop_submatrix(qm_ham->el->Fao_alp, tmp, subset);

        }
        else{

          vector<int> subset(qm_ham->el->Fao_alp->n_cols);
          for(int i=0;i<qm_ham->el->Fao_alp->n_cols;i++){ subset[i] = i; }  
          push_submatrix(tmp, qm_ham->el->Fao_alp, subset);

        }

        *ham_dia += *tmp;

        delete tmp;


        if(0){
          cout<<"internal Fao = \n"<<*qm_ham->el->Fao_alp<<endl;
          cout<<"diabatic Ham = \n"<<*ham_dia<<endl;
        }
      }// method == 0 - orbital picture

      else if(method_option==1){ // Excitonic case

        // don't do anything here - for atomistic system we only use adiabatic basis anyways
      }


    }// QM
    //==========================================================

    // class - specific computations: MM, semi-empirical, HF, etc   

    //==========================================================

    // Set status flag 
    status_dia = 1;

  }//   status_dia == 0

  //cout<<"finishing Hamiltonian_Atomistic::compute_diabatic()\n";

}



void Hamiltonian_Atomistic::compute_adiabatic(){
/**
  This function computes adiabatic PESs and derivative couplings

  For MM - the actual computations are already done in the compute_diabatic() (called from here)
  For QM - here is gonna be SCF solver, to get adiabatic states
  
  Note: at this point, no diabatic calculations are performed for QM in compute_diabatic - for QM 
  we use only adiabatic so far. Moreover, we use the total energies/many-electron determinants rather
  than 1-electron orbitals and their energies - so see only the section which comes with 
  method_option = 1 (method_option = 0 is there, but is never reached, since this parameter is 
  set inside of the code - you'd need to recompile the code to get to tht part). That part is left in the 
  code for future developments.
*/

  //cout<<"in Hamiltonian_Atomistic::compute_adiabatic()\n";

  Timer tim1;
  compute_diabatic();

  if(status_adi == 0){
    // Zero forces   
    _syst->zero_forces_and_torques();

    *ham_adi = 0.0;
    for(int n=0;n<nnucl;n++){    *d1ham_adi[n] = 0.0;    }


    if(ham_types[0]==1){  // MM hamiltonian

      *ham_adi += *ham_dia;  
      
      // Just use diabatic forces
      for(int n=0;n<nnucl;n++){       *d1ham_adi[n] += *d1ham_dia[n];      }// for n


    }// ham_type[0]==1

    if(ham_types[1]==1){

      int method_option = 1; // 0 - simple case: energies are the orbital energies
                             // 1 - excitonic case: energies are the total energies for different density matrices

      if(method_option==0){
        /// method_option == 0
        /// For this case ham_dia and d1ham_dia already contain (possible) contribution from MM part, so
        /// we just substitute the final ham_adi and d1ham_adi with the obtained result, no need to keep track of the 
        /// above adiabatic MM contributions (which is seemingly disregarded, but is actually included via diabatic part)
 
        MATRIX* S; S = new MATRIX(nelec, nelec);  *S = 0.0;
        MATRIX* C; C = new MATRIX(nelec, nelec);  *C = 0.0;

        if(nelec < qm_ham->el->Sao->n_cols){      
          vector<int> subset(nelec); 
          for(int i=0;i<nelec;i++){ subset[i] = i; }  
          pop_submatrix(qm_ham->el->Sao, S, subset);
        }
        else{
          vector<int> subset(qm_ham->el->Sao->n_cols); 
          for(int i=0;i<qm_ham->el->Sao->n_cols;i++){ subset[i] = i; }  
          push_submatrix(S, qm_ham->el->Sao, subset);

          for(int i=qm_ham->el->Sao->n_cols;i<nelec;i++){
            S->set(i,i,1.0);
          }
        }
      
        // Debug:
        //cout<<"Sao = \n"<<*qm_ham->el->Sao<<endl;
        //cout<<"S= \n"<<*S<<endl;
        //cout<<"internal MO = \n"<<*qm_ham->el->C_alp<<endl;

        // Transformation to adiabatic basis
        solve_eigen(ham_dia, S, ham_adi, C, 0);  // H_dia * C = S * C * H_adi

        // Debug:
        //cout<<"new C =\n"<<*C<<endl;
        //cout<<"internal E_alp = \n"<<*qm_ham->el->E_alp<<endl;
        //cout<<"adiabatic Ham = \n"<<*ham_adi<<endl;

        // Now compute the derivative couplings (off-diagonal, multiplied by energy difference) and adiabatic gradients (diagonal)
        for(int n=0;n<nnucl;n++){   *d1ham_adi[n] = (*C).T() * (*d1ham_dia[n]) * (*C);    }// for n

        delete S;
        delete C;

      }// method_option == 0

      else if(method_option==1){ 
        /// method_option==1: This is true case - excitonic one

        int x_period = 0;  int y_period = 0;  int z_period = 0; 
        VECTOR t1, t2, t3;

        tim1.start();

        int i = 0;
        _syst->zero_forces_and_torques();

        /// Compute the ground state energy first
        double E_i = qm_ham->energy_and_forces(*_syst); // ground state energy
        int Norb = qm_ham->el->Norb;

        // Additive contribution - so to potentially include MM part
        ham_adi->M[i*nelec+i] += E_i;

        cout<<"Ground state calculations, time = "<<tim1.stop()<<endl;


        tim1.start();
      
        MATRIX* Dmo_a_x; Dmo_a_x = new MATRIX(Norb,Norb);
        MATRIX* Dmo_a_y; Dmo_a_y = new MATRIX(Norb,Norb);
        MATRIX* Dmo_a_z; Dmo_a_z = new MATRIX(Norb,Norb);
        MATRIX* Dmo_b_x; Dmo_b_x = new MATRIX(Norb,Norb);
        MATRIX* Dmo_b_y; Dmo_b_y = new MATRIX(Norb,Norb);
        MATRIX* Dmo_b_z; Dmo_b_z = new MATRIX(Norb,Norb);


        for(int n=0;n<_syst->Number_of_atoms;n++){         

          // Additive contribution - so to potentially include MM part      
          d1ham_adi[3*n  ]->M[i*nelec+i] -= _syst->Atoms[n].Atom_RB.rb_force.x;
          d1ham_adi[3*n+1]->M[i*nelec+i] -= _syst->Atoms[n].Atom_RB.rb_force.y;
          d1ham_adi[3*n+2]->M[i*nelec+i] -= _syst->Atoms[n].Atom_RB.rb_force.z;
            
          //-------- Also, here we want to compute NACs ----------------

          // ----------- Derivative couplings are set once and for all -------------------
          derivative_couplings1(*qm_ham->el, *_syst, qm_ham->basis_ao, qm_ham->prms, qm_ham->modprms,
             qm_ham->atom_to_ao_map, qm_ham->ao_to_atom_map, *qm_ham->el->Hao, *qm_ham->el->Sao,
             qm_ham->el->Norb, n, x_period, y_period, z_period, t1, t2, t3,
             *Dmo_a_x, *Dmo_a_y, *Dmo_a_z, *Dmo_b_x, *Dmo_b_y, *Dmo_b_z
          );

//          cout<<"Derivative couplings1 for n = "<<n<<endl;
          //exit(0);
          if(qm_ham->basis_ex.size()<nelec){
            cout<<"Error: in Hamiltonian_Atomistic::compute_adiabatic()...\n";
            cout<<"The number of excitations ("<<qm_ham->basis_ex.size()<<") is less than the dimensionality of the Hamiltonian ("<<nelec<<")\n";
            cout<<"Use the ham.add_excitation(int,int,int,int) function to add the appropriate number of excitations\n";
            cout<<"Exiting now...\n";
            exit(0);
          }
          else if(qm_ham->basis_ex.size()>nelec){
            cout<<"Warning: in Hamiltonian_Atomistic::compute_adiabatic()...\n";
            cout<<"The number of excitations ("<<qm_ham->basis_ex.size()<<") is larger than the dimensionality of the Hamiltonian ("<<nelec<<")\n";
            cout<<"Only first "<<nelec<<" states will be handled. Other excitations are disregarded\n";
          }

 
//          cout<<"nelec = "<<nelec<<endl;

          for(int I=0;I<nelec;I++){    // over all excitons

            // [0] - implies we only work with single excitations!
            int I_fo = qm_ham->basis_ex[I].from_orbit[0];
            int I_fs = qm_ham->basis_ex[I].from_spin[0];
            int I_to = qm_ham->basis_ex[I].to_orbit[0];
            int I_ts = qm_ham->basis_ex[I].to_spin[0];


            for(int J=0;J<nelec;J++){  // over all excitons

              int J_fo = qm_ham->basis_ex[J].from_orbit[0];
              int J_fs = qm_ham->basis_ex[J].from_spin[0];
              int J_to = qm_ham->basis_ex[J].to_orbit[0];
              int J_ts = qm_ham->basis_ex[J].to_spin[0];



              if(I!=J){

                double dc_x = 0.0;
                double dc_y = 0.0;
                double dc_z = 0.0;


                //========== Here we need to work out the conditions very carefully =========
                // but for now we will just compare the final (only) spins               


                if(I_ts==J_ts){  // final spins must be the same 


                  if(I_ts==1){ // alpha channel

                    if(I_fo==J_fo && I_to != J_to){  // electronic coupling
                    
                      int glob_I_to_indx = I_to + qm_ham->el->Nocc_alp; // index of the obital
                      int glob_J_to_indx = J_to + qm_ham->el->Nocc_alp; // index of the obital

                      dc_x = Dmo_a_x->get(glob_I_to_indx,glob_J_to_indx);
                      dc_y = Dmo_a_y->get(glob_I_to_indx,glob_J_to_indx);
                      dc_z = Dmo_a_z->get(glob_I_to_indx,glob_J_to_indx);

                    }
                    if(I_to==J_to && I_fo != J_fo){  // hole coupling

                      int glob_I_fo_indx = I_fo + qm_ham->el->Nocc_alp; // index of the obital
                      int glob_J_fo_indx = J_fo + qm_ham->el->Nocc_alp; // index of the obital

                      dc_x = Dmo_a_x->get(glob_I_fo_indx,glob_J_fo_indx);
                      dc_y = Dmo_a_y->get(glob_I_fo_indx,glob_J_fo_indx);
                      dc_z = Dmo_a_z->get(glob_I_fo_indx,glob_J_fo_indx);
        
                    }
                    else{ ;; }  // in this case couplings are zero

                  }// alpha


                  else if(I_ts==-1){ // beta channel

                    if(I_fo==J_fo && I_to != J_to){  // electronic coupling

                      int glob_I_to_indx = I_to + qm_ham->el->Nocc_alp; // index of the obital
                      int glob_J_to_indx = J_to + qm_ham->el->Nocc_alp; // index of the obital

                      dc_x = Dmo_b_x->get(glob_I_to_indx,glob_J_to_indx);
                      dc_y = Dmo_b_y->get(glob_I_to_indx,glob_J_to_indx);
                      dc_z = Dmo_b_z->get(glob_I_to_indx,glob_J_to_indx);

                    }
                    if(I_to==J_to && I_fo != J_fo){  // hole coupling

                      int glob_I_fo_indx = I_fo + qm_ham->el->Nocc_alp; // index of the obital
                      int glob_J_fo_indx = J_fo + qm_ham->el->Nocc_alp; // index of the obital

                      dc_x = Dmo_b_x->get(glob_I_fo_indx,glob_J_fo_indx);
                      dc_y = Dmo_b_y->get(glob_I_fo_indx,glob_J_fo_indx);
                      dc_z = Dmo_b_z->get(glob_I_fo_indx,glob_J_fo_indx);
        
                    }
                    else{ ;; }  // in this case couplings are zero


                  }// beta



                }// I_ts == J_ts
                else{ ;; } // in this case couplings are zero


//                cout<<"I= "<<I<<" J= "<<J<<endl;
                
                d1ham_adi[3*n  ]->set(I,J, dc_x);
                d1ham_adi[3*n+1]->set(I,J, dc_y);
                d1ham_adi[3*n+2]->set(I,J, dc_z);

              }// I!=J

            }// for J
          }// for I


//          cout<<"Excitations are done:\n";

        }// for n

        delete Dmo_a_x;  delete Dmo_a_y;  delete Dmo_a_z;
        delete Dmo_b_x;  delete Dmo_b_y;  delete Dmo_b_z;


        cout<<"Derivative couplings computations for all atoms, time = "<<tim1.stop()<<endl;


        /// Then repeat the ground-state calculations with different MO occupation schemes (for excitations)
        //================ Now excitations ===========================
        vector< pair<int,double> > occ_alp_grnd(Norb,pair<int,double>(0,0.0)); occ_alp_grnd = qm_ham->el->occ_alp;
        vector< pair<int,double> > occ_bet_grnd(Norb,pair<int,double>(0,0.0)); occ_bet_grnd = qm_ham->el->occ_bet;

        tim1.start();

        for(int i=1;i<qm_ham->basis_ex.size();i++){  // over all excitons in this basis
        
          
          cout<<"excitation "<<i<<"  \n";
          qm_ham->el->occ_alp = occ_alp_grnd;
          qm_ham->el->occ_bet = occ_bet_grnd;

          libcalculators::excite(Norb, qm_ham->basis_ex[i], qm_ham->el->Nocc_alp, qm_ham->el->occ_alp, 
                                            qm_ham->el->Nocc_bet, qm_ham->el->occ_bet); // ground state excitation

          // And recompute density matrix
          compute_density_matrix(qm_ham->el->occ_alp, qm_ham->el->C_alp, qm_ham->el->P_alp);
          compute_density_matrix(qm_ham->el->occ_bet, qm_ham->el->C_bet, qm_ham->el->P_bet);

          // Also update the total density matrix:
          *qm_ham->el->P = *qm_ham->el->P_alp + *qm_ham->el->P_bet;

          // Update Fock
          Hamiltonian_Fock(qm_ham->el, *_syst, qm_ham->basis_ao, qm_ham->prms, qm_ham->modprms, qm_ham->atom_to_ao_map, qm_ham->ao_to_atom_map);


          // Finally, compute energy
          E_i = (energy_elec(qm_ham->el->P_alp, qm_ham->el->Hao, qm_ham->el->Fao_alp) +
                 energy_elec(qm_ham->el->P_bet, qm_ham->el->Hao, qm_ham->el->Fao_bet) 
                );

          cout<<"E_i = "<<E_i<<endl;
 

          // Forces and nuclear contributions:

          for(int n=0;n<_syst->Number_of_atoms;n++){

            _syst->Atoms[n].Atom_RB.rb_force = 
            force(*qm_ham->el, *_syst, qm_ham->basis_ao, qm_ham->prms, qm_ham->modprms,
                  qm_ham->atom_to_ao_map, qm_ham->ao_to_atom_map, *qm_ham->el->Hao, *qm_ham->el->Sao,
                  qm_ham->el->Norb, n, x_period, y_period, z_period, t1, t2, t3);

          }// for n


          // - nuclear-nuclear repulsion
          vector<double> Zeff;
          vector<VECTOR> G;
          vector<VECTOR> R;
          for(int n=0;n<_syst->Number_of_atoms;n++){
            R.push_back(_syst->Atoms[n].Atom_RB.rb_cm);
            Zeff.push_back(qm_ham->modprms.PT[_syst->Atoms[n].Atom_element].Zeff);
          }// for n

          double Enucl = energy_nucl(R, Zeff, G);
          for(int n=0;n<_syst->Number_of_atoms;n++){   _syst->Atoms[n].Atom_RB.rb_force -= G[n];   }

          E_i += Enucl;

                    

          // Additive contribution - so to potentially include MM part
          ham_adi->M[i*nelec+i] += E_i;


          for(int n=0;n<_syst->Number_of_atoms;n++){         

            // Additive contribution - so to potentially include MM part
            d1ham_adi[3*n  ]->M[i*nelec+i] -= _syst->Atoms[n].Atom_RB.rb_force.x;
            d1ham_adi[3*n+1]->M[i*nelec+i] -= _syst->Atoms[n].Atom_RB.rb_force.y;
            d1ham_adi[3*n+2]->M[i*nelec+i] -= _syst->Atoms[n].Atom_RB.rb_force.z;

          }
                  

        }// for i

        
        cout<<"Excited states energies and forces, time = "<<tim1.stop()<<endl;


        tim1.start();

        /// When excitations are computed, reset the density matrix to current electronic state and recompute
        /// electronic structure properties: Fock matrix, forces, derivative couplings, etc
        //=================== Reset current variables to be those of the ground state ===============
        i = 0;
        cout<<"excitation "<<i<<"  \n";
        qm_ham->el->occ_alp = occ_alp_grnd;
        qm_ham->el->occ_bet = occ_bet_grnd;

        libcalculators::excite(Norb, qm_ham->basis_ex[i], qm_ham->el->Nocc_alp, qm_ham->el->occ_alp, 
                                          qm_ham->el->Nocc_bet, qm_ham->el->occ_bet); // ground state excitation

        // And recompute density matrix
        compute_density_matrix(qm_ham->el->occ_alp, qm_ham->el->C_alp, qm_ham->el->P_alp);
        compute_density_matrix(qm_ham->el->occ_bet, qm_ham->el->C_bet, qm_ham->el->P_bet);

        // Also update the total density matrix:
        *qm_ham->el->P = *qm_ham->el->P_alp + *qm_ham->el->P_bet;

        // Update Fock
        Hamiltonian_Fock(qm_ham->el, *_syst, qm_ham->basis_ao, qm_ham->prms, qm_ham->modprms, qm_ham->atom_to_ao_map, qm_ham->ao_to_atom_map);

        /*
        // Finally, compute energy
        E_i = (energy_elec(qm_ham->el->P_alp, qm_ham->el->Hao, qm_ham->el->Fao_alp) +
               energy_elec(qm_ham->el->P_bet, qm_ham->el->Hao, qm_ham->el->Fao_bet) 
              );

        cout<<"E_i = "<<E_i<<endl;
 

        cout<<"Resetting to ground state, time = "<<tim1.stop()<<endl;


        tim1.start();
        // Forces and nuclear contributions:
        for(int n=0;n<_syst->Number_of_atoms;n++){

          _syst->Atoms[n].Atom_RB.rb_force = 
          force(*qm_ham->el, *_syst, qm_ham->basis_ao, qm_ham->prms, qm_ham->modprms,
                qm_ham->atom_to_ao_map, qm_ham->ao_to_atom_map, *qm_ham->el->Hao, *qm_ham->el->Sao,
                qm_ham->el->Norb, n, x_period, y_period, z_period, t1, t2, t3);

        }// for n


        // - nuclear-nuclear repulsion
        vector<double> Zeff;
        vector<VECTOR> G;
        vector<VECTOR> R;
        for(n=0;n<_syst->Number_of_atoms;n++){
          R.push_back(_syst->Atoms[n].Atom_RB.rb_cm);
          Zeff.push_back(qm_ham->modprms.PT[_syst->Atoms[n].Atom_element].Zeff);
        }// for n

        double Enucl = energy_nucl(R, Zeff, G);
        for(n=0;n<_syst->Number_of_atoms;n++){   _syst->Atoms[n].Atom_RB.rb_force -= G[n];   }

        E_i += Enucl;


        // Additive contribution - so to potentially include MM part
        ham_adi->M[i*nelec+i] += E_i;


        for(int n=0;n<_syst->Number_of_atoms;n++){     

          // Additive contribution - so to potentially include MM part
          d1ham_adi[3*n  ]->M[i*nelec+i] -= _syst->Atoms[n].Atom_RB.rb_force.x;
          d1ham_adi[3*n+1]->M[i*nelec+i] -= _syst->Atoms[n].Atom_RB.rb_force.y;
          d1ham_adi[3*n+2]->M[i*nelec+i] -= _syst->Atoms[n].Atom_RB.rb_force.z;


        }
        */
        cout<<"Nuclear-nuclear calculations, time = "<<tim1.stop()<<endl;
      
      
 
      }// method_option == 1

    }// ham_type[1]==1


    // Set status flag
    status_adi = 1;

  }// status_adi == 0


  //cout<<"finishing Hamiltonian_Atomistic::compute_adiabatic()\n";
  
}





void Hamiltonian_Atomistic::init_qm_Hamiltonian(std::string ctrl_filename){
/**
  Initialize the QM Hamiltonian using control settings in the file
  Note the file will also refer to another file - the one with Hamiltonian parameters, so
  actually we would need to prepare 2 files: one - the control parameters file, ctrl_filename, 
  the other - the parameters file.
  This function will read the file, initialize all necessary variables (except for the coordinates)
  and performe an SCF calculations, to compute electronic structure and related properties. 

  Note: the molecular (chemical) system object we are studying must be created beforehand (e.g. loaded
  from yet another file) and the object must be bound the the atomistic Hamiltonian object - this all
  has to be done before calling the present function

  \param[in] ctrl_filename The file containing all the control settings for the QM calculations
*/

  if(_syst==NULL){ cout<<"Error: System is not initialized yet\n"; exit(0); }
  
  else{
    qm_ham->init(ctrl_filename, *_syst);
  }

}


void Hamiltonian_Atomistic::excite_alp(int I, int J){
/**
  \param[in] I the index of the source orbital involved in the alpha-excitation
  \param[in] J the index of the target orbital involved in the alpha-excitation

  This function creates an excitation - the electronic excited state - the basis state for NA-MD
  This excitation is for the alpha electron going from orbital with index I to that with index J (I-->J)  
*/

  qm_ham->excite_alp(I,J);

  status_dia = 0;
  status_adi = 0;

}

void Hamiltonian_Atomistic::excite_bet(int I, int J){
/**
  \param[in] I the index of the source orbital involved in the beta-excitation
  \param[in] J the index of the target orbital involved in the beta-excitation

  This function creates an excitation - the electronic excited state - the basis state for NA-MD
  This excitation is for the beta electron going from orbital with index I to that with index J (I-->J)  
*/

                 
  qm_ham->excite_bet(I,J);

  status_dia = 0;
  status_adi = 0;

}

MATRIX3x3 Hamiltonian_Atomistic::get_stress(std::string opt){

  MATRIX3x3 res; res = 0.0;
  
  if(ham_types[0]==1){
    int nint = mm_ham->active_interactions.size();

    if(opt=="at"){
      for(int i=0;i<nint;i++){
        if(mm_ham->interactions[mm_ham->active_interactions[i]].get_status()){
          res += mm_ham->interactions[mm_ham->active_interactions[i]].stress_at; 
        }
      }
    }// at

    else if(opt=="fr"){
      for(int i=0;i<nint;i++){
        if(mm_ham->interactions[mm_ham->active_interactions[i]].get_status()){
          res += mm_ham->interactions[mm_ham->active_interactions[i]].stress_fr; 
        }
      }
    }// fr

    else if(opt=="ml"){
      for(int i=0;i<nint;i++){
        if(mm_ham->interactions[mm_ham->active_interactions[i]].get_status()){
          res += mm_ham->interactions[mm_ham->active_interactions[i]].stress_ml; 
        }
      }
    }// ml

  }// MM Hamiltonians


  return res;
}




}// namespace libatomistic

}// liblibra



