/*********************************************************************************
* Copyright (C) 2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_variables_electronic.cpp
  \brief The file implements the methods for electronic dyn vars
*/

#include "../util/libutil.h"
#include "../converters/libconverters.h"
#include "dyn_variables.h"
#include "dyn_ham.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libutil;
using namespace libconverters;


/// libdyn namespace
namespace libdyn{

namespace bp = boost::python;



CMATRIX transform_amplitudes(int rep_in, int rep_out, CMATRIX& C, nHamiltonian& ham){
/**
  This function converts the amplitudes from one representation to another

  The reason: we may be solving TD-SE (computing forces) in one representation
  but compute the hopping probabilities in another one.

  This function assumes we already have the basis transformation matrix in ham object 
  computed/updated

*/

  int nst = C.n_rows;    
  int ntraj = C.n_cols;

  CMATRIX Coeff(nst,ntraj); 

  /// Depending on the basis, select which   
  /// C - the basis in which the electron-nuclear propagation is done
  /// Coeff - the basis in which SH is done

  // Input in the diabatic basis
  if(rep_in==0){                   
    // Output in the diabatic basis too
    if(rep_out==0){    Coeff = C;   }

    // Output in the adiabatic basis
    else if(rep_out==1){  ham.ampl_dia2adi(C, Coeff, 0, 1);  }
  }

  // Input in the adiabatic basis 
  else if(rep_in==1){   
    // Output in the diabatic basis 
    if(rep_out==0){  ham.ampl_adi2dia(Coeff, C, 0, 1);   } 

    // Output in the diabatic basis too
    else if(rep_out==1){     Coeff = C;    } 
  }

  return Coeff;
}


void dyn_variables::update_amplitudes(dyn_control_params& dyn_params, nHamiltonian& ham){

  // The diabatic rep is the major, so we update the adiabatic amplitudes
  if(dyn_params.rep_tdse==0){   ham.ampl_dia2adi(ampl_dia, ampl_adi, 0, 1);  }

  // The adiabatic rep is the major, so we update the diabatic amplitudes
  else if(dyn_params.rep_tdse==1){   ham.ampl_adi2dia(ampl_dia, ampl_adi, 0, 1);  }

}

void dyn_variables::update_amplitudes(bp::dict dyn_params, nHamiltonian& ham){

  dyn_control_params _prms;
  _prms.set_parameters(dyn_params);

  update_amplitudes(_prms, ham);

}


void dyn_variables::update_amplitudes(dyn_control_params& dyn_params, bp::object compute_model, bp::dict model_params){
/**
    """
    Args:
        dyn_params ( Python dictionary ): control of the dynamics
        compute_model ( Python function ): the function that does the calculations
        model_params ( Python dictionary ): parameters of the computational model

    Returns:
        None: updates ampl_dia or ampl_adi, depending on `dyn_params.rep_tdse` variable:

        rep

    """
*/

    //# Prepare the Hamiltonian's hierarchy
    nHamiltonian ham(ndia, nadi, ndof);
    ham.add_new_children(ndia, nadi, ndof, ntraj);
    ham.init_all(2,1);

    //# Compute the Hamiltonian transformation
    //update_Hamiltonian_q(dyn_params, *q, ham, compute_model, model_params);
    update_Hamiltonian_variables(dyn_params, *this, ham, ham, compute_model, model_params, 0);


    // The diabatic rep is the major, so we update the adiabatic amplitudes
    if(dyn_params.rep_tdse==0){   ham.ampl_dia2adi(ampl_dia, ampl_adi, 0, 1);  }

    // The adiabatic rep is the major, so we update the diabatic amplitudes
    else if(dyn_params.rep_tdse==1){   ham.ampl_adi2dia(ampl_dia, ampl_adi, 0, 1);  }

}


void dyn_variables::update_amplitudes(bp::dict dyn_params, bp::object compute_model, bp::dict model_params){

  dyn_control_params _prms;
  _prms.set_parameters(dyn_params);

  update_amplitudes(_prms, compute_model, model_params);

}




void dyn_variables::update_density_matrix(dyn_control_params& dyn_params, nHamiltonian& ham, int lvl){
/**
    """

    Compute the trajectory-averaged density matrices in diabatic
    or adiabatic representations

    Args: 
        ham ( nHamiltonian ): object that handles Hamiltonian-related calculations with many trajectories
        Cdia ( CMATRIX(ndia, ntraj) ): amplitudes of diabatic states in the TD wavefunction expansion
        Cadi ( CMATRIX(ndia, ntraj) ): amplitudes of adiabatic states in the TD wavefunction expansion
        projectors ( list of CMATRIX(nst, nst)): dynamically-consistent corrections   
        rep ( int ): a selector of which representation is considered main (being propagated)
            E.g. if rep = 0 - that means we propagate the diabatic coefficients, that is the calculation 
            of the diabatic density matrix is straightforward, but we need to involve some transformations 
            to compute the adiabatic density matrix and vice versa, if rep = 1, the propagation is done according
            to the adiabatic properties and we'd need to convert to the diabatic representation in the end
             
            - 0: diabatic
            - 1: adiabatic

        lvl ( int ): The level of the Hamiltonian that treats the transformations:
            - 0: ham is the actual Hamiltonian to use (use with single trajectory),
            - 1: ham is the parent of the Hamiltonians to use (use with multiple trajectories)

        isNBRA ( int ): The flag for NBRA type calculations:
            - 0: the Hamiltonian related properties are computed for all of the trajectories [default]
            - 1: the Hamiltonian related properties are computed only for one trajectory 

    Returns:
        tuple: ( dm_dia, dm_adi ):

            * dm_dia ( CMATRIX(ndia, ndia) ): the trajectory-averaged density matrix in
                the diabatic representation. Here, ndia - is the number of diabatic basis
                states
            * dm_adi ( CMATRIX(nadi, nadi) ): the trajectory-averaged density matrix in
                the adiabatic representation. Here, nadi - is the number of adiabatic basis
                states

    """
*/

  CMATRIX su(ndia, nadi);
  CMATRIX S(ndia, ndia);
  CMATRIX U(ndia, nadi);
  CMATRIX cd(ndia, 1);
  CMATRIX ca(nadi, 1);

  //cout<<"In update density matrix\n";
  vector<int> indx;    
  if(lvl==0){ indx = vector<int>(1, 0); }
  else if(lvl==1) { indx = vector<int>(2, 0); } 

  for(int traj=0; traj<ntraj; traj++){

    if(lvl==1){  indx[1] = traj; }
  
    S = ham.get_ovlp_dia(indx);
    U = ham.get_basis_transform(indx);

    /// Diabatic wfc representation
    if(dyn_params.rep_tdse==0){
      cd = ampl_dia->col(traj);
      *dm_dia[traj] = S * (cd * cd.H()) * S; 
      *dm_adi[traj] = U.H() * (*dm_dia[traj]) * U;
    }  
    /// Adiabatic wfc representation
    else if(dyn_params.rep_tdse==1){  
      ca = ampl_adi->col(traj);
      *dm_adi[traj] = ca * ca.H();            
      su = S * U; 
      *dm_dia[traj] =  su * (*dm_adi[traj]) * su.H();       
    }
    /// Diabatic DM representation
    else if(dyn_params.rep_tdse==2){
      *dm_adi[traj] = U.H() * (*dm_dia[traj]) * U;
    }
    /// Adiabatic DM representation
    else if(dyn_params.rep_tdse==3){
      su = S * U;
      *dm_dia[traj] =  su * (*dm_adi[traj]) * su.H();
    }

  }// for traj

  
}

void dyn_variables::update_density_matrix(bp::dict dyn_params, nHamiltonian& ham, int lvl){

  dyn_control_params _prms;
  _prms.set_parameters(dyn_params);

  update_density_matrix(_prms, ham, lvl);

}



void dyn_variables::update_density_matrix(dyn_control_params& dyn_params, bp::object compute_model, bp::dict model_params, int lvl){

  //# Prepare the Hamiltonian's hierarchy
  nHamiltonian ham(ndia, nadi, ndof);
  ham.add_new_children(ndia, nadi, ndof, ntraj);
  ham.init_all(2,1);

  //# Compute the Hamiltonian transformation
  update_Hamiltonian_variables(dyn_params, *this, ham, ham, compute_model, model_params, 0);
  //update_Hamiltonian_q(dyn_params, *q, ham, compute_model, model_params);
  update_density_matrix(dyn_params, ham, lvl);

}


void dyn_variables::update_density_matrix(bp::dict dyn_params, bp::object compute_model, bp::dict model_params, int lvl){

  dyn_control_params _prms;
  _prms.set_parameters(dyn_params);


  //# Prepare the Hamiltonian's hierarchy
  nHamiltonian ham(ndia, nadi, ndof);
  ham.add_new_children(ndia, nadi, ndof, ntraj);
  ham.init_all(2,1);

  //# Compute the Hamiltonian transformation
  //update_Hamiltonian_q(_prms, *q, ham, compute_model, model_params);
  update_Hamiltonian_variables(_prms, *this, ham, ham, compute_model, model_params, 0);

  update_density_matrix(_prms, ham, lvl);

}


//bp::dict dyn_variables::init_electronic_dyn_var(bp::dict& _params, Random& rnd){
void dyn_variables::init_amplitudes(bp::dict _params, Random& rnd){
/**
    """
    Args:

        params ( dictionary ): control parameters
 
            * **params["init_type"]** ( int ): the type of sampling of electronic DOFs
     
                - 0 : initialize all states according to "istate" and sets all
                    amplitudes to 1.0  [ default ]
 
                - 1 : initialize all states according to "istate" but sets 
                    amplitudes to exp(-2*pi*i*rnd), where rnd is a random 
                    number uniformly distributed on the [0, 1] interval

                - 2 : initialize all states according to "istates" - the 
                    integer indices are selected randomly according to populations provided in
                    variable "istates", the amplitudes are set to be sqrt(istates[i]), but
                    their phases are identical (like in the option 0)

                - 3 : initialize all states according to "istates" - the 
                    integer indices are selected randomly according to populations provided in
                    variable "istates", the amplitudes are set to be sqrt(istates[i]), but
                    their phases are set to exp(-2*pi*i*rnd), where rnd is a random 
                    number uniformly distributed on the [0, 1] interval (line in option 1)

            * **params["nstates"]** ( int ): the number of electronic states in the basis
                [ default: 1 ]

            * **params["istate"]** ( int ): the index of the initial electronic state, used 
                only when **params["init_type"]** is 0 or 1, in which case it defines on 
                which state the amplitudes will be initialized [ default: 0]

            * **params["istates"]** ( list of ints ): the list of the populations on all electronic
                states included in the basis, used only when **params["init_type"]** is 2 or 3, 
                in which case it defines on which states the amplitudes will be initialized.
                The length of this list should be consistent with ```nstates``` variable. And the sum
                of all entries should be 1.0 [ default: [1.0], meaning that only the lowest state
                is occupied ]

            * **params["rep"]** ( int ): defines for which repersentation we generate the amplitudes
 
                - 0 : diabatic, the corresponding matrix is initialized, the adiabatic amplitudes are zeroes
                - 1 : adiabatic, the corresponding matrix is initialized, the diabatic amplitudes are zeroes [ default ]

            * **params["ntraj"]** ( int ): the number of trajectories - the parameter defines the
                number of columns the output matrices will have [ default: 1]

            * **params["is_nbra"]** (int): A flag for NBRA type calculations. If it is set to 1 then
                                          the Hamiltonian related properties are only computed for one trajectory [ default : 0]

        rnd ( Random ): random numbers generator object


    Returns:
        Cdia, Cadi, states:  where:

            * Cdia ( CMATRIX(nstates, ntraj) ) : amplitudes on all diabatic states for all trajectories
            * Cadi ( CMATRIX(nstates, ntraj) ) : amplitudes on all adiabatic states for all trajectories
            * states ( list of ints ) : state indices for each trajectory

    """
*/

  //# Read the parameters
  bp::list critical_params; 
  bp::dict default_params;
  bp::dict params(_params);

  default_params["init_type"] = 0;
  default_params["istate"] = 0;
  bp::list lst; lst.append(1.0);
  default_params["istates"] = lst;
  default_params["rep"] = 1;
  default_params["is_nbra"] = 0;
  default_params["verbosity"] = 0;

  check_input(params, default_params, critical_params);


  int i, traj;
  double ksi;
  int init_type;
  int istate;
  vector<double> istates;
  int rep;
  int is_nbra; 
  int verbosity;

  std::string key;
  for(int i=0;i<len(params.values());i++){
    key = bp::extract<std::string>(params.keys()[i]);

    ///================= Computing Hamiltonian-related properties ====================
    if(key=="init_type") {  init_type = bp::extract<int>(params.values()[i]); }
    else if(key=="istate") {  istate = bp::extract<int>(params.values()[i]); }
    else if(key=="istates") {  istates = liblibra::libconverters::Py2Cpp<double>( bp::extract< bp::list >(params.values()[i]) ); }
    else if(key=="rep") {  rep = bp::extract<int>(params.values()[i]); }
    else if(key=="is_nbra") {  is_nbra = bp::extract<int>(params.values()[i]); }
    else if(key=="verbosity") {  verbosity = bp::extract<int>(params.values()[i]); }

  }


  ///# Sanity check
  if(! ((rep==0) || (rep==1)) ){
    cout<<"WARNINIG in init_amplitudes:\
           the rep = "<<rep<<" is not known. Allowed values are: [0, 1]\n";
  }

  if(! (init_type==0 || init_type==1 || init_type==2 || init_type==3)){
    cout<<"WARNINIG in init_amplitudes:\
           the init_type = "<<init_type<<" is not known. Allowed values are: [0, 1, 2, 3]\n";
  }
  
  if(init_type==0 || init_type==1){
    if(istate>=ndia && rep==0){
      cout<<"ERROR in init_amplitudes: the istate is= "<<istate<<", but should be less than "<<ndia<<"\n";
      exit(0);
    }
    if(istate>=nadi && rep==1){
      cout<<"ERROR in init_amplitudes: the istate is= "<<istate<<", but should be less than "<<nadi<<"\n";
      exit(0);
    }
  }

  if(init_type==2 || init_type==3){
    if(istates.size()!=ndia && rep==0){
      cout<<"ERROR in init_amplitudes: the istates array is of length = "<<istates.size()<<", but should be of length "<<ndia<<"\n";
      exit(0);
    }
    if(istates.size()!=nadi && rep==1){
      cout<<"ERROR in init_amplitudes: the istates array is of length = "<<istates.size()<<", but should be of length "<<nadi<<"\n";
      exit(0);
    }

    double summ = 0.0;
    for(i=0; i<istates.size(); i++){  summ += istates[i]; }
    if( fabs(summ - 1.0) > 1e-5 ){
      cout<<"ERROR in init_amplitudes: the sum of the entries in the istates array is "<<summ<<", but should be 1.0\n";
      exit(0);
    }
  }

  //# Dynamical variables
  for(traj=0; traj<ntraj; traj++){

    //========================== First, let's set up the amplitudes =======================================
    if(init_type==0){
      if(verbosity > 0){
        cout<<"======= Initialization type is "<<init_type<<" ========\n";
        cout<<"setting representation "<<rep<<" coefficient C_"<<istate<<" to 1.0\n";
      }

      if(rep==0){ ampl_dia->set(istate, traj, complex<double>(1.0, 0.0) ); }  
      else if(rep==1){ ampl_adi->set(istate, traj, complex<double>(1.0, 0.0) );  }

    }// init_type==0

    else if(init_type==1){
      if(verbosity > 0){
        cout<<"======= Initialization type is "<<init_type<<" ========\n";
        cout<<"setting representation "<<rep<<" coefficient to complex C_"<<istate<<" such that |C_"<<istate<<"}|^2 = 1.0\n";
      }

      ksi = rnd.uniform(0.0, 1.0);
      complex<double> ampl(cos(2.0*M_PI*ksi), sin(2.0*M_PI*ksi) );

      if(rep==0){ ampl_dia->set(istate, traj, ampl );   }
      else if(rep==1){ ampl_adi->set(istate, traj, ampl );  }

    }// init_type==1

    else if(init_type==2){
      if(verbosity > 0){
        cout<<"======= Initialization type is "<<init_type<<" ========\n";
        cout<<"setting representation "<<rep<<" coefficients C_i for all i to sqrt( target populations) \n";

      }

      for(i=0; i<istates.size(); i++){
          double ampl = sqrt( istates[i] );

          if(rep==0){ ampl_dia->set(i, traj, complex<double>(ampl, 0.0) );  }
          else if(rep==1){ ampl_adi->set(i, traj, complex<double>(ampl, 0.0) );  }
      }

    }// init_type==2

    else if(init_type==3){
      if(verbosity > 0){
        cout<<"======= Initialization type is "<<init_type<<" ========\n";
        cout<<"setting representation "<<rep<<" coefficients C_i for all i to complex numbers such that |C_i|^2  = target populations \n";

      }

      for(i=0; i<istates.size(); i++){
          ksi = rnd.uniform(0.0, 1.0);
          complex<double> ampl(cos(2.0*M_PI*ksi), sin(2.0*M_PI*ksi) );
          ampl = ampl * sqrt( istates[i] );

          if(rep==0){ ampl_dia->set(i, traj, ampl );  }
          else if(rep==1){ ampl_adi->set(i, traj, ampl);  }
      }

    }// init_type==2
  }// for traj

  if(verbosity > 1){
      cout<<"========== ampl_dia ===============\n";
      ampl_dia->show_matrix();
      cout<<"========== ampl_adi ===============\n";
      ampl_adi->show_matrix();
      cout<<"========== states ===============\n";
  }



  //return params;

}



void dyn_variables::init_density_matrix(bp::dict _params){

  //# Read the parameters
  bp::list critical_params; 
  bp::dict default_params;
  bp::dict params(_params);

  default_params["init_dm_type"] = 0; /// 0 - based on amplitudes; 1 - using direct input
  default_params["verbosity"] = 0;
  default_params["extern_dm"] = CMATRIX(nadi, nadi); /// 0 - based on amplitudes; 1 - using direct input


  check_input(params, default_params, critical_params);


  int i, traj;
  double ksi;
  int init_dm_type;
  int rep;
  int verbosity;
  //CMATRIX extern_dm_dia();

  std::string key;
  for(int i=0;i<len(params.values());i++){
    key = bp::extract<std::string>(params.keys()[i]);

    ///================= Computing Hamiltonian-related properties ====================
    if(key=="init_dm_type") {  init_dm_type = bp::extract<int>(params.values()[i]); }
    else if(key=="rep") {  rep = bp::extract<int>(params.values()[i]); }
    else if(key=="verbosity") {  verbosity = bp::extract<int>(params.values()[i]); }

  }


  ///================= Sanity check ================================
  if(! ((rep==0) || (rep==1)) ){
    cout<<"WARNINIG in init_density_matrix:\
           the rep = "<<rep<<" is not known. Allowed values are: [0, 1]\n";
  }

  if(! (init_dm_type==0 || init_dm_type==1 )){
    cout<<"WARNINIG in init_density_matrix:\
           the init_dm_type = "<<init_dm_type<<" is not known. Allowed values are: [0, 1]\n";
  }




  ///================= Actual calculations  ================================

  //# Dynamical variables
  for(traj=0; traj<ntraj; traj++){

    // From the amplitudes
    if(init_dm_type==0){
      if(verbosity > 0){
        cout<<"======= Initialization type is "<<init_dm_type<<" ========\n";
        cout<<"setting DM in representation "<<rep<<" using wfc amplitudes \n";
      }

      if(rep==0){ 
        CMATRIX C(ndia, 1); 
        C = ampl_dia->col(traj);
        *dm_dia[traj] = C * C.H(); 
      }  
      else if(rep==1){ 
        CMATRIX C(nadi, 1); 
        C = ampl_adi->col(traj);
        *dm_adi[traj] = C * C.H(); 
      }

    }// init_dm_type==0


    if(init_dm_type==1){
      if(verbosity > 0){
        cout<<"======= Initialization type is "<<init_dm_type<<" ========\n";
        cout<<"setting DM in representation "<<rep<<" using provided input \n";
      }

      if(rep==0){ 
        for(int i=0;i<len(params.values());i++){
          key = bp::extract<std::string>(params.keys()[i]);
          if(key=="extern_dm") {  *dm_dia[traj] = bp::extract<CMATRIX>(params.values()[i]); }
        }// for i
      }// rep == 0 
      else if(rep==1){ 
        for(int i=0;i<len(params.values());i++){
          key = bp::extract<std::string>(params.keys()[i]);
          if(key=="extern_dm") {  *dm_adi[traj] = bp::extract<CMATRIX>(params.values()[i]); }
        }// for i
      }// rep == 1 

    }// init_dm_type == 1
 
  }// for traj
  

  //return params;
}


void dyn_variables::init_active_states(bp::dict _params, Random& rnd){

  //# Read the parameters
  bp::list critical_params; 
  bp::dict default_params;
  bp::dict params(_params);


  default_params["init_type"] = 0;
  default_params["istate"] = 0;
  bp::list lst; lst.append(1.0);
  default_params["istates"] = lst;
  default_params["rep"] = 1;
  default_params["verbosity"] = 0;

  check_input(params, default_params, critical_params);

  int i, traj;
  double ksi;
  int init_type;
  int istate;
  vector<double> istates;
  int rep;
  int verbosity;


  std::string key;
  for(int i=0;i<len(params.values());i++){
    key = bp::extract<std::string>(params.keys()[i]);

    ///================= Computing Hamiltonian-related properties ====================
    if(key=="init_type") {  init_type = bp::extract<int>(params.values()[i]); }
    else if(key=="istate") {  istate = bp::extract<int>(params.values()[i]); }
    else if(key=="istates") {  istates = liblibra::libconverters::Py2Cpp<double>( bp::extract< bp::list >(params.values()[i]) ); }
    else if(key=="rep") {  rep = bp::extract<int>(params.values()[i]); }
    else if(key=="verbosity") {  verbosity = bp::extract<int>(params.values()[i]); }

  }


  ///================= Sanity check ================================
  if(! ((rep==0) || (rep==1)) ){
    cout<<"WARNINIG in init_active_states:\
           the rep = "<<rep<<" is not known. Allowed values are: [0, 1]\n";
  }

  if(! (init_type==0 || init_type==1 || init_type==2 || init_type==3)){
    cout<<"WARNINIG in init_active_states:\
           the init_type = "<<init_type<<" is not known. Allowed values are: [0, 1, 2, 3]\n";
  }
  
  if(init_type==0 || init_type==1){
    if(istate>=ndia && rep==0){
      cout<<"ERROR in init_active_states: the istate is= "<<istate<<", but should be less than "<<ndia<<"\n";
      exit(0);
    }
    if(istate>=nadi && rep==1){
      cout<<"ERROR in init_active_states: the istate is= "<<istate<<", but should be less than "<<nadi<<"\n";
      exit(0);
    }
  }

  if(init_type==2 || init_type==3){
    if(istates.size()!=ndia && rep==0){
      cout<<"ERROR in init_active_states: the istates array is of length = "<<istates.size()<<", but should be of length "<<ndia<<"\n";
      exit(0);
    }
    if(istates.size()!=nadi && rep==1){
      cout<<"ERROR in init_active_states: the istates array is of length = "<<istates.size()<<", but should be of length "<<nadi<<"\n";
      exit(0);
    }

    double summ = 0.0;
    for(i=0; i<istates.size(); i++){  summ += istates[i]; }
    if( fabs(summ - 1.0) > 1e-5 ){
      cout<<"ERROR in init_active_states: the sum of the entries in the istates array is "<<summ<<", but should be 1.0\n";
      exit(0);
    }
  }



  ///================= Actual calculations  ================================

  act_states.clear();  

  //# Dynamical variables
  for(traj=0; traj<ntraj; traj++){

    if(init_type==0 || init_type==1){ 
      act_states.push_back(istate);
    }
    else if(init_type==2 || init_type==3){
      ksi = rnd.uniform(0.0, 1.0);
      act_states.push_back( liblibra::libspecialfunctions::set_random_state(istates, ksi) );
    }

  }// for traj

  //return params;

}


void dyn_variables::init_electronic_dyn_var(bp::dict _params, Random& rnd){

//  bp::dict res;

//  return res;

}


//vector<int> update_active_states(vector<int>& act_states, vector<CMATRIX*>& T){
void dyn_variables:: update_active_states(){
/*
   act_states[itraj] - index of the active adiabatic state for the trajectory `itraj`
                       in terms of the energy-ordered (instantaneous) eign-states
*/

  //cout<<"========= in update_active_states=============\n";
  //int ntraj = act_states.size();
  //int nst = T[0]->n_cols;
  //vector<int> res(ntraj, 0);

  int nst = nadi;

  CMATRIX t2(nst,nst);
  CMATRIX t2_half(nst, nst);
  CMATRIX t2_i_half(nst, nst);
  CMATRIX coeff(nst, 1);
  complex<double> max_val;

  for(int itraj=0; itraj<ntraj; itraj++){
    //cout<<"itraj = "<<itraj<<endl;

    CMATRIX t(*proj_adi[itraj]);
    t2 = t.H() * t;
    sqrt_matrix(t2, t2_half, t2_i_half);
    t2 = t * t2_i_half;

    coeff = 0.0;
    coeff.set(act_states[itraj], 0, complex<double>(1.0, 0.0));
    //cout<<"initial coeff = \n"; coeff.show_matrix();
    //cout<<"t2 = \n"; t2.show_matrix();
    coeff = t2 * coeff;
    //cout<<"new coeff = \n"; coeff.show_matrix();
    coeff.max_col_elt(0, max_val, act_states[itraj]);
    //cout<<" max_val = "<<max_val<<" in position = "<<act_states[itraj]<<endl;
    //cout<<"itraj = "<<itraj<<"  "<<act_states[itraj]<<endl; //#<<" -> "<<res[itraj]<<endl;

  }// for itraj

  //cout<<"============ done with the update ================\n";

//  return res;
}

CMATRIX orthogonalized_T(CMATRIX& T){
/*
*/

  int nst = T.n_cols;

  CMATRIX t2(nst,nst);
  CMATRIX t2_half(nst, nst);
  CMATRIX t2_i_half(nst, nst);
    
  t2 = T.H() * T;
  sqrt_matrix(t2, t2_half, t2_i_half);
  t2 = T * t2_i_half;

  return t2;

}

CMATRIX dyn_variables::compute_average_dm(int rep){

  int sz;
  if(rep==0 || rep==2){ sz = ndia; }
  else if(rep==1 || rep==3){ sz = nadi; }

  CMATRIX res(sz, sz);

  for(int traj=0; traj<ntraj; traj++){
    if(rep==0 || rep==2){   res += *dm_dia[traj]; }
    else if(rep==1 || rep==3){ res += *dm_adi[traj]; }
  }

  res = res / (float)ntraj;

  return res;
}

  
vector<double> dyn_variables::compute_average_se_pop(int rep){

  int sz; 
  if(rep==0 || rep==2){ sz = ndia; }
  else if(rep==1 || rep==3){ sz = nadi; }

  MATRIX ave(sz, sz);
  ave = compute_average_dm(rep).real();  

  vector<double> res(sz, 0.0);

  for(int i=0; i<sz; i++){ res[i] = ave.get(i,i); }

  return res;

}
  

vector<double> dyn_variables::compute_average_sh_pop(){

  int sz = nadi; 

  //cout<<"In compute_average_sh_pop...\n";

  vector<double> res(sz, 0.0);
//  vector<int> effective_states( update_active_states(act_states, proj_adi) );
  vector<int> effective_states( act_states );

  for(int traj=0; traj<ntraj; traj++){

    int i = effective_states[traj];
    res[i] += 1.0;
    //cout<<" traj = "<<traj<<" act_state = "<<act_states[traj]<<" eff state = "<<i<<endl;
  }
  
//  cout<<" Populations : ";
  for(int j=0; j<sz; j++){   
    res[j] = res[j] / (float)ntraj; 
  //  cout<<res[j]<<" ";
  }
  //cout<<endl;

  return res;
}


void dyn_variables::save_curr_dm_into_prev(){

  if(fssh2_vars_status==1){
    for(int i=0; i<ntraj; i++){
      *dm_dia_prev[i] = *dm_dia[i];
      *dm_adi_prev[i] = *dm_adi[i];
    }
  }

}


}// namespace libdyn
}// liblibra

