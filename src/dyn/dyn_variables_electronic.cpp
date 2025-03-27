/*********************************************************************************
* Copyright (C) 2022-2024 Alexey V. Akimov
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

void dyn_variables::update_basis_transform(nHamiltonian& Ham){
/**

*/
  for(int itraj=0; itraj<ntraj; itraj++){ 
    nHamiltonian* ham = Ham.children[itraj];
    *basis_transform[itraj] = *(ham->basis_transform);

  }// for itraj
}


void dyn_variables::update_amplitudes(dyn_control_params& dyn_params){
/**
           |\psi_adi > = | \psi_dia > U

   |\Psi> = |\psi_adi > C_adi = |\psi_dia> U C_adi = |\psi_dia> C_dia, so:

   C_dia = U C_adi

   C_adi = U.H() C_dia
*/
  int i,j;
  CMATRIX cd(ndia, 1);
  CMATRIX ca(nadi, 1);

  for(int traj=0; traj<ntraj; traj++){

    CMATRIX& U = (*basis_transform[traj]);

    /// Diabatic wfc representation
    if(dyn_params.rep_tdse==0){
      cd = ampl_dia->col(traj);
      ca = U.H() * cd;
      for(i=0;i<nadi;i++){ ampl_adi->set(i, traj, ca.get(i,0) ); }
    }
    /// Adiabatic wfc representation
    else if(dyn_params.rep_tdse==1){
      ca = ampl_adi->col(traj);
      cd = U * ca; 
      for(i=0;i<ndia;i++){ ampl_dia->set(i, traj, cd.get(i,0) ); }
    }

  }// for traj

}


void dyn_variables::update_density_matrix(dyn_control_params& dyn_params){
/**
    """

    Update the density matrices in diabatic and adiabatic representations

    Args:
        * **dyn_params** - contains the parameters of transformation

        rep ( int ): a selector of which representation is considered main (being propagated)
            E.g. if rep = 0 - that means we propagate the diabatic coefficients, that is the calculation
            of the diabatic density matrix is straightforward, but we need to involve some transformations
            to compute the adiabatic density matrix and vice versa, if rep = 1, the propagation is done according
            to the adiabatic properties and we'd need to convert to the diabatic representation in the end

            - 0: diabatic
            - 1: adiabatic

    Notes:
           Since the basis_transform corresponding to the transformation of the diabatic-to-raw-adiabatic states: 
           |\psi_adi (t) > = | \psi_dia (t)> U
           |\psi\tilde_adi (t) > = |\psi_adi (t) > T =  | \psi_dia (t)> U * T, so the effective U is U_eff = U * T
          
           T is the basis re-projection matrix

    This function should be used only after the proj_adi = T has been computed and the basis_transform has been set

    """
*/

  CMATRIX cd(ndia, 1);
  CMATRIX ca(nadi, 1);

  // Assuming that in the diabatic representation, the overlap 
  // matrix is identity S = I
  for(int traj=0; traj<ntraj; traj++){

    CMATRIX& U = (*basis_transform[traj]); // * (*proj_adi[indx]); // U_eff - let's wait with this

    /// Diabatic wfc representation
    if(dyn_params.rep_tdse==0){
      cd = ampl_dia->col(traj);
      *dm_dia[traj] = (cd * cd.H());
      *dm_adi[traj] = U.H() * (*dm_dia[traj]) * U;
    }
    /// Adiabatic wfc representation
    else if(dyn_params.rep_tdse==1){
      ca = ampl_adi->col(traj);
      *dm_adi[traj] = ca * ca.H();
      *dm_dia[traj] =  U * (*dm_adi[traj]) * U.H();
    }
    /// Diabatic DM representation
    else if(dyn_params.rep_tdse==2){
      *dm_adi[traj] = U.H() * (*dm_dia[traj]) * U;
    }
    /// Adiabatic DM representation
    else if(dyn_params.rep_tdse==3){
      *dm_dia[traj] =  U * (*dm_adi[traj]) * U.H();
    }

  }// for traj

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

                - 4 : initialize all states according to Voronoi-style sampling 
                    Designed for MASH dynamics
                    Uses "istates" instead of the "istate", similar to the types 2 and 3

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

  if(! (init_type==0 || init_type==1 || init_type==2 || init_type==3 || init_type==4)){
    cout<<"WARNINIG in init_amplitudes:\
           the init_type = "<<init_type<<" is not known. Allowed values are: [0, 1, 2, 3, 4]\n";
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

    }// init_type==3

    else if(init_type==4){
      if(verbosity > 0){
        cout<<"======= Initialization type is "<<init_type<<" ========\n";
        cout<<"setting representation "<<rep<<" coefficients C_i for all i to complex numbers such that populations are sampled evenly \n";
      }

      istates = rnd.voron(istates.size(), istate);
      for(i=0; i<istates.size(); i++){
        cout<<"Voron_mag["<<i<<"]= "<<istates[i]<<endl;
        ksi = rnd.uniform(0.0, 1.0);
        complex<double> ampl(cos(2.0*M_PI*ksi), sin(2.0*M_PI*ksi) );
        ampl = ampl * sqrt( istates[i] );
        cout<<"Voron_amp["<<i<<"]= "<<ampl<<endl;

        if(rep==0){ ampl_dia->set(i, traj, ampl );  }
        else if(rep==1){ ampl_adi->set(i, traj, ampl);  }
      }
    }

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

  if(! (init_type==0 || init_type==1 || init_type==2 || init_type==3 || init_type==4) ){
    cout<<"WARNINIG in init_active_states:\
           the init_type = "<<init_type<<" is not known. Allowed values are: [0, 1, 2, 3, 4]\n";
  }
  
  if(init_type==0 || init_type==1 || init_type==4){
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

  //================ For adiabatic - set up them directly =================
  for(traj=0; traj<ntraj; traj++){
    if(init_type==0 || init_type==1 || init_type==4){ 
      act_states.push_back(istate);
    }
    else if(init_type==2 || init_type==3){
      ksi = rnd.uniform(0.0, 1.0);
      act_states.push_back( liblibra::libspecialfunctions::set_random_state(istates, ksi) );
    }
  }// for traj

  //============= For adiabatic: more complicated =========================
  if(rep==0){
    vector<int> dia_act_states_temp(act_states);

    CMATRIX pop_adi(nadi, nadi);
    CMATRIX pop_dia(ndia, ndia);

    for(traj=0; traj<ntraj; traj++){
      i = dia_act_states_temp[traj]; // active diabatic state
      pop_dia *= 0.0; pop_dia.set(i, i, complex<double>(1.0, 0.0) );

      // The following transformation is correct only for S_dia = 1
      pop_adi = (*basis_transform[traj]).H() * pop_dia * (*basis_transform[traj]);

      vector<double> res(nadi, 0.0);
      for(i=0; i<nadi;i++){ res[i] = pop_adi.get(i,i).real(); }

      // Now select an adiabatic state according to the probabilities in the res
      ksi = rnd.uniform(0.0, 1.0);
      act_states[traj] = hop(res, ksi);  // this is the adiabatic state that corresponds to
                                         // the distribution of the diabatic populations


    }// for traj

  }// rep == 0

}


void dyn_variables::init_active_states_dia(bp::dict _params, Random& rnd){

  //# Read the parameters
  bp::list critical_params; 
  bp::dict default_params;
  bp::dict params(_params);


  default_params["init_type"] = 0;
  default_params["istate"] = 0;
  bp::list lst; lst.append(1.0);
  default_params["istates"] = lst;
  default_params["rep"] = 0;
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
    cout<<"WARNINIG in init_active_states_dia:\
           the rep = "<<rep<<" is not known. Allowed values are: [0, 1]\n";
  }

  if(! (init_type==0 || init_type==1 || init_type==2 || init_type==3 || init_type==4) ){
    cout<<"WARNINIG in init_active_states_dia:\
           the init_type = "<<init_type<<" is not known. Allowed values are: [0, 1, 2, 3, 4]\n";
  }
  
  if(init_type==0 || init_type==1 || init_type==4){
    if(istate>=ndia && rep==0){
      cout<<"ERROR in init_active_states_dia: the istate is= "<<istate<<", but should be less than "<<ndia<<"\n";
      exit(0);
    }
    if(istate>=ndia && rep==1){
      cout<<"ERROR in init_active_states_dia: the istate is= "<<istate<<", but should be less than "<<ndia<<"\n";
      exit(0);
    }
  }

  if(init_type==2 || init_type==3){
    if(istates.size()!=ndia && rep==0){
      cout<<"ERROR in init_active_states_dia: the istates array is of length = "<<istates.size()<<", but should be of length "<<ndia<<"\n";
      exit(0);
    }
    if(istates.size()!=ndia && rep==1){
      cout<<"ERROR in init_active_states_dia: the istates array is of length = "<<istates.size()<<", but should be of length "<<ndia<<"\n";
      exit(0);
    }

    double summ = 0.0;
    for(i=0; i<istates.size(); i++){  summ += istates[i]; }
    if( fabs(summ - 1.0) > 1e-5 ){
      cout<<"ERROR in init_active_states_dia: the sum of the entries in the istates array is "<<summ<<", but should be 1.0\n";
      exit(0);
    }
  }



  ///================= Actual calculations  ================================

  act_states_dia.clear();  

  //================ For diabatic - set up them directly =================
  for(traj=0; traj<ntraj; traj++){
    if(init_type==0 || init_type==1 || init_type==4){ 
      act_states_dia.push_back(istate);
    }
    else if(init_type==2 || init_type==3){
      ksi = rnd.uniform(0.0, 1.0);
      act_states_dia.push_back( liblibra::libspecialfunctions::set_random_state(istates, ksi) );
    }
  }// for traj

  //============= For diabatic: more complicated =========================
  if(rep==1){
    vector<int> act_states_temp(act_states_dia);

    CMATRIX pop_adi(nadi, nadi);
    CMATRIX pop_dia(ndia, ndia);

    for(traj=0; traj<ntraj; traj++){
      i = act_states_temp[traj]; // active diabatic state
      pop_adi *= 0.0; pop_adi.set(i, i, complex<double>(1.0, 0.0) );

      // The following transformation is correct only for S_dia = 1
      pop_dia = (*basis_transform[traj]) * pop_adi * (*basis_transform[traj]).H();

      vector<double> res(ndia, 0.0);
      for(i=0; i<ndia;i++){ res[i] = pop_dia.get(i,i).real(); }

      // Now select an adiabatic state according to the probabilities in the res
      ksi = rnd.uniform(0.0, 1.0);
      act_states_dia[traj] = hop(res, ksi);  // this is the adiabatic state that corresponds to
                                         // the distribution of the diabatic populations


    }// for traj

  }// rep == 1

}


void dyn_variables::init_electronic_dyn_var(bp::dict _params, Random& rnd){

//  bp::dict res;

//  return res;

}


//vector<int> update_active_states(vector<int>& act_states, vector<CMATRIX*>& T){
void dyn_variables:: update_active_states(int direction, int property){
/*
  |psi_adi_tilde> C_adi_tilde = |psi_adi> C_adi

  |psi_adi_tilde> = |psi_adi> T, so:

  T C_adi_tilde = C_adi and  C_adi_tilde = T^+ C_adi
 
   direction = 1 - forward:  C_adi_tilde = T^+ C_adi
              -1 - backward  C_adi = T C_adi_tilde

   property  = 0 - active state only
             = 1 - adiabatic amplitudes only 
             = 2 - both active states and adiabatic amplitudes

   act_states[itraj] - index of the active adiabatic state for the trajectory `itraj`
                       in terms of the energy-ordered (instantaneous) eign-states
*/

  int i;
  int nst = nadi;
  CMATRIX coeff(nst, 1);
  MATRIX rho_tmp(nst,nst);
  complex<double> max_val;

  for(int itraj=0; itraj<ntraj; itraj++){

    //================ Active states ===================
    if(property==0 || property==2){
      coeff = 0.0;
      coeff.set(act_states[itraj], 0, complex<double>(1.0, 0.0));

      if(direction==1){    coeff = proj_adi[itraj]->H() * coeff; }
      else if(direction==-1){  coeff = (*proj_adi[itraj]) * coeff; }

      //coeff.max_col_elt(0, max_val, act_states[itraj]);
      rho_tmp = (coeff * coeff.H()).real();
      //double max_val = fabs(rho_tmp.get(0,0)); 
      int max_val_indx = act_states[itraj];
      double max_val = fabs(rho_tmp.get(max_val_indx, max_val_indx));
      for(int j=0;j<nst;j++){ 
        // >= is important !!! not just >
        // on 2/15/2025:  actually, we need ">"
        if( fabs(rho_tmp.get(j,j)) > max_val){  max_val_indx = j; max_val = fabs(rho_tmp.get(j,j)); } 
      }
      act_states[itraj] = max_val_indx;
    }

    //============ Amplitudes of states ===================
    if(property==1 || property==2){
      coeff = ampl_adi->col(itraj);

      if(direction==1){    coeff = proj_adi[itraj]->H() * coeff; }
      else if(direction==-1){  coeff = (*proj_adi[itraj]) * coeff; }

      for(i=0; i<nst; i++){ ampl_adi->set(i, itraj,  coeff.get(i, 0) ); }
    }

  }// for itraj
}

void dyn_variables::update_active_states(){
/**
  for backward-compatability
*/
   update_active_states(1, 0);
}



void dyn_variables:: set_active_states_diff_rep(int rep_sh, Random& rnd){
/*
  Set the active state in the other representation through the transformation matrix

  rep_sh = 1 - Use the adiabatic active states to set the diabatic ones
           0 - Use the diabatic active states to set the adiabatic ones

  The process for determining the active state is based on the transformation of the density matrix,
  whose diagonal elements are replaced by the active state information, in similar to the following.
  Tempelaar, R.; Reichman, D. R. Generalization of Fewest-Switches Surface Hopping for Coherences. 
  The Journal of Chemical Physics 2018, 148 (10), 102309. https://doi.org/10.1063/1.5000843
*/
  int sz, i, j, traj;
  double ksi; 
  if(rep_sh==0){ sz = ndia; }
  else if(rep_sh==1){ sz = nadi; }
  else{ 
    cout<<"Wrong representation = "<<rep_sh<<"\nExiting...\n";
    exit(0);
  }
  
  // Transform the active states between representations based on the transformation matrix.
  CMATRIX pop_adi(nadi, nadi);
  CMATRIX pop_dia(ndia, ndia);
  CMATRIX U(ndia, nadi);
  
  // Obtain the diabatic active states from the adiabatic ones
  if(rep_sh==1){

    for(traj=0; traj<ntraj; traj++){
      i = act_states[traj]; // active adiabatic state
      pop_adi = *dm_adi[traj];
      for(j=0;j<nadi; j++){ pop_adi.set(j,j, complex<double>(0.0, 0.0) ); }
      pop_adi.set(i, i, complex<double>(1.0, 0.0) );

      // The following transformation is correct only for S_dia = 1
      U = (*basis_transform[traj]);// * (*proj_adi[traj]);
      
      pop_dia = U * pop_adi * U.H(); 
      
      vector<double> diag_els(ndia, 0.0);
      for(i=0;i<ndia; i++){ diag_els[i] = pop_dia.get(i,i).real(); }

      ksi = rnd.uniform(0.0, 1.0);
      act_states_dia[traj] = hop(diag_els, ksi); 
    }
  }

  // Obtain the adiabatic active states from the diabatic ones
  else{
    for(traj=0; traj<ntraj; traj++){
      i = act_states_dia[traj]; // active diabatic state
      pop_dia = *dm_dia[traj];
      for(j=0;j<nadi; j++){ pop_dia.set(j,j, complex<double>(0.0, 0.0) ); }
      pop_dia.set(i, i, complex<double>(1.0, 0.0) );

      // The following transformation is correct only for S_dia = 1
      U = (*basis_transform[traj]);// * (*proj_adi[traj]);
      
      pop_adi = U.H() * pop_dia * U; 
      
      vector<double> diag_els(nadi, 0.0);
      for(i=0;i<nadi; i++){ diag_els[i] = pop_adi.get(i,i).real(); }

      ksi = rnd.uniform(0.0, 1.0);
      act_states[traj] = hop(diag_els, ksi); 
    }
  }

}


/*
void dyn_variables::update_ampl_reproject(int option){

 changes all the adiabatic properties of the Hamiltonian by the matrix T

  option = 0 - use matrix U = T
         = 1 - use matrix U = T.H()



  CMATRIX U(*T); // option == 0
  if (option==1){  U = CMATRIX(T->H()); }
  else if(option==2 || option==4){

}

*/





CMATRIX orthogonalized_T(CMATRIX& T){
/**
   Return a unitary equivalent of the matrix T

  |psi_adi_tilde(t+dt)> = |psi_adi(t+dt)> T(t+dt)

  <psi_adi_tilde(t)|psi_adi_tilde(t+dt)> = <psi_adi(t)|psi_adi(t+dt)> T(t+dt)
        (should be close to I)                    (call it P (t, t+dt)  )

  So:

  I = P(t, t+dt) * T(t+dt)

  So, ideally T(t+dt) should be P^{-1}, but we also want it to be orthonormal

  T = P^{-1} * N

  We want:  T^+ T = N^+ (P^{-1})^+ P^(-1) N = I

  N =  (X^+ X)^-{1/2}, where X = P^(-1), indeed:

  N^+ (P^{-1})^+ P^(-1) N = ((X^+ X)^{-1/2})^+  (X^+ X) (X^+ X )^{-1/2} = 

  = [((X^+ X)^{-1/2})^+ ] [ (X^+ X )^{1/2}] = (X X^+)^{-1/2} (X^+ X)^{1/2} = 
 
  = (X^+)^{-1/2} X^{-1/2}   X^{1/2} (X^+)^{1/2} = (X^+)^{-1/2} (X^+)^{1/2} = I
  
  The input to this function should be P^{-1} 

*/

  int nst = T.n_cols;

  CMATRIX t2(nst,nst);
  CMATRIX t2_half(nst, nst);
  CMATRIX t2_i_half(nst, nst);
    
  t2 = T.H() * T;
  sqrt_matrix(t2, t2_half, t2_i_half);

  CMATRIX id(nst,nst); id.identity();
  if( std::abs( (t2_half * t2_i_half - id).max_elt()) > 1e-5 ){
    cout<<"Error in orthogonlized_T: S^{1/2} * S^{-1/2} is not an identity matrix\nExiting...\n";
    exit(0);
  }
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


vector<double> dyn_variables::compute_average_mash_pop(int rep){
/**
  Computing the MASH population estimators based on:

  (1) E. Runeson, J.; P. Fay, T.; E. Manolopoulos, D. Exciton Dynamics from the Mapping Approach to Surface Hopping: Comparison with Förster and Redfield Theories. Physical Chemistry Chemical Physics 2024, 26 (6), 4929–4938. https://doi.org/10.1039/D3CP05926J.

*/

  //======= First, compute the SE populations ========
  int sz,i;

  if(rep==0 || rep==2){ sz = ndia; }
  else if(rep==1 || rep==3){ sz = nadi; }

  MATRIX ave(sz, sz);
  ave = compute_average_dm(rep).real();

  vector<double> res(sz, 0.0);

  for(i=0; i<sz; i++){ res[i] = ave.get(i,i); }

  //====== Now transform them according to MASH estimators formula =====

  double sum = 0.0;
  for(i=1; i<=sz; i++){  sum += 1.0/i; }  sum -= 1.0;
  double alpha_N = (sz - 1.0)/sum;
  double beta_N = (1.0 - alpha_N)/sz;

  for(i=0; i<sz; i++){  res[i] = alpha_N * res[i] + beta_N; }

  return res;

}



vector<double> dyn_variables::compute_average_sh_pop(int rep){
/**
  Computing the SH population with the active state of the input rep.
  If rep is different from rep_sh, the resultant SH population is computed indirectly through the sampling based on diagonal elments of the density matrix
*/

  int sz, i, j, traj; 
  if(rep==0){ sz = ndia; }
  else if(rep==1){ sz = nadi; }
  else{ 
    cout<<"Can not compute SH population for representation = "<<rep<<"\nExiting...\n";
    exit(0);
  }

  vector<double> res(sz, 0.0);

  if(rep==0){
    vector<int> effective_states( act_states_dia );
    for(traj=0; traj<ntraj; traj++){ i = effective_states[traj]; res[i] += 1.0;  }
  }// rep == 0

  else if(rep==1){
    vector<int> effective_states( act_states );
    for(traj=0; traj<ntraj; traj++){ i = effective_states[traj]; res[i] += 1.0;  }
  }// rep == 1

  for(j=0; j<sz; j++){   res[j] = res[j] / (float)ntraj;   }

  return res;
}



vector<double> dyn_variables::compute_average_sh_pop_TR(int rep){
/**
  Computing the SH population based on the following work:
  Tempelaar, R.; Reichman, D. R. Generalization of Fewest-Switches Surface Hopping for Coherences. 
  The Journal of Chemical Physics 2018, 148 (10), 102309. https://doi.org/10.1063/1.5000843
*/
  int sz, i, j, traj; 
  if(rep==0){ sz = ndia; }
  else if(rep==1){ sz = nadi; }
  else{ 
    cout<<"Can not compute SH population for representation = "<<rep<<"\nExiting...\n";
    exit(0);
  }

  vector<double> res(sz, 0.0);

  if(rep==0){
    vector<int> effective_states( act_states );

    CMATRIX pop_adi(nadi, nadi);
    CMATRIX pop_dia(ndia, ndia);
    CMATRIX U(ndia, nadi);

    for(traj=0; traj<ntraj; traj++){
      i = effective_states[traj]; // active adiabatic state
      pop_adi = *dm_adi[traj];
      for(j=0;j<nadi; j++){ pop_adi.set(j,j, complex<double>(0.0, 0.0) ); }
      pop_adi.set(i, i, complex<double>(1.0, 0.0) );

      // The following transformation is correct only for S_dia = 1
      U = (*basis_transform[traj]);// * (*proj_adi[traj]);
      
      pop_dia = U * pop_adi * U.H(); 
      for(j=0; j<ndia; j++){ res[j] += pop_dia.get(j,j).real(); }
    }
  }// rep == 0

  else if(rep==1){
    vector<int> effective_states( act_states_dia );

    CMATRIX pop_adi(nadi, nadi);
    CMATRIX pop_dia(ndia, ndia);
    CMATRIX U(ndia, nadi);

    for(traj=0; traj<ntraj; traj++){
      i = effective_states[traj]; // active diabatic state
      pop_dia = *dm_dia[traj];
      for(j=0;j<ndia; j++){ pop_dia.set(j,j, complex<double>(0.0, 0.0) ); }
      pop_dia.set(i, i, complex<double>(1.0, 0.0) );

      // The following transformation is correct only for S_dia = 1
      U = (*basis_transform[traj]);// * (*proj_adi[traj]);
      
      pop_adi = U.H() * pop_dia * U; 
      for(j=0; j<nadi; j++){ res[j] += pop_adi.get(j,j).real(); }
    }
  }// rep == 1

  for(j=0; j<sz; j++){   res[j] = res[j] / (float)ntraj;   }

  return res;
}


MATRIX dyn_variables::compute_coherence_indicator(int rep){

  int sz;
  if(rep==0 || rep==2){ sz = ndia; }
  else if(rep==1 || rep==3){ sz = nadi; }

  MATRIX res(sz, sz);
  
  CMATRIX dm(sz, sz); MATRIX temp(sz, ntraj);

  // Making a temporary matrix collecting trajectory-wise population elements
  for(int traj=0; traj<ntraj; traj++){
    if(rep==0 || rep==2){   dm = *dm_dia[traj]; }
    else if(rep==1 || rep==3){ dm = *dm_adi[traj]; }
  
    for(int i=0; i<sz; i++){temp.set(i, traj, dm.get(i,i).real() );}
  }

  res = temp * temp.T() / (float)ntraj;

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

