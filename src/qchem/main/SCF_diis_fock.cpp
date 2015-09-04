#include "SCF.h"

namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_qm{


double scf_diis_fock(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM){


/// This function implements the SCF based on the DIIS with Fock matrix mixing(extrapolation)
/// if prms.use_diis = 0 => routine reduces to standard SCF
/// 
/// I tried following the recipe from:
/// [1] Pulay, P. "Convergence Acceleration of Iterative Sequences. The Case of SCF Iteration" Chem. Phys. Lett. 1980, 73, 393
/// [2] Pulay, P. "Improved SCF Convergence Acceleration" J. Comp. Chem. 1982, 3, 556


  int i;
  double lamb_min;
  std::string eigen_method="generalized";


  //----------- Control parameters ---------
  int iter = 0;
  int Niter = prms.Niter;

  int Norb = el->Norb;
  int Nocc_alp = el->Nocc_alp;
  int Nocc_bet = el->Nocc_bet;

  double den_tol = prms.den_tol;  
  double den_err = 2.0*den_tol;

  double ene_tol = prms.etol;
  double Eelec_prev = 0.0;
  double Eelec = 0.0;
  double dE = 2.0*ene_tol;

  vector<Timer> bench_t2(4);


  Electronic* el_tmp;   el_tmp = new Electronic(el);


  // DIIS auxiliary variables
  int N_diis = 0;
  int stop_diis = 0;
  int start_diis = 0;
  int diis_delay = 0;

  DIIS* diis_alp;
  DIIS* diis_bet;
  MATRIX* err_alp;
  MATRIX* err_bet;

  err_alp = new MATRIX(Norb,Norb);
  err_bet = new MATRIX(Norb,Norb);  

  if(prms.use_diis){

    diis_alp = new DIIS(prms.diis_max,Norb);
    diis_bet = new DIIS(prms.diis_max,Norb);

  }// use_diis



  // Initialization:
  *el_tmp->P     = *el->P;
  *el_tmp->P_alp = *el->P_alp;
  *el_tmp->P_bet = *el->P_bet;
  *el_tmp->C_alp = *el->C_alp;
  *el_tmp->C_bet = *el->C_bet;


  // Initialization:
  // F_0 = F(D_0), F~_0 = F(D~_0) = F_0
  Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix

  Eelec_prev = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao,el_tmp->Fao_alp);
  Eelec_prev+= ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao,el_tmp->Fao_bet);


  
  //=========================== Now enter main SCF cycle ===========================================
  ofstream f1("energy.txt",ios::out);

  cout<<"----------------------- Entering main SCF cycle for RHF calculations --------------------\n"; 

  do{

    cout<<"===============Iteration# "<<iter<<" =====================\n";   

    if(iter>prms.diis_start_iter){ start_diis = 1; }


    // Step 1: Assemble a new Fock matrix at a given density matrix F_{k} = F(D_{k})
    Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el_tmp - now contains an updated Fock matrix


    // Step 2: Compute error matrices
    // Here one should use the same density that was used to construct Fock matix, not the density produced from the Fock matrix
    *err_alp = (*el_tmp->Fao_alp) * (*el_tmp->P_alp) * (*el_tmp->Sao) - (*el_tmp->Sao) * (*el_tmp->P_alp) * (*el_tmp->Fao_alp);
    *err_bet = (*el_tmp->Fao_bet) * (*el_tmp->P_bet) * (*el_tmp->Sao) - (*el_tmp->Sao) * (*el_tmp->P_bet) * (*el_tmp->Fao_bet);

    // Transform errors into MO basis
    // Note: as indicated in [Hamilton,Pulay, JCP 84, 5728 (1986) ], the error matrices should be expressed in same basis
    // either AO, or !!fixed!! MO basis - e.g. from first iteraction
    // Thus, we use AO basis from now on, the commented lines below show how to convert to MO basis, but this is not good
    // because MO basis is not fixed
    // In fact, we need to transform to MO basis to get OE - EO matrix, otherwise a converged result will always have non-zero error in density matrix
    *err_alp = ((*el_tmp->C_alp).T()) * *err_alp * (*el_tmp->C_alp);
    *err_bet = ((*el_tmp->C_bet).T()) * *err_bet * (*el_tmp->C_bet);

    // Scale by a common factor:
//    double alp_err = ((*err_alp).T() * (*err_alp)).tr();
//    double bet_err = ((*err_bet).T() * (*err_bet)).tr();
//    double tot_err = 0.5*(alp_err + bet_err);
//    if(alp_err>0.0) { *err_alp *= (tot_err/alp_err); }
//    if(bet_err>0.0) { *err_bet *= (tot_err/bet_err); }



    // Error in commutation:  max_elt |FPS - SPF|_{alp} + max_elt |FPS - SPF|_{bet}
    den_err = 0.5*(fabs(err_alp->max_elt()) +  fabs(err_bet->max_elt()));   
    cout<<"commutation error = "<<den_err<<endl;


    // Step 3: Either terminate or continue with the search
    if(den_err<den_tol && fabs(dE)<ene_tol){  ;;  }  
    else{

      if(prms.use_diis && start_diis){

        // DIIS Step 1: Store error matrix and current Fock matrices
        diis_alp->add_diis_matrices(el_tmp->Fao_alp, err_alp);
        diis_bet->add_diis_matrices(el_tmp->Fao_bet, err_bet);  


        if(diis_delay >= diis_alp->N_diis_max){

          // DIIS Setp 2: Compute extrapolated Fock matrix and replace the current Fock matrix by the extrapolated one
          diis_alp->extrapolate_matrix(el_tmp->Fao_alp);
          diis_bet->extrapolate_matrix(el_tmp->Fao_bet);

        }

        diis_delay++;

      }// use_diis && start_diis

      //---------- Obtain a new density for this iteration -------------------------
      // Step 4: Diagonalize F_k, assemble corresponding new density - D_{k+1} via aufbau or Fermi populations
      Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, prms.pop_opt, el_tmp->Fao_alp, el_tmp->Sao, el_tmp->C_alp, el_tmp->E_alp, el_tmp->bands_alp, el_tmp->occ_alp, el_tmp->P_alp, bench_t2);
      Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, prms.pop_opt, el_tmp->Fao_bet, el_tmp->Sao, el_tmp->C_bet, el_tmp->E_bet, el_tmp->bands_bet, el_tmp->occ_bet, el_tmp->P_bet, bench_t2);


      *el_tmp->P = *el_tmp->P_alp + *el_tmp->P_bet;

      cout<<"at this point max|FC-SCE|= "<<(*el_tmp->Fao_alp * *el_tmp->C_alp - *el_tmp->Sao * *el_tmp->C_alp * *el_tmp->E_alp).max_elt()<<endl;

      cout<<"|F|= "<<el_tmp->Fao_alp->max_elt()<<endl;
      cout<<"|S|= "<<el_tmp->Sao->max_elt()<<endl;
      cout<<"|C|= "<<el_tmp->C_alp->max_elt()<<endl;
      cout<<"|E|= "<<el_tmp->E_alp->max_elt()<<endl;


      
    }// else: den_err>=den_tol

   
    Eelec = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao,el_tmp->Fao_alp);
    Eelec+= ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao,el_tmp->Fao_bet);

    dE = Eelec - Eelec_prev;
    Eelec_prev = Eelec;


    show_bands(el_tmp->Norb, el_tmp->Nocc_alp, el_tmp->bands_alp, el_tmp->occ_alp);
    show_bands(el_tmp->Norb, el_tmp->Nocc_bet, el_tmp->bands_bet, el_tmp->occ_bet);
        
    cout<<"Total electronic energy = "<<Eelec<<endl;
    cout<<"Mulliken scf orbital populations:\n";
    for(i=0;i<el_tmp->Norb;i++){   cout<<"i= "<<i<<"  pop(gross)= "<<el_tmp->Mull_orb_pop_gross[i]<<"  pop(net)= "<<el_tmp->Mull_orb_pop_net[i]<<endl;  }
    cout<<"Mulliken scf charges:\n";
    for(i=0;i<mol.Nnucl;i++){   cout<<"i= "<<i<<"  q(gross)= "<<mol.Mull_charges_gross[i]<<"  q(net)= "<<mol.Mull_charges_net[i]<<endl;  }

    if(0){
      cout<<"Hao=\n"<<*el_tmp->Hao<<endl;
      cout<<"P=\n"<<*el_tmp->P<<endl;
      cout<<"F_alp=\n"<<*el_tmp->Fao_alp<<endl;
    }

    f1 << iter<<" Eelec= "<<Eelec<<" dE= "<<dE<<" den_err = "<<den_err<<endl; //<<" Tr(S*D)= "<<(*el_tmp->Sao * *el_tmp->P).tr()<<endl;

    iter++;    

 
  }while(iter<Niter && (den_err>den_tol || fabs(dE)>ene_tol) );

  f1.close();


  if(fabs(den_err)>den_tol){
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    cout<<"!!!!! Error: Convergence in density is not achieved after "<<Niter<<" iterations\n den_err = "<<den_err<<" !!!!!\n";
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    exit(0);
  }

  if(fabs(dE)>ene_tol){
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    cout<<"!!!!! Error: Convergence in energy is not achieved after "<<Niter<<" iterations\n dE = "<<dE<<" !!!!!\n";
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    exit(0);
  }

  if(prms.do_annihilate==1){ annihilate(Nocc_alp,Nocc_bet,el_tmp->P_alp,el_tmp->P_bet); }

  *el->P_alp = *el_tmp->P_alp;
  *el->P_bet = *el_tmp->P_bet;
  *el->P     = *el->P_alp + *el->P_bet;

  // Compute Fock matrix at final density
  Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el, el0,mem); // el - now contains an updated Fock matrix

  // Update eigenvalues and eigenvectors of final Fock matrix, but do not modify the density matrix:
  Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, 0, el->Fao_alp, el->Sao, el->C_alp, el->E_alp, el->bands_alp, el->occ_alp, el_tmp->P_alp, bench_t2);  
  Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, 0, el->Fao_bet, el->Sao, el->C_bet, el->E_bet, el->bands_bet, el->occ_bet, el_tmp->P_bet, bench_t2);  


  // Recompute energy using converged density and corresponding Fock matrix
  Eelec = ::energy_elec(Norb,el->P_alp,el->Hao,el->Fao_alp) + ::energy_elec(Norb,el->P_bet,el->Hao,el->Fao_bet);


  // Clean up the memory
  el_tmp->~Electronic();
  delete err_alp;
  delete err_bet;


  return Eelec;
}


}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian



