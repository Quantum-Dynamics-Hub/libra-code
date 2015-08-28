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

/****************************************************************************//**
 \file SCF.cpp
 \brief Implementation of the self-consistent field methods

 Details:    

  This file contains the following functions:

  double scf_oda(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                 vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                 Electronic* el,Electronic* el0, Memory*)

  double scf_diis_fock(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                       vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                       Electronic* el,Electronic* el0, Memory*)

  double scf_diis_dm(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                     vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                     Electronic* el,Electronic* el0, Memory*)




****************************************************************************/
#include "DIIS.h"
#include "SCF.h"
#include "Engine.h"
#include "Hamiltonian.h"
#include "units.h"
#include "Timer.h"


double scf(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
           vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
           Electronic* el,Electronic* el0, Memory* mem){

  double res = 0.0;

  if(prms.scf_algo=="oda"){
    res = scf_oda(prms,modprms,mol,fragment,basis_fo,basis_ao,at_orbitals,el,el0,mem);
//    res = scf_oda_disk(prms,modprms,mol,fragment,basis_fo,basis_ao,at_orbitals,el,el0,mem);
  }
//  else if(prms.scf_algo=="oda_disk"){
//    res = scf_oda(prms,modprms,mol,fragment,basis_fo,basis_ao,at_orbitals,el,el0,mem);
//    res = scf_oda_disk(prms,modprms,mol,fragment,basis_fo,basis_ao,at_orbitals,el,el0,mem);
//  }
  else if(prms.scf_algo=="diis_fock"){
    res = scf_diis_fock(prms,modprms,mol,fragment,basis_fo,basis_ao,at_orbitals,el,el0,mem);
  }
  else if(prms.scf_algo=="diis_dm"){
    res = scf_diis_dm(prms,modprms,mol,fragment,basis_fo,basis_ao,at_orbitals,el,el0,mem);
  }

  return res;

}


double scf_oda(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
               vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
               Electronic* el,Electronic* el0, Memory* mem){

/// This function implements the SCF based on the optimal damping algorithm (ODA)
/// which uses fractional occupation numbers, leading to robust convergence in difficult cases
/// See more details in: 
/// [1] Kudin K.N.; Scuseria, G.E.; Cances, E. J. Chem. Phys. 116, 8255 (2002)
/// [2] Cances J. Chem. Phys. 114, 10616 (2001) 

  int BM = 1; // bemchmark flag

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

  vector<Timer> bench_t(10); // timers for different type of operations
  vector<Timer> bench_t2(4);


  if(BM){ bench_t[5].start(); }
  MATRIX* dP;           dP          = new MATRIX(Norb,Norb);  // dP = P_k+1 - P_k
  MATRIX* temp;         temp        = new MATRIX(Norb,Norb);
  MATRIX* P_old_alp;    P_old_alp   = new MATRIX(Norb,Norb);
  MATRIX* P_old_bet;    P_old_bet   = new MATRIX(Norb,Norb);
  MATRIX* P_old;        P_old       = new MATRIX(Norb,Norb);

  MATRIX* P_alp;        P_alp       = new MATRIX(Norb,Norb);
  MATRIX* P_bet;        P_bet       = new MATRIX(Norb,Norb);
  MATRIX* P;            P           = new MATRIX(Norb,Norb);

  MATRIX* P_til_alp;    P_til_alp   = new MATRIX(Norb,Norb);
  MATRIX* P_til_bet;    P_til_bet   = new MATRIX(Norb,Norb);
  MATRIX* P_til;        P_til       = new MATRIX(Norb,Norb);

  MATRIX* Fao_alp;      Fao_alp     = new MATRIX(Norb,Norb);
  MATRIX* Fao_bet;      Fao_bet     = new MATRIX(Norb,Norb);


//  MATRIX* dFao_alp_dP_alp;  dFao_alp_dP_alp = new MATRIX(Norb,Norb); *dFao_alp_dP_alp = 0.0;
//  MATRIX* dFao_alp_dP_bet;  dFao_alp_dP_bet = new MATRIX(Norb,Norb); *dFao_alp_dP_bet = 0.0;
//  MATRIX* dFao_bet_dP_alp;  dFao_bet_dP_alp = new MATRIX(Norb,Norb); *dFao_bet_dP_alp = 0.0;
//  MATRIX* dFao_bet_dP_bet;  dFao_bet_dP_bet = new MATRIX(Norb,Norb); *dFao_bet_dP_bet = 0.0;



  MATRIX* Fao_til_alp;  Fao_til_alp = new MATRIX(Norb,Norb);
  MATRIX* Fao_til_bet;  Fao_til_bet = new MATRIX(Norb,Norb);


  Electronic* el_tmp;   el_tmp = new Electronic(el);

  if(BM){ bench_t[5].stop(); }



  // Interface
  *P = *el->P;
  *P_alp = *el->P_alp;
  *P_bet = *el->P_bet;
  *Fao_alp = *el->Fao_alp;
  *Fao_bet = *el->Fao_bet;
  
  // Old
  *P_old_alp = *P_alp;
  *P_old_bet = *P_bet;
  *P_old = *P;

  // Tilda
  // D~_0 = D_0
  *P_til_alp = *P_alp;
  *P_til_bet = *P_bet;
  *P_til = *P;
  

  // Initialization:
  // F_0 = F(D_0), F~_0 = F(D~_0) = F_0
  if(BM){ bench_t[1].start(); }
  Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
  if(BM){ bench_t[1].stop(); }
  *Fao_til_alp = *el_tmp->Fao_alp;
  *Fao_til_bet = *el_tmp->Fao_bet;


  if(BM){ bench_t[0].start(); }
  Eelec_prev = ::energy_elec(Norb,el->P_alp,el->P_bet, el->Hao, el->Hao, el->Fao_alp, el->Fao_bet,
                 el->dFao_alp_dP_alp,el->dFao_alp_dP_bet,el->dFao_bet_dP_alp,el->dFao_bet_dP_bet,
                 temp);
  if(BM){ bench_t[0].stop(); }


// Debug:
//  cout<<"In SCF: initital density matrix: P = \n"<<*P<<endl; 

  //=========================== Now enter main SCF cycle ===========================================
  ofstream f1("energy.txt",ios::out);

  cout<<"----------------------- Entering main SCF cycle for RHF calculations --------------------\n"; 

  do{

    cout<<"===============Iteration# "<<iter<<" =====================\n";   

    //---------- Obtain a new density for this iteration -------------------------
    // ODA Step 1: Diagonalize F~_k, assemble D_{k+1} via aufbau (so forcibly set prms.pop_opt = 0) 
    if(BM){ bench_t[2].start(); }
    if(prms.use_damping==0){
      Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, 0, Fao_til_alp, el_tmp->Sao, el_tmp->C_alp, el_tmp->E_alp, el_tmp->bands_alp, el_tmp->occ_alp, P_alp, bench_t2);
      Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, 0, Fao_til_bet, el_tmp->Sao, el_tmp->C_bet, el_tmp->E_bet, el_tmp->bands_bet, el_tmp->occ_bet, P_bet, bench_t2);
      *P = *P_alp + *P_bet;
    }
    else if(prms.use_damping==1){
      // This is only specific to SC-EHT, since dF/dP doesn't depend on P, so we can use it computed at any stage
      // dFao_alp_dP_alp

//  This is original!!!
      Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, prms.pop_opt, Fao_til_alp, el_tmp->Sao, el_tmp->C_alp, el_tmp->E_alp, el_tmp->bands_alp, el_tmp->occ_alp, P_alp, bench_t2);
      Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, prms.pop_opt, Fao_til_bet, el_tmp->Sao, el_tmp->C_bet, el_tmp->E_bet, el_tmp->bands_bet, el_tmp->occ_bet, P_bet, bench_t2);

      // This is corrected
//      *temp->Fao_alp = *Fao_til_alp + P_til_alp * el_tmp->dFao_alp_dP_alp + P_til_bet * el_tmp->dFao_alp_dP_bet;
//      *temp->Fao_bet = *Fao_til_bet + P_til_alp * el_tmp->dFao_bet_dP_alp + P_til_bet * el_tmp->dFao_bet_dP_bet;
//      Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, prms.pop_opt, temp->Fao_alp, el_tmp->Sao, el_tmp->C_alp, el_tmp->E_alp, el_tmp->bands_alp, el_tmp->occ_alp, P_alp, bench_t2);
//      Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, prms.pop_opt, temp->Fao_bet, el_tmp->Sao, el_tmp->C_bet, el_tmp->E_bet, el_tmp->bands_bet, el_tmp->occ_bet, P_bet, bench_t2);

      *P = *P_alp + *P_bet;
    }
    if(BM){ bench_t[2].stop(); }


// Debug:
//    cout<<"D_{k+1} = D(F~_k) = \n"<<*P<<endl;
//    cout<<"D~_{k} = \n"<<*P_til<<endl;

    if(BM){ bench_t[3].start(); }
      cout<<"Pmax = "<<P->max_elt()<<endl;    
      cout<<"Pold_max = "<<P_old->max_elt()<<endl; 

//      den_err = fabs((*P - *P_old).max_elt());    //  !!!
      den_err = fabs((*P_til - *P_old).max_elt());      
      cout<<"den_err = "<<den_err<<endl;
    if(BM){ bench_t[3].stop(); }


    if(BM){ bench_t[3].start(); }
//      *P_old_alp = *P_alp;
//      *P_old_bet = *P_bet;
//      *P_old = *P;
      *P_old_alp = *P_til_alp;
      *P_old_bet = *P_til_bet;
      *P_old = *P_til;

    if(BM){ bench_t[3].stop(); }


    //--------- ODA Step 2: Either terminate or continue with the search ------------
    // D_{k+1} - D~_k
    if(den_err<den_tol && fabs(dE)<ene_tol){  ;;  }  
    else{

      //------- ODA Step 3: Assemble F_{k+1} = F(D_{k+1}) ---------
      if(BM){ bench_t[3].start(); }
        *dP = *P - *P_til;
        *el_tmp->P_alp = *P_alp;
        *el_tmp->P_bet = *P_bet;
        *el_tmp->P = *P;
      if(BM){ bench_t[3].stop(); }


      if(BM){ bench_t[1].start(); }
        Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
      if(BM){ bench_t[1].stop(); }


      if(BM){ bench_t[3].start(); }
        *Fao_alp = *el_tmp->Fao_alp;   
        *Fao_bet = *el_tmp->Fao_bet;
      if(BM){ bench_t[3].stop(); }

  
      //-----------  ODA Step 4: Solve the line search problem (via interpolation) or use fixed step -------------
      lamb_min = 0.0;
      if(prms.use_damping){

        if(iter<=prms.damping_start){  lamb_min = 1.0; }
        else{ lamb_min = prms.damping_const;   }
        cout<<"Using constant lamb_min = "<<lamb_min<<endl;

      }else{

        cout<<"Line search:\n";   

        //-------------- lamb = 1 ------------------
        if(BM){ bench_t[3].start(); }
          double lamb = 1.0;        
          *el_tmp->P     = *P_til     + lamb * *dP;
          *el_tmp->P_alp = *P_til_alp + lamb * (0.5* *dP);
          *el_tmp->P_bet = *P_til_bet + lamb * (0.5* *dP);
        if(BM){ bench_t[3].stop(); }

        if(BM){ bench_t[1].start(); }
          Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
        if(BM){ bench_t[1].stop(); }

        if(BM){ bench_t[0].start();}
          double en1 = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->P_bet, el_tmp->Hao, el_tmp->Hao, el_tmp->Fao_alp, el_tmp->Fao_bet,
                         el_tmp->dFao_alp_dP_alp,el_tmp->dFao_alp_dP_bet,el_tmp->dFao_bet_dP_alp,el_tmp->dFao_bet_dP_bet,
                         temp);
        if(BM){ bench_t[0].stop();}


        //-------------- lamb = 0.5 ------------------
        if(BM){ bench_t[3].start(); }
          lamb = 0.5;        
          *el_tmp->P     = *P_til     + lamb * *dP;
          *el_tmp->P_alp = *P_til_alp + lamb * (0.5* *dP);
          *el_tmp->P_bet = *P_til_bet + lamb * (0.5* *dP);
        if(BM){ bench_t[3].stop(); }

        if(BM){ bench_t[1].start(); }
          Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
        if(BM){ bench_t[1].stop(); }

        if(BM){ bench_t[0].start(); }
          double en2 = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->P_bet, el_tmp->Hao, el_tmp->Hao, el_tmp->Fao_alp, el_tmp->Fao_bet,
                         el_tmp->dFao_alp_dP_alp,el_tmp->dFao_alp_dP_bet,el_tmp->dFao_bet_dP_alp,el_tmp->dFao_bet_dP_bet,
                         temp);
        if(BM){ bench_t[0].stop(); }


        //-------------- lamb = 0.0 ------------------
        if(BM){ bench_t[3].start(); }
          lamb = 0.0;        
          *el_tmp->P     = *P_til     + lamb * *dP;
          *el_tmp->P_alp = *P_til_alp + lamb * (0.5* *dP);
          *el_tmp->P_bet = *P_til_bet + lamb * (0.5* *dP);
        if(BM){ bench_t[3].stop(); }

        if(BM){ bench_t[1].start(); }
          Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
        if(BM){ bench_t[1].stop(); }

        if(BM){ bench_t[0].start(); }
          double en0 = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->P_bet, el_tmp->Hao, el_tmp->Hao, el_tmp->Fao_alp, el_tmp->Fao_bet,
                         el_tmp->dFao_alp_dP_alp,el_tmp->dFao_alp_dP_bet,el_tmp->dFao_bet_dP_alp,el_tmp->dFao_bet_dP_bet,
                         temp);
        if(BM){ bench_t[0].stop(); }

     
        // After this operation all Fock matrices (in el_tmp!) are as if no intermediate calculations were performed
        // for Fao_* is actually F(P_*_tmp) = F~_k
        // At this point Fao_ contains F~_k
        cout<<"E(0)= "<<en0<<endl;
        cout<<"E(1/2)= "<<en2<<endl;
        cout<<"E(1)= "<<en1<<endl;

        double _c = en0;
        double _b = 4.0*en2 - en1 - 3.0*en0;
        double _a = en1 - en0 - _b;

        cout<<"Interpolation polynomial = "<<_a<<" * lamb^2 + "<<_b<<" * lamb + "<<_c<<endl;        
        lamb_min = 0.0; 
        if(fabs(_a)>1e-10){ lamb_min = -_b/(2.0*_a); }

        cout<<"lamb_min = "<<lamb_min<<endl;

        if(0<lamb_min && lamb_min<1){
          Eelec = _a*lamb_min*lamb_min + _b*lamb_min + _c;
          cout<<"Functional minimum = "<<Eelec<<endl;
        }
        else{
          lamb_min = (en0<en1)?0.0:1.0;
          Eelec = ((en0<en1)?en0:en1);
          cout<<"infinum at = "<<((en0<en1)?"lamb_min = 0.0":"lamb_min = 1.0")<<" value = "<<Eelec<<endl;
        }

      }// do not use damping

      // ODA Step 5:
      // P~{k+1} = P~{k} + lamb_min * dP = (1 - lamb_min)*P~_k + lamb_min * P_{k+1}
      // dP = P(k+1) - P~(k)
      if(BM){ bench_t[3].start(); }
        *P_til_alp   = (1.0 - lamb_min) * (*P_til_alp) + lamb_min * (*P_alp);
        *P_til_bet   = (1.0 - lamb_min) * (*P_til_bet) + lamb_min * (*P_bet);    
        *P_til       = *P_til_alp + *P_til_bet;

        // F~{k+1} = (1 - lamb_min)*F~_k + lamb_min * F_{k+1}   - this is original approach, but less general
        // *Fao_til_alp = (1.0 - lamb_min) * (*Fao_til_alp) + lamb_min * (*Fao_alp);
        // *Fao_til_bet = (1.0 - lamb_min) * (*Fao_til_bet) + lamb_min * (*Fao_bet);
  
        // in fact, F~{k+1} = F(D~_{k+1})  - this is more general approach
        *el_tmp->P     = *P_til;    
        *el_tmp->P_alp = *P_til_alp;
        *el_tmp->P_bet = *P_til_bet;
      if(BM){ bench_t[3].stop(); }


      if(BM){ bench_t[1].start(); }
        Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
      if(BM){ bench_t[1].stop(); }

        *Fao_til_alp = *el_tmp->Fao_alp;
        *Fao_til_bet = *el_tmp->Fao_bet;

      
    }// else: den_err>=den_tol

   
    //------ Recompute current energy using extrapolated (or old) density matrix ---------------
    if(BM){ bench_t[3].start(); }
      *el_tmp->P     = *P_til;//     + lamb_min * *dP;
      *el_tmp->P_alp = *P_til_alp;// + lamb_min * (0.5* *dP);
      *el_tmp->P_bet = *P_til_bet;// + lamb_min * (0.5* *dP);
    if(BM){ bench_t[3].stop(); }

    if(BM){ bench_t[1].start(); }
      Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
    if(BM){ bench_t[1].stop(); }

    if(BM){ bench_t[0].start(); }
      Eelec = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->P_bet, el_tmp->Hao, el_tmp->Hao, el_tmp->Fao_alp, el_tmp->Fao_bet,
                el_tmp->dFao_alp_dP_alp,el_tmp->dFao_alp_dP_bet,el_tmp->dFao_bet_dP_alp,el_tmp->dFao_bet_dP_bet,
                temp);
    if(BM){ bench_t[0].stop(); }


    //------- Compute energy change --------------------
    dE = Eelec - Eelec_prev;
    Eelec_prev = Eelec;

     
    //--------------- Some output ---------------------
    if(BM){ bench_t[4].start(); }    
      if(0){  // Suppress excessive output - it takes time

        show_bands(el_tmp->Norb, el_tmp->Nocc_alp, el_tmp->bands_alp, el_tmp->occ_alp);
        show_bands(el_tmp->Norb, el_tmp->Nocc_bet, el_tmp->bands_bet, el_tmp->occ_bet);
        
        cout<<"Total electronic energy = "<<Eelec<<endl;
        cout<<"Mulliken scf orbital populations:\n";
        for(i=0;i<el_tmp->Norb;i++){   cout<<"i= "<<i<<"  pop(gross)= "<<el_tmp->Mull_orb_pop_gross[i]<<"  pop(net)= "<<el_tmp->Mull_orb_pop_net[i]<<endl;  }
        cout<<"Mulliken scf charges:\n";
        for(i=0;i<mol.Nnucl;i++){   cout<<"i= "<<i<<"  q(gross)= "<<mol.Mull_charges_gross[i]<<"  q(net)= "<<mol.Mull_charges_net[i]<<endl;  }
      }

      f1 << iter<<" Eelec= "<<Eelec<<" dE= "<<dE<<" den_err = "<<den_err<<endl;
      //<<" Tr(S*D)= "<<(*el_tmp->Sao * *P_til).tr()<<endl;
    if(BM){ bench_t[4].stop(); }    


    //------------- Continue iterative process -------------
    iter++;    

 
  }while(iter<Niter && (den_err>den_tol || fabs(dE)>ene_tol) );

  f1.close();




  if(prms.do_annihilate==1){ annihilate(Nocc_alp,Nocc_bet,P_til_alp,P_til_bet); }

  if(BM){ bench_t[3].start(); }
  *el->P_alp = *P_til_alp;
  *el->P_bet = *P_til_bet;
  *el->P     = *el->P_alp + *el->P_bet;

  *el->Fao_alp = *Fao_til_alp;
  *el->Fao_bet = *Fao_til_bet;
  if(BM){ bench_t[3].stop(); }


  if(BM){ bench_t[1].start(); }
  Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el, el0,mem); // el - now contains an updated Fock matrix
  if(BM){ bench_t[1].stop(); }

  // Test: At this point F(P~) = F~, so it is valid to use P~ to construct F~ via normal rules
  //cout<<"Difference of Fock matrices: \n"<<*el->Fao_alp - *Fao_til_alp<<endl;

  // Update eigenvalues and eigenvectors of final Fock matrix, but do not modify the density matrix:
  if(BM){ bench_t[2].start(); }
  Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, 0, el->Fao_alp, el->Sao, el->C_alp, el->E_alp, el->bands_alp, el->occ_alp, P_alp, bench_t2);  
  Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, 0, el->Fao_bet, el->Sao, el->C_bet, el->E_bet, el->bands_bet, el->occ_bet, P_bet, bench_t2);  
  if(BM){ bench_t[2].stop(); }

  bench_t[0].start();

// Old!!!
//  Eelec = ::energy_elec(Norb,el->P_alp,el->Hao,el->Fao_alp) + ::energy_elec(Norb,el->P_bet,el->Hao,el->Fao_bet);

// New!!!
/*
    *temp = *el->Fao_alp;
    if(el->dFao_alp_dP_alp->max_elt()>1e-10){   *temp += *el->P_alp * *el->dFao_alp_dP_alp;    }
    if(el->dFao_alp_dP_bet->max_elt()>1e-10){   *temp += *el->P_bet * *el->dFao_alp_dP_bet;    }
    Eelec = ::energy_elec(Norb,el->P_alp,el->Hao, temp);


    *temp = *el->Fao_bet;
    if(el->dFao_alp_dP_alp->max_elt()>1e-10){   *temp += *el->P_alp * *el->dFao_alp_dP_alp;    }
    if(el->dFao_alp_dP_bet->max_elt()>1e-10){   *temp += *el->P_bet * *el->dFao_alp_dP_bet;    }

//    *temp = (*el->Fao_bet + *el->P_alp * *el->dFao_bet_dP_alp + *el->P_bet * *el->dFao_bet_dP_bet); - less efficient
    Eelec+= ::energy_elec(Norb,el->P_bet,el->Hao, temp);
*/
    Eelec = ::energy_elec(Norb,el->P_alp,el->P_bet, el->Hao, el->Hao, el->Fao_alp, el->Fao_bet,
              el->dFao_alp_dP_alp,el->dFao_alp_dP_bet,el->dFao_bet_dP_alp,el->dFao_bet_dP_bet,
              temp);



  bench_t[0].stop();

/*  DEBUG
  cout<<"In SCF:\n";
  cout<<"F*C = \n"<<*el->Fao_alp * *el->C_alp<<endl;
  cout<<"S*C*E = \n"<<*el->Sao * *el->C_alp * *el->E_alp<<endl;
  cout<<"MO orthogonality:\n";

  // MO
  for(i=0;i<el->Norb;i++){
    for(int j=0;j<el->Norb;j++){

      double res = 0.0;

      int ii = el->bands_alp[i].first;
      int jj = el->bands_alp[j].first;
       

      for(int a=0;a<el->Norb;a++){    // over basis of fragment A at time t
        for(int b=0;b<el->Norb;b++){  // over basis of fragment B at time t

          res += el->C_alp->M[a*el->Norb+ii] * el->C_alp->M[b*el->Norb+jj] * el->Sao->M[a*el->Norb+b];
                
        }// for b
      }// for a
      cout<<res<<"  ";
    }
    cout<<endl;
  }
*/


  // Clean up the memory
  if(BM){ bench_t[5].start(); }
  delete dP;
  delete temp;
  delete P_old_alp;  
  delete P_old_bet;  
  delete P_old;      

  delete P_alp;     
  delete P_bet;     
  delete P;         

  delete P_til_alp;   
  delete P_til_bet;   
  delete P_til;       

  delete Fao_alp;     
  delete Fao_bet;     

  delete Fao_til_alp; 
  delete Fao_til_bet; 

  el_tmp->~Electronic();
  if(BM){ bench_t[5].stop(); }


  if(BM){
    cout<<"Time for energy calculation = "<<bench_t[0].show()<<endl;
    cout<<"Time for Fock matrix formaion = "<<bench_t[1].show()<<endl;
    cout<<"Time for Fock diagonalization and density matrix formation = "<<bench_t[2].show()<<endl;
    cout<<"   - eigensolver    = "<<bench_t2[0].show()<<endl;
    cout<<"   - sorting        = "<<bench_t2[1].show()<<endl;
    cout<<"   - populate       = "<<bench_t2[2].show()<<endl;
    cout<<"   - density matrix = "<<bench_t2[3].show()<<endl;
    cout<<"Time for matrix operations = "<<bench_t[3].show()<<endl;
    cout<<"Time for output = "<<bench_t[4].show()<<endl;
    cout<<"Time for allocation/deallocation = "<<bench_t[5].show()<<endl;
 
   

  }//

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



  return Eelec;

}


double scf_oda_disk(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
               vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
               Electronic* el,Electronic* el0, Memory* mem){

/// This function implements the SCF based on the optimal damping algorithm (ODA)
/// which uses fractional occupation numbers, leading to robust convergence in difficult cases
/// See more details in: 
/// [1] Kudin K.N.; Scuseria, G.E.; Cances, E. J. Chem. Phys. 116, 8255 (2002)
/// [2] Cances J. Chem. Phys. 114, 10616 (2001) 
///
/// In this version we will try using as few temporary matrices as possible, the rest will be stored on the disk via I/O

  int BM = 1; // bemchmark flag

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
  double dE = 0.0;

  vector<Timer> bench_t(10); // timers for different type of operations
  vector<Timer> bench_t2(4);


  if(BM){ bench_t[5].start(); }

//  exit(0);

  // Only 3 auxiliary matrices
  MATRIX* aux1;          aux1 = new MATRIX(Norb,Norb);
  MATRIX* aux2;          aux2 = new MATRIX(Norb,Norb);
  MATRIX* aux3;          aux3 = new MATRIX(Norb,Norb);

/*
  MATRIX* dP;           dP          = new MATRIX(Norb,Norb);  // dP = P_k+1 - P_k
  MATRIX* temp;         temp        = new MATRIX(Norb,Norb);
  MATRIX* P_old_alp;    P_old_alp   = new MATRIX(Norb,Norb);
  MATRIX* P_old_bet;    P_old_bet   = new MATRIX(Norb,Norb);
  MATRIX* P_old;        P_old       = new MATRIX(Norb,Norb);

  MATRIX* P_alp;        P_alp       = new MATRIX(Norb,Norb);
  MATRIX* P_bet;        P_bet       = new MATRIX(Norb,Norb);
  MATRIX* P;            P           = new MATRIX(Norb,Norb);

  MATRIX* P_til_alp;    P_til_alp   = new MATRIX(Norb,Norb);
  MATRIX* P_til_bet;    P_til_bet   = new MATRIX(Norb,Norb);
  MATRIX* P_til;        P_til       = new MATRIX(Norb,Norb);

  MATRIX* Fao_alp;      Fao_alp     = new MATRIX(Norb,Norb);
  MATRIX* Fao_bet;      Fao_bet     = new MATRIX(Norb,Norb);

  MATRIX* Fao_til_alp;  Fao_til_alp = new MATRIX(Norb,Norb);
  MATRIX* Fao_til_bet;  Fao_til_bet = new MATRIX(Norb,Norb);
*/

  Electronic* el_tmp;   el_tmp = new Electronic(el);


  if(BM){ bench_t[5].stop(); }



  // Interface
  el->P->bin_dump("job__P.bin");             
  //*P = *el->P;

  el->P_alp->bin_dump("job__P_alp.bin");     
  //*P_alp = *el->P_alp;

  el->P_bet->bin_dump("job__P_bet.bin");     
  //*P_bet = *el->P_bet;

  el->Fao_alp->bin_dump("job__Fao_alp.bin"); 
  //*Fao_alp = *el->Fao_alp;

  el->Fao_bet->bin_dump("job__Fao_bet.bin"); 
  //*Fao_bet = *el->Fao_bet;
  
  // Old
  el->P_alp->bin_dump("job__P_old_alp.bin"); 
  //*P_old_alp = *P_alp;

  el->P_bet->bin_dump("job__P_old_bet.bin"); 
  //*P_old_bet = *P_bet;

  el->P->bin_dump("job__P_old.bin");         
  //*P_old = *P;

  // Tilda
  // D~_0 = D_0
  el->P_alp->bin_dump("job__P_til_alp.bin"); 
  //*P_til_alp = *P_alp;

  el->P_bet->bin_dump("job__P_til_bet.bin"); 
  //*P_til_bet = *P_bet;

  el->P->bin_dump("job__P_til.bin");         
  //*P_til = *P;
  

  // Initialization:
  // F_0 = F(D_0), F~_0 = F(D~_0) = F_0
  if(BM){ bench_t[1].start(); }
  Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
  if(BM){ bench_t[1].stop(); }
  el_tmp->Fao_alp->bin_dump("job__Fao_til_alp.bin"); 
  //*Fao_til_alp = *el_tmp->Fao_alp;

  el_tmp->Fao_bet->bin_dump("job__Fao_til_bet.bin"); 
  //*Fao_til_bet = *el_tmp->Fao_bet;

  if(BM){ bench_t[0].start(); }

// Old!!!
//  Eelec_prev = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao,el_tmp->Fao_alp);
//  Eelec_prev+= ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao,el_tmp->Fao_bet);

// New!!!
  *aux1 = (*el_tmp->Fao_alp + *el_tmp->P_alp * *el_tmp->dFao_alp_dP_alp + *el_tmp->P_bet * *el_tmp->dFao_alp_dP_bet);
  Eelec_prev = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao, aux1);
  *aux1 = (*el_tmp->Fao_bet + *el_tmp->P_alp * *el_tmp->dFao_bet_dP_alp + *el_tmp->P_bet * *el_tmp->dFao_bet_dP_bet);
  Eelec_prev+= ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao, aux1);


//
  if(BM){ bench_t[0].stop(); }


// Debug:
//  cout<<"In SCF: initital density matrix: P = \n"<<*P<<endl; 

  //=========================== Now enter main SCF cycle ===========================================
  ofstream f1("energy.txt",ios::out);

  cout<<"----------------------- Entering main SCF cycle for RHF calculations --------------------\n"; 

  do{

    cout<<"===============Iteration# "<<iter<<" =====================\n";   

    //---------- Obtain a new density for this iteration -------------------------
    // ODA Step 1: Diagonalize F~_k, assemble D_{k+1} via aufbau (so forcibly set prms.pop_opt = 0) 
    if(BM){ bench_t[2].start(); }
    if(prms.use_damping==0){
/*
      Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, 0, Fao_til_alp, el_tmp->Sao, el_tmp->C_alp, el_tmp->E_alp, el_tmp->bands_alp, el_tmp->occ_alp, P_alp, bench_t2);
      Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, 0, Fao_til_bet, el_tmp->Sao, el_tmp->C_bet, el_tmp->E_bet, el_tmp->bands_bet, el_tmp->occ_bet, P_bet, bench_t2);
      *P = *P_alp + *P_bet;
*/
      exit(0);
    }
    else if(prms.use_damping==1){


      aux1->bin_load("job__Fao_til_alp.bin");
//      aux2->bin_load("job__P_alp.bin");
//      *Fao_til_alp = *aux1;
//      *aux1 = *Fao_til_alp;
      Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, prms.pop_opt, aux1, el_tmp->Sao, el_tmp->C_alp, el_tmp->E_alp, el_tmp->bands_alp, el_tmp->occ_alp, aux2, bench_t2);
      aux2->bin_dump("job__P_alp.bin");
//      Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, prms.pop_opt, Fao_til_alp, el_tmp->Sao, el_tmp->C_alp, el_tmp->E_alp, el_tmp->bands_alp, el_tmp->occ_alp, P_alp, bench_t2);


      aux1->bin_load("job__Fao_til_bet.bin");
//      aux3->bin_load("job__P_bet.bin");
      Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, prms.pop_opt, aux1, el_tmp->Sao, el_tmp->C_bet, el_tmp->E_bet, el_tmp->bands_bet, el_tmp->occ_bet, aux3, bench_t2);
      aux3->bin_dump("job__P_bet.bin");
//      Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, prms.pop_opt, Fao_til_bet, el_tmp->Sao, el_tmp->C_bet, el_tmp->E_bet, el_tmp->bands_bet, el_tmp->occ_bet, P_bet, bench_t2);


      *aux2 += *aux3;
      aux2->bin_dump("job__P.bin");
//      *P = *P_alp + *P_bet;



    }
    if(BM){ bench_t[2].stop(); }


//  exit(0);

// Debug:
//    cout<<"D_{k+1} = D(F~_k) = \n"<<*P<<endl;
//    cout<<"D~_{k} = \n"<<*P_til<<endl;
    if(BM){ bench_t[3].start(); }


    aux1->bin_load("job__P.bin");    cout<<"Pmax = "<<aux1->max_elt()<<endl;
    aux2->bin_load("job__P_old.bin");cout<<"Pold_max = "<<aux2->max_elt()<<endl;
    *aux1 -= *aux2;
    den_err = fabs(aux1->max_elt());


//cout<<"Pmax = "<<P->max_elt()<<endl;    
//cout<<"Pold_max = "<<P_old->max_elt()<<endl;

//    den_err = fabs((*P - *P_old).max_elt());      
    if(BM){ bench_t[3].stop(); }
    cout<<"den_err = "<<den_err<<endl;

//  exit(0);


    if(BM){ bench_t[3].start(); }
    aux1->bin_load("job__P_alp.bin"); aux1->bin_dump("job__P_old_alp.bin");  
    //*P_old_alp = *P_alp;

    aux1->bin_load("job__P_bet.bin"); aux1->bin_dump("job__P_old_bet.bin");  
    //*P_old_bet = *P_bet;

    aux1->bin_load("job__P.bin");     aux1->bin_dump("job__P_old.bin");      
    //*P_old = *P;
    if(BM){ bench_t[3].stop(); }


    // ODA Step 2: Either terminate or continue with the search
    // D_{k+1} - D~_k
    if(den_err<den_tol && fabs(dE)<ene_tol){  ;;  }  
    else{

      if(BM){ bench_t[3].start(); }

      aux1->bin_load("job__P.bin");
      aux2->bin_load("job__P_til.bin");
      *aux1 -= *aux2;
      aux1->bin_dump("job__dP.bin");
//      *dP = *P - *P_til;

      // ODA Step 3: Assemble F_{k+1} = F(D_{k+1})
      el_tmp->P_alp->bin_load("job__P_alp.bin");  
      //*el_tmp->P_alp = *P_alp;

      el_tmp->P_bet->bin_load("job__P_bet.bin");  
      //*el_tmp->P_bet = *P_bet;

      el_tmp->P->bin_load("job__P.bin");          
      //*el_tmp->P = *P;
      if(BM){ bench_t[3].stop(); }


      if(BM){ bench_t[1].start(); }
      Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
      if(BM){ bench_t[1].stop(); }


      if(BM){ bench_t[3].start(); }
      el_tmp->Fao_alp->bin_dump("job__Fao_alp.bin"); 
      //*Fao_alp = *el_tmp->Fao_alp;   

      el_tmp->Fao_bet->bin_dump("job__Fao_bet.bin"); 
      //*Fao_bet = *el_tmp->Fao_bet;
      if(BM){ bench_t[3].stop(); }

  
      // ODA Step 4: Solve the line search problem (via interpolation) or use fixed step
      lamb_min = 0.0;
      if(prms.use_damping){

        if(iter<=prms.damping_start){  lamb_min = 1.0; }
        else{
          lamb_min = prms.damping_const;
        }
        cout<<"Using constant lamb_min = "<<lamb_min<<endl;

      }else{

        cout<<"In oda_scf:  SCF with ODA and storage of matrices to disk is not implemented yet\nExiting...\n";
        exit(0);
/*

      !!!!!!!!!!!!!  Temporary disable !!!!!!!!!!!!!!!!

        cout<<"Line search:\n";   

        double lamb = 1.0;        
        if(BM){ bench_t[3].start(); }
        *el_tmp->P     = *P_til     + lamb * *dP;
        *el_tmp->P_alp = *P_til_alp + lamb * (0.5* *dP);
        *el_tmp->P_bet = *P_til_bet + lamb * (0.5* *dP);
        if(BM){ bench_t[3].stop(); }


        if(BM){ bench_t[1].start(); }
        Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
        if(BM){ bench_t[1].stop(); }

        bench_t[0].start();
        double en1  = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao,el_tmp->Fao_alp);
               en1 += ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao,el_tmp->Fao_bet);
        bench_t[0].stop();



        lamb = 0.5;        
        if(BM){ bench_t[3].start(); }
        *el_tmp->P     = *P_til     + lamb * *dP;
        *el_tmp->P_alp = *P_til_alp + lamb * (0.5* *dP);
        *el_tmp->P_bet = *P_til_bet + lamb * (0.5* *dP);
        if(BM){ bench_t[3].stop(); }


        if(BM){ bench_t[1].start(); }
        Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
        if(BM){ bench_t[1].stop(); }

        bench_t[0].start();
        double en2  = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao,el_tmp->Fao_alp);
               en2 += ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao,el_tmp->Fao_bet);
        bench_t[0].stop();



        if(BM){ bench_t[3].start(); }
        lamb = 0.0;        
        *el_tmp->P     = *P_til     + lamb * *dP;
        *el_tmp->P_alp = *P_til_alp + lamb * (0.5* *dP);
        *el_tmp->P_bet = *P_til_bet + lamb * (0.5* *dP);
        if(BM){ bench_t[3].stop(); }


        if(BM){ bench_t[1].start(); }
        Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
        if(BM){ bench_t[1].stop(); }

        bench_t[0].start();
        double en0  = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao,el_tmp->Fao_alp);
               en0 += ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao,el_tmp->Fao_bet);
        bench_t[0].stop();


     
        // After this operation all Fock matrices are as if no intermediate calculations were performed
        // for Fao_* is actually F(P_*_tmp) = F~_k
        // At this point Fao_ contains F~_k
        cout<<"E(0)= "<<en0<<endl;
        cout<<"E(1/2)= "<<en2<<endl;
        cout<<"E(1)= "<<en1<<endl;

        double _c = en0;
        double _b = 4.0*en2 - en1 - 3.0*en0;
        double _a = en1 - en0 - _b;

        cout<<"Interpolation polynomial = "<<_a<<" * lamb^2 + "<<_b<<" * lamb + "<<_c<<endl;        
        lamb_min = 0.0; 
        if(fabs(_a)>1e-10){ lamb_min = -_b/(2.0*_a); }

        cout<<"lamb_min = "<<lamb_min<<endl;

        if(0<lamb_min && lamb_min<1){
          Eelec = _a*lamb_min*lamb_min + _b*lamb_min + _c;
          cout<<"Functional minimum = "<<Eelec<<endl;
        }
        else{
          lamb_min = (en0<en1)?0.0:1.0;
          Eelec = ((en0<en1)?en0:en1);
          cout<<"infinum at = "<<((en0<en1)?"lamb_min = 0.0":"lamb_min = 1.0")<<" value = "<<Eelec<<endl;
        }

*/
      }// do not use damping

      // ODA Step 5:
      // D~{k+1} = D~{k} + lamb_min * d_k = (1 - lamb_min)*D~_k + lamb_min * D_{k+1}
      if(BM){ bench_t[3].start(); }

      aux1->bin_load("job__P_til_alp.bin"); *aux1 *= (1.0 - lamb_min);
      aux2->bin_load("job__P_alp.bin"); *aux2 *= lamb_min;
      *aux1 += *aux2;
      aux1->bin_dump("job__P_til_alp.bin");

//      *P_til_alp   = (1.0 - lamb_min) * (*P_til_alp) + lamb_min * (*P_alp);


      aux1->bin_load("job__P_til_bet.bin"); *aux1 *= (1.0 - lamb_min);
      aux2->bin_load("job__P_bet.bin"); *aux2 *= lamb_min;
      *aux1 += *aux2;
      aux1->bin_dump("job__P_til_bet.bin");

//      *P_til_bet   = (1.0 - lamb_min) * (*P_til_bet) + lamb_min * (*P_bet);    


      aux1->bin_load("job__P_til_alp.bin"); 
      aux2->bin_load("job__P_til_bet.bin");
      *aux1 += *aux2;
      aux1->bin_dump("job__P_til.bin");

//      *P_til       = *P_til_alp + *P_til_bet;

      // F~{k+1} = (1 - lamb_min)*F~_k + lamb_min * F_{k+1}   - this is original approach, but less general
//      *Fao_til_alp = (1.0 - lamb_min) * (*Fao_til_alp) + lamb_min * (*Fao_alp);
//      *Fao_til_bet = (1.0 - lamb_min) * (*Fao_til_bet) + lamb_min * (*Fao_bet);


      // in fact, F~{k+1} = F(D~_{k+1})  - this is more general approach
      el_tmp->P->bin_load("job__P_til.bin");
      //*el_tmp->P     = *P_til;    

      el_tmp->P_alp->bin_load("job__P_til_alp.bin");
      //*el_tmp->P_alp = *P_til_alp;

      el_tmp->P_bet->bin_load("job__P_til_bet.bin");
      //*el_tmp->P_bet = *P_til_bet;
      if(BM){ bench_t[3].stop(); }


      if(BM){ bench_t[1].start(); }
      Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
      if(BM){ bench_t[1].stop(); }


      el_tmp->Fao_alp->bin_dump("job__Fao_til_alp.bin"); 
//      *Fao_til_alp = *el_tmp->Fao_alp;

      el_tmp->Fao_bet->bin_dump("job__Fao_til_bet.bin"); 
//      *Fao_til_bet = *el_tmp->Fao_bet;


      
    }// else: den_err>=den_tol

   
    // Recompute current energy using extrapolated density matrix
    if(BM){ bench_t[3].start(); }
    el_tmp->P->bin_load("job__P_til.bin");  
    //*el_tmp->P     = *P_til;//     + lamb_min * *dP;

    el_tmp->P_alp->bin_load("job__P_til_alp.bin");
    //*el_tmp->P_alp = *P_til_alp;// + lamb_min * (0.5* *dP);

    el_tmp->P_bet->bin_load("job__P_til_bet.bin");
    //*el_tmp->P_bet = *P_til_bet;// + lamb_min * (0.5* *dP);
    if(BM){ bench_t[3].stop(); }


    if(BM){ bench_t[1].start(); }
    Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix
    if(BM){ bench_t[1].stop(); }

    bench_t[0].start();
// Old !!!
//    Eelec = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao,el_tmp->Fao_alp);
//    Eelec+= ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao,el_tmp->Fao_bet);

// New!!!
    *aux1 = (*el_tmp->Fao_alp + *el_tmp->P_alp * *el_tmp->dFao_alp_dP_alp + *el_tmp->P_bet * *el_tmp->dFao_alp_dP_bet);
    Eelec = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao, aux1);

    *aux1 = (*el_tmp->Fao_bet + *el_tmp->P_alp * *el_tmp->dFao_bet_dP_alp + *el_tmp->P_bet * *el_tmp->dFao_bet_dP_bet);
    Eelec+= ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao, aux1);


    bench_t[0].stop();



    dE = Eelec - Eelec_prev;
    Eelec_prev = Eelec;

     
    if(BM){ bench_t[4].start(); }    

    if(0){  // Suppress excessive output - it takes time

      show_bands(el_tmp->Norb, el_tmp->Nocc_alp, el_tmp->bands_alp, el_tmp->occ_alp);
      show_bands(el_tmp->Norb, el_tmp->Nocc_bet, el_tmp->bands_bet, el_tmp->occ_bet);
        
      cout<<"Total electronic energy = "<<Eelec<<endl;
      cout<<"Mulliken scf orbital populations:\n";
      for(i=0;i<el_tmp->Norb;i++){   cout<<"i= "<<i<<"  pop(gross)= "<<el_tmp->Mull_orb_pop_gross[i]<<"  pop(net)= "<<el_tmp->Mull_orb_pop_net[i]<<endl;  }
      cout<<"Mulliken scf charges:\n";
      for(i=0;i<mol.Nnucl;i++){   cout<<"i= "<<i<<"  q(gross)= "<<mol.Mull_charges_gross[i]<<"  q(net)= "<<mol.Mull_charges_net[i]<<endl;  }
    }

    f1 << iter<<" Eelec= "<<Eelec<<" dE= "<<dE<<" den_err = "<<den_err<<endl;
    //<<" Tr(S*D)= "<<(*el_tmp->Sao * *P_til).tr()<<endl;
    if(BM){ bench_t[4].stop(); }    

    iter++;    

 
  }while(iter<Niter && (den_err>den_tol || fabs(dE)>ene_tol) );

  f1.close();


//  exit(0);



  if(prms.do_annihilate==1){ 
    aux1->bin_load("job__P_til_alp.bin");
    aux2->bin_load("job__P_til_bet.bin");
    annihilate(Nocc_alp,Nocc_bet,aux1,aux2);
    //annihilate(Nocc_alp,Nocc_bet,P_til_alp,P_til_bet);
  }

  if(BM){ bench_t[3].start(); }

  el->P_alp->bin_load("job__P_til_alp.bin");  //*el->P_alp = *P_til_alp;
  el->P_bet->bin_load("job__P_til_bet.bin");  //*el->P_bet = *P_til_bet;
  *el->P     = *el->P_alp + *el->P_bet;

  el->Fao_alp->bin_load("job__Fao_til_alp.bin"); //*el->Fao_alp = *Fao_til_alp;
  el->Fao_bet->bin_load("job__Fao_til_bet.bin"); //*el->Fao_bet = *Fao_til_bet;
  if(BM){ bench_t[3].stop(); }


  if(BM){ bench_t[1].start(); }
  Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el, el0,mem); // el - now contains an updated Fock matrix
  if(BM){ bench_t[1].stop(); }

  // Test: At this point F(P~) = F~, so it is valid to use P~ to construct F~ via normal rules
  //cout<<"Difference of Fock matrices: \n"<<*el->Fao_alp - *Fao_til_alp<<endl;

  // Update eigenvalues and eigenvectors of final Fock matrix, but do not modify the density matrix:
  if(BM){ bench_t[2].start(); }
  Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, 0, el->Fao_alp, el->Sao, el->C_alp, el->E_alp, el->bands_alp, el->occ_alp, aux1, bench_t2);  
  Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, 0, el->Fao_bet, el->Sao, el->C_bet, el->E_bet, el->bands_bet, el->occ_bet, aux1, bench_t2);  
  if(BM){ bench_t[2].stop(); }

  bench_t[0].start();

// Old!!!
//  Eelec = ::energy_elec(Norb,el->P_alp,el->Hao,el->Fao_alp) + ::energy_elec(Norb,el->P_bet,el->Hao,el->Fao_bet);

// New!!!
    *aux1 = (*el->Fao_alp + *el->P_alp * *el->dFao_alp_dP_alp + *el->P_bet * *el->dFao_alp_dP_bet);
    Eelec = ::energy_elec(Norb,el->P_alp,el->Hao, aux1);

    *aux1 = (*el->Fao_bet + *el->P_alp * *el->dFao_bet_dP_alp + *el->P_bet * *el->dFao_bet_dP_bet);
    Eelec+= ::energy_elec(Norb,el->P_bet,el->Hao, aux1);



  bench_t[0].stop();

/*  DEBUG
  cout<<"In SCF:\n";
  cout<<"F*C = \n"<<*el->Fao_alp * *el->C_alp<<endl;
  cout<<"S*C*E = \n"<<*el->Sao * *el->C_alp * *el->E_alp<<endl;
  cout<<"MO orthogonality:\n";

  // MO
  for(i=0;i<el->Norb;i++){
    for(int j=0;j<el->Norb;j++){

      double res = 0.0;

      int ii = el->bands_alp[i].first;
      int jj = el->bands_alp[j].first;
       

      for(int a=0;a<el->Norb;a++){    // over basis of fragment A at time t
        for(int b=0;b<el->Norb;b++){  // over basis of fragment B at time t

          res += el->C_alp->M[a*el->Norb+ii] * el->C_alp->M[b*el->Norb+jj] * el->Sao->M[a*el->Norb+b];
                
        }// for b
      }// for a
      cout<<res<<"  ";
    }
    cout<<endl;
  }
*/


  // Clean up the memory
  if(BM){ bench_t[5].start(); }
/*
  delete dP;
  delete temp;
  delete P_old_alp;  
  delete P_old_bet;  
  delete P_old;      

  delete P_alp;     
  delete P_bet;     
  delete P;         

  delete P_til_alp;   
  delete P_til_bet;   
  delete P_til;       

  delete Fao_alp;     
  delete Fao_bet;     

  delete Fao_til_alp; 
  delete Fao_til_bet; 
*/

  delete aux1;
  delete aux2;
  delete aux3;
 

  el_tmp->~Electronic();
  if(BM){ bench_t[5].stop(); }


  if(BM){
    cout<<"Time for energy calculation = "<<bench_t[0].show()<<endl;
    cout<<"Time for Fock matrix formaion = "<<bench_t[1].show()<<endl;
    cout<<"Time for Fock diagonalization and density matrix formation = "<<bench_t[2].show()<<endl;
    cout<<"   - eigensolver    = "<<bench_t2[0].show()<<endl;
    cout<<"   - sorting        = "<<bench_t2[1].show()<<endl;
    cout<<"   - populate       = "<<bench_t2[2].show()<<endl;
    cout<<"   - density matrix = "<<bench_t2[3].show()<<endl;
    cout<<"Time for matrix operations = "<<bench_t[3].show()<<endl;
    cout<<"Time for output = "<<bench_t[4].show()<<endl;
    cout<<"Time for allocation/deallocation = "<<bench_t[5].show()<<endl;
 
   

  }//

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



  return Eelec;

} // oda_disc




double scf_diis_fock(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                     vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                     Electronic* el,Electronic* el0, Memory* mem){

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


// !!!TEMPORARY FORGET ABOUT THIS FUNCTION!!!

double scf_diis_dm(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                   vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                   Electronic* el,Electronic* el0, Memory* mem){

  exit(0);

/// This function implements the SCF based on the DIIS with density matrix mixing(extrapolation)
/// This can eventually be along ODA lines
/// if prms.use_diis = 0 => routine reduces to standard SCF
/// 
/// I tried following the recipe from:
/// [1] Pulay, P. "Convergence Acceleration of Iterative Sequences. The Case of SCF Iteration" Chem. Phys. Lett. 1980, 73, 393
/// [2] Pulay, P. "Improved SCF Convergence Acceleration" J. Comp. Chem. 1982, 3, 556


  int i;
  double lamb_min;
  std::string eigen_method="generalized";
//  if(prms.hamiltonian=="indo"){ eigen_method = "standard"; }


  //----------- Control parameters ---------
  int iter = 0;
  int Niter = prms.Niter;

  int Norb = el->Norb;
  int Nocc_alp = el->Nocc_alp;
  int Nocc_bet = el->Nocc_bet;

  double den_tol = prms.den_tol;  
  double den_err = 2.0*den_tol;
  double Eelec = 0.0;

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
//    *err_alp = ((*el_tmp->C_alp).T()) * *err_alp * (*el_tmp->C_alp);
//    *err_bet = ((*el_tmp->C_bet).T()) * *err_bet * (*el_tmp->C_bet);

    // Scale by a common factor:
//    double alp_err = ((*err_alp).T() * (*err_alp)).tr();
//    double bet_err = ((*err_bet).T() * (*err_bet)).tr();
//    double tot_err = 0.5*(alp_err + bet_err);
//    if(alp_err>0.0) { *err_alp *= (tot_err/alp_err); }
//    if(bet_err>0.0) { *err_bet *= (tot_err/bet_err); }



    // Error in commutation:  max_elt |FPS - SPF|_{alp} + max_elt |FPS - SPF|_{bet}
    den_err = err_alp->max_elt() +  err_bet->max_elt();   
    cout<<"commutation error = "<<den_err<<endl;


    // Step 3: Either terminate or continue with the search
    if(den_err<den_tol){  ;;  }  
    else{

      if(prms.use_diis && start_diis){

        // DIIS Step 1: Store error matrix and current Fock matrices
        diis_alp->add_diis_matrices(el_tmp->P_alp, err_alp);
        diis_bet->add_diis_matrices(el_tmp->P_bet, err_bet);  


        if(diis_delay >= diis_alp->N_diis_max){

          // DIIS Setp 2: Compute extrapolated Fock matrix and replace the current Fock matrix by the extrapolated one
          diis_alp->extrapolate_matrix(el_tmp->P_alp);
          diis_bet->extrapolate_matrix(el_tmp->P_bet);
          *el_tmp->P = *el_tmp->P_alp + *el_tmp->P_bet;

          Hamiltonian_Fock(prms,modprms,mol,fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el_tmp - now contains an updated Fock matrix

        }

        diis_delay++;

      }// use_diis && start_diis

      //---------- Obtain a new density for this iteration -------------------------
      // Step 4: Diagonalize F_k, assemble corresponding new density - D_{k+1} via aufbau or Fermi populations
      Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, prms.pop_opt, el_tmp->Fao_alp, el_tmp->Sao, el_tmp->C_alp, el_tmp->E_alp, el_tmp->bands_alp, el_tmp->occ_alp, el_tmp->P_alp, bench_t2);
      Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, prms.pop_opt, el_tmp->Fao_bet, el_tmp->Sao, el_tmp->C_bet, el_tmp->E_bet, el_tmp->bands_bet, el_tmp->occ_bet, el_tmp->P_bet, bench_t2);
      *el_tmp->P = *el_tmp->P_alp + *el_tmp->P_bet;


      
    }// else: den_err>=den_tol

   
    Eelec = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao,el_tmp->Fao_alp);
    Eelec+= ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao,el_tmp->Fao_bet);

        
    cout<<"Total electronic energy = "<<Eelec<<endl;
    cout<<"Mulliken scf orbital populations:\n";
    for(i=0;i<el_tmp->Norb;i++){   cout<<"i= "<<i<<"  pop(gross)= "<<el_tmp->Mull_orb_pop_gross[i]<<"  pop(net)= "<<el_tmp->Mull_orb_pop_net[i]<<endl;  }
    cout<<"Mulliken scf charges:\n";
    for(i=0;i<mol.Nnucl;i++){   cout<<"i= "<<i<<"  q(gross)= "<<mol.Mull_charges_gross[i]<<"  q(net)= "<<mol.Mull_charges_net[i]<<endl;  }

    f1 << iter<<"  "<<Eelec<<" Tr(S*D)= "<<(*el_tmp->Sao * *el_tmp->P).tr()<<endl;
    iter++;    

 
  }while(iter<Niter && den_err>den_tol);

  f1.close();


  if(fabs(den_err)>den_tol){
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    cout<<"!!!!! Error: Convergence is not achieved after "<<Niter<<" iterations\n den_err = "<<den_err<<" !!!!!\n";
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

