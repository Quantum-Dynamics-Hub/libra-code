
void compute_vertical_ip(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                         vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                         Electronic* el,Electronic* el0, Memory* mem){

//
// This function relies on the variables set up in main scf cycle

  int e,a,b,i,ii1,ii2;
  VECTOR dIdA, dIdB, TV, r_ab;
  VECTOR mu; mu = 0.0;



  Electronic* el_tmp;   el_tmp = new Electronic(el);
  Nuclear mol_tmp;      mol_tmp = mol;

  int Norb = el_tmp->Norb;
  double E_ground = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao,el_tmp->Fao_alp)
                  + ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao,el_tmp->Fao_bet);
  cout<<"E_ground = "<<E_ground<<endl;

  // Now compute vertical excitation energies
  vector< pair<int,double> > occ_alp_grnd(Norb,pair<int,double>(0,0.0)); occ_alp_grnd = el_tmp->occ_alp;
  vector< pair<int,double> > occ_bet_grnd(Norb,pair<int,double>(0,0.0)); occ_bet_grnd = el_tmp->occ_bet;
  vector< pair<int,double> > docc(Norb,pair<int,double>(0,0.0)); 


  // Important!!! Do not recompute density from MO populations - because the final density may in fact be
  // an ensemble density (of last two few iterates - e.g. if DIIS or damping is used)
  // So we only compute changes of density matrix that correspond to change of MO populations within reasonable limits
  MATRIX* dP_;    dP_ = new MATRIX(Norb,Norb); *dP_ = 0.0;


  // For oscillator strength calculations
  MATRIX* mux;  mux = new MATRIX(Norb,Norb); *mux = 0.0;
  MATRIX* muy;  muy = new MATRIX(Norb,Norb); *muy = 0.0;
  MATRIX* muz;  muz = new MATRIX(Norb,Norb); *muz = 0.0;
  VECTOR mu_nucl;

  compute_dipole_moments(prms,modprms,mol_tmp, fragment, basis_fo, basis_ao, at_orbitals, el_tmp, mem, mu_nucl, mux, muy, muz);




  // !!! === Compute vertical excitations === !!!
  for(e=0;e<prms.num_excitations;e++){ 

    cout<<"excitation "<<e<<"  \n";
    el_tmp->occ_alp = occ_alp_grnd;
    el_tmp->occ_bet = occ_bet_grnd;

    excite(Norb, prms.excitations[e], el_tmp->Nocc_alp, el_tmp->occ_alp, el_tmp->Nocc_bet, el_tmp->occ_bet); // ground state excitation


    cout<<"Excitation: "<<prms.excitations[e].from_orbit[0]<<"  "<<prms.excitations[e].from_spin[0]<<"  -->  "
                        <<prms.excitations[e].to_orbit[0]<<"  "<<prms.excitations[e].to_spin[0]<<endl;


    // Update delta of density matrix
    //!!! But keep in mind that WFC is not always equivalent to density matrix
    for(i=0;i<Norb;i++){ docc[i] = pair<int,double>(i, el_tmp->occ_alp[i].second - occ_alp_grnd[i].second); }
    compute_density_matrix(docc, el_tmp->C_alp, dP_);    
    *el_tmp->P_alp = *el->P_alp  +  *dP_;
        
    for(i=0;i<Norb;i++){ docc[i] = pair<int,double>(i, el_tmp->occ_bet[i].second - occ_bet_grnd[i].second); }
    compute_density_matrix(docc, el_tmp->C_bet, dP_);
    *el_tmp->P_bet = *el->P_bet + *dP_;

    *el_tmp->P = *el_tmp->P_alp + *el_tmp->P_bet;  // this is an excited state density matrix (up to our approximation)


    // Update Fock matrix
    Hamiltonian_Fock(prms,modprms,mol_tmp, fragment, basis_fo,basis_ao,at_orbitals,el_tmp,el0,mem); // el - now contains an updated Fock matrix



    // Compute energy
    double E_ext = ::energy_elec(Norb,el_tmp->P_alp,el_tmp->Hao,el_tmp->Fao_alp)
                 + ::energy_elec(Norb,el_tmp->P_bet,el_tmp->Hao,el_tmp->Fao_bet);
    prms.excitations[e].Energy = E_ext - E_ground;


    // Compute oscillator strength
    // Compute transition dipole moment for given excitation:
    //double x; 
    //!!!!!!!!!! ASSUME: Excitaions are only on alpha-channel !!!!!!!!!!!!
    // Alp component
    ii1 = el_tmp->Nocc_alp + prms.excitations[e].from_orbit[0] - 1; // index of the source MO
    ii2 = el_tmp->Nocc_alp + prms.excitations[e].to_orbit[0] - 1;   // index of the target MO

    for(a=0;a<Norb;a++){
      for(b=0;b<Norb;b++){

        double p = el_tmp->C_alp->M[a*Norb+ii1] *  el_tmp->C_alp->M[b*Norb+ii2];

        mu.x += mux->M[a*Norb+b]*p;
        mu.y += muy->M[a*Norb+b]*p;
        mu.z += muz->M[a*Norb+b]*p;

      }// for b
    }// a

    double f =  (2.0/3.0)*(prms.excitations[e].Energy)*mu.length2(); // a.u.
    prms.excitations[e].f = f;

    cout<<"Excitation a= "<<e<<" E-E0,eV= "<<prms.excitations[e].Energy/eV<<" oscillator strength= "<<f<<" Tr(S*D)= "<<(*el_tmp->Sao * *el_tmp->P).tr()<<endl;



  }// for e



  // Free memory allocated for temporary matrices
  delete dP_;
  delete mux;
  delete muy;
  delete muz;

  el_tmp->~Electronic();

//  exit(0);

}

