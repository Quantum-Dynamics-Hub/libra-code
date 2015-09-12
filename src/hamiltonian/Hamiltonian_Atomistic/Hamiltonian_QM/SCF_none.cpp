
#include "SCF.h"


namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_qm{



double scf_none(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM){

/// This function implements the SCF based on the optimal damping algorithm (ODA)
/// which uses fractional occupation numbers, leading to robust convergence in difficult cases
/// See more details in: 
/// [1] Kudin K.N.; Scuseria, G.E.; Cances, E. J. Chem. Phys. 116, 8255 (2002)
/// [2] Cances J. Chem. Phys. 114, 10616 (2001) 
/// BM - benchmark flag

  int DF = 0;
  int Norb = el->Norb;
  int Nocc_alp = el->Nocc_alp;
  int Nocc_bet = el->Nocc_bet;
  std::string eigen_method="generalized";

  vector<Timer> bench_t(10); // timers for different type of operations
  vector<Timer> bench_t2(4);


  MATRIX* P_alp_old;        P_alp_old       = new MATRIX(Norb,Norb);
  MATRIX* P_bet_old;        P_bet_old       = new MATRIX(Norb,Norb);


  //===============  Initialization =======================

  Hamiltonian_Fock(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

  if(DF){
    cout<<"Fock matrix at first iteration (alp)\n";
    el->get_Fao_alp().show_matrix();
  }

//  double E = (energy_elec(el->get_P_alp(), el->get_Hao(), el->get_Fao_alp())
//            + energy_elec(el->get_P_bet(), el->get_Hao(), el->get_Fao_bet()));
  double E = (energy_elec(el->P_alp, el->Hao, el->Fao_alp) + energy_elec(el->P_bet, el->Hao, el->Fao_bet)  );

  cout<<"Initial energy = "<< E<<endl;

  double E_old = E;
  double e_err = 2.0*prms.etol;
  double d_err = 2.0*prms.den_tol;

  *P_alp_old = *el->P_alp;
  *P_bet_old = *el->P_bet;


  //===============  Now to SCF iterations =======================
  int i = 0;
  int run = 1;
  while(run){
    

    Fock_to_P(Norb, Nocc_alp, 1, Nocc_alp, eigen_method, prms.pop_opt, el->Fao_alp, el->Sao, el->C_alp, el->E_alp, el->bands_alp, el->occ_alp, el->P_alp, bench_t2);
    Fock_to_P(Norb, Nocc_bet, 1, Nocc_bet, eigen_method, prms.pop_opt, el->Fao_bet, el->Sao, el->C_bet, el->E_bet, el->bands_bet, el->occ_bet, el->P_bet, bench_t2);
    *el->P = *el->P_alp + *el->P_bet;

    Hamiltonian_Fock(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);


    d_err = fabs((*el->P_alp - *P_alp_old).max_elt()) + fabs((*el->P_bet - *P_bet_old).max_elt());

    *P_alp_old = *el->P_alp;
    *P_bet_old = *el->P_bet;

    E_old = E;
    E = (energy_elec(el->P_alp, el->Hao, el->Fao_alp) + energy_elec(el->P_bet, el->Hao, el->Fao_bet)  );

    e_err = fabs(E_old - E);

    cout<<"Iteration "<<i<<" e_err = "<<e_err<<" d_err = "<<d_err<<" E_el = "<<E<<endl;

    if(i>prms.Niter){
        run = 0;
        cout<<"Convergence is not achieved in "<<prms.Niter<<" iterations\n";
    }
    if(e_err<prms.etol && d_err<prms.den_tol){
        run = 0;
        cout<<"Success: Convergence is achieved\n";
        cout<<"Electronic energy = "<<E<<endl;
    }

    i = i + 1;
  }// while

}


double scf_none(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM
){

  scf_none(&el,syst,basis_ao,  prms,modprms,  atom_to_ao_map,ao_to_atom_map, BM);
}




}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


    

