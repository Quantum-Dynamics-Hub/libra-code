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
/**
  \file Excitations.cpp
  \brief The file implements functions for creation of the excitation objects for excited state calculations
    
*/

#include "Excitations.h"
#include "Bands.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libcalculators namespace
namespace libcalculators{

void excite(int I, int J, vector< pair<int,double> >& occ_ini, vector< pair<int,double> >& occ_fin){
/**
  \brief Create a list representing a I --> J ecxitation

  basically, vector< pair<int,double> > datastructure represents Slater determinant (by selecting ordering of orbitals in
  the active space)

  \param[in]  occ_ini  initial occupation list 
  \param[out] occ_fin  final occupation list
*/


  int Norb = occ_ini.size();

  if(I<0){ std::cout<<"Error: source orbital index("<<I<<") can not be negative\n";  exit(0); }
  if(J<0){ std::cout<<"Error: target orbital index("<<J<<") can not be negative\n";  exit(0); }
  if(I>=Norb){ std::cout<<"Error: source orbital index("<<I<<") exceeds the number of orbitals in the active space("<<Norb<<")\n";  exit(0); }
  if(J>=Norb){ std::cout<<"Error: target orbital index("<<J<<") exceeds the number of orbitals in the active space("<<Norb<<")\n";  exit(0); }


  if(occ_fin.size()>0){  occ_fin.clear(); }                    // clean result

  for(int i=0;i<Norb;i++){  occ_fin.push_back(occ_ini[i]);  }  // copy initial to final

  // Swap populations I-th and J-th orbitals
  double pop = occ_fin[I].second; 
  occ_fin[I].second = occ_fin[J].second;
  occ_fin[J].second = pop;
   
}

boost::python::list excite(int I, int J, boost::python::list occ_ini){
/**
  \brief Create a list representing a I --> J ecxitation (Python-friendly)

  basically, vector< pair<int,double> > datastructure represents Slater determinant (by selecting ordering of orbitals in
  the active space)
  Returns the final occupation list

  \param[in]  occ_ini  initial occupation list 
*/

  vector< pair<int,double> > occ_i;
  vector< pair<int,double> > occ_f;

  convert_1(occ_ini, occ_i);

  excite(I,J,occ_i, occ_f);

  return convert_2(occ_f);

}



void excite(int Norb, excitation& ex, 
            int Nocc_alp, vector< pair<int,double> >& occ_alp,
            int Nocc_bet, vector< pair<int,double> >& occ_bet){
// N-order excitation
// 
// from[0] -> to[0]
// from[1] -> to[0] 
// ...
// from[N] -> to[N]
// Example: from = [ 0A, 0B] to = [1A, 1B] - double excitation of 0A -> 1A & 0B -> 1B
// Here 0 = HOMO, -1 = HOMO-1, 1 = LUMO, 2 = LUMO+1, etc.

// This function annihilates on <from> whatever is there and adds it to <to> whatever exists there already
// This function changes occ_alp and occ_bet - these are both input and output variables

  cout<<"In excite(...)\n";
  cout<<"Excitation size = "<<ex.size<<endl;
  for(int i=0;i<ex.size;i++){

    int source_indx = 0;
    int target_indx = 0;
    double source_pop = 0.0;
    double target_pop = 0.0;

    if(ex.from_spin[i]==1){
      source_pop = occ_alp[Nocc_alp + ex.from_orbit[i] - 1].second;

      cout<<"Source is alpha\n";
      cout<<"Source indx = "<<Nocc_alp + ex.from_orbit[i] - 1<<endl;
      cout<<"Source population = "<<source_pop<<endl;
    }
    else if(ex.from_spin[i]==-1){
      source_pop = occ_bet[Nocc_bet + ex.from_orbit[i] - 1].second;

      cout<<"Source is beta\n";
      cout<<"Source indx = "<<Nocc_bet + ex.from_orbit[i] - 1<<endl;
      cout<<"Source population = "<<source_pop<<endl;

    }


    if(ex.to_spin[i]==1){
      target_pop = occ_alp[Nocc_alp + ex.to_orbit[i] - 1].second;

      cout<<"Target is alpha\n";
      cout<<"Target indx = "<<Nocc_alp + ex.to_orbit[i] - 1<<endl;
      cout<<"Target population = "<<target_pop<<endl;

    }
    else if(ex.to_spin[i]==-1){
      target_pop = occ_bet[Nocc_bet + ex.to_orbit[i] - 1].second;

      cout<<"Target is beta\n";
      cout<<"Target indx = "<<Nocc_bet + ex.to_orbit[i] - 1<<endl;
      cout<<"Target population = "<<target_pop<<endl;

    }


    if(source_pop<=(1.0-target_pop)){    // There is plenty of room in target
      cout<<"Population to transfer is: "<<source_pop<<endl;
      // Annihilate population on source orbitals
      if(ex.from_spin[i]==1){
        occ_alp[Nocc_alp + ex.from_orbit[i] - 1].second -= source_pop; // take all population off     
      }
      else if(ex.from_spin[i]==-1){
        occ_bet[Nocc_bet + ex.from_orbit[i] - 1].second -= source_pop; // take all population off     
      }

      // Create population on target orbitals
      if(ex.to_spin[i]==1){
        occ_alp[Nocc_alp + ex.to_orbit[i] - 1].second += source_pop;
      }
      else if(ex.to_spin[i]==-1){
        occ_bet[Nocc_bet + ex.to_orbit[i] - 1].second += source_pop;
      }

    }
    else{  // In this case there is more that can be take from the source than can be put to the target
      cout<<"Population to transfer is: "<<(1.0-target_pop)<<endl;

      // Annihilate population on source orbitals
      if(ex.from_spin[i]==1){
        occ_alp[Nocc_alp + ex.from_orbit[i] - 1].second -= (1.0-target_pop); // take all population off     
      }
      else if(ex.from_spin[i]==-1){
        occ_bet[Nocc_bet + ex.from_orbit[i] - 1].second -= (1.0-target_pop); // take all population off     
      }

      // Create population on target orbitals
      if(ex.to_spin[i]==1){
        occ_alp[Nocc_alp + ex.to_orbit[i] - 1].second += (1.0-target_pop);
      }
      else if(ex.to_spin[i]==-1){
        occ_bet[Nocc_bet + ex.to_orbit[i] - 1].second += (1.0-target_pop);
      }

    }
    
    
  }// for i

} // excite(...)



/*

void compute_excitations(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                         vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                         Electronic* el,Electronic* el0, Memory* mem){

//
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
*/



}// namespace libcalculators

}// liblibra


