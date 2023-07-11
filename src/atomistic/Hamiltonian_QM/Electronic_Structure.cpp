/*********************************************************************************
* Copyright (C) 2014-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Electronic_Structure.cpp
  \brief This file implements the Electronic_Structure class - a container for electronic variables.
  Some methods for handling the Electronic_Structure objects are also described here

*/

#include "Electronic_Structure.h"
#include "../../calculators/libcalculators.h"

/// liblibra namespace
namespace liblibra{

using namespace libcalculators;

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{



/// Implementation of the copy constructor of the class Electronic
Electronic_Structure::Electronic_Structure(const Electronic_Structure& obj){

    Norb = obj.Norb;
    Nocc_alp = obj.Nocc_alp;
    Nocc_bet = obj.Nocc_bet;
    Nelec = obj.Nelec;

    P_alp = new MATRIX(Norb,Norb);   *P_alp   = *obj.P_alp;
    P_bet = new MATRIX(Norb,Norb);   *P_bet   = *obj.P_bet;
    P     = new MATRIX(Norb,Norb);   *P       = *obj.P;
    C_alp = new MATRIX(Norb,Norb);   *C_alp   = *obj.C_alp;
    C_bet = new MATRIX(Norb,Norb);   *C_bet   = *obj.C_bet;
    Sao   = new MATRIX(Norb,Norb);   *Sao     = *obj.Sao;
    Hao   = new MATRIX(Norb,Norb);   *Hao     = *obj.Hao;
    Fao_alp = new MATRIX(Norb,Norb); *Fao_alp = *obj.Fao_alp;
    Fao_bet = new MATRIX(Norb,Norb); *Fao_bet = *obj.Fao_bet;
    dFao_alp_dP_alp = new MATRIX(Norb,Norb);  *dFao_alp_dP_alp = *obj.dFao_alp_dP_alp;
    dFao_alp_dP_bet = new MATRIX(Norb,Norb);  *dFao_alp_dP_bet = *obj.dFao_alp_dP_bet;
    dFao_bet_dP_alp = new MATRIX(Norb,Norb);  *dFao_bet_dP_alp = *obj.dFao_bet_dP_alp;
    dFao_bet_dP_bet = new MATRIX(Norb,Norb);  *dFao_bet_dP_bet = *obj.dFao_bet_dP_bet;

    E_alp = new MATRIX(Norb,Norb);   *E_alp   = *obj.E_alp;
    E_bet = new MATRIX(Norb,Norb);   *E_bet   = *obj.E_bet;

    Mull_orb_pop_net = obj.Mull_orb_pop_net;
    Mull_orb_pop_gross = obj.Mull_orb_pop_gross;

    bands_alp = obj.bands_alp;
    bands_bet = obj.bands_bet;
    occ_alp = obj.occ_alp;    
    occ_bet = obj.occ_bet;    

}

Electronic_Structure& Electronic_Structure::operator=(const Electronic_Structure& ob){

  if(this == &ob){  return *this;    }
  else{
    *this = Electronic_Structure(ob);
    return *this;
  }

}



/// Implementation of the constructor from the address pointing to an existing object of class Electronic
Electronic_Structure::Electronic_Structure(Electronic_Structure* obj){

    Norb = obj->Norb;
    Nocc_alp = obj->Nocc_alp;
    Nocc_bet = obj->Nocc_bet;
    Nelec = obj->Nelec;

    P_alp = new MATRIX(Norb,Norb);   *P_alp   = *obj->P_alp;
    P_bet = new MATRIX(Norb,Norb);   *P_bet   = *obj->P_bet;
    P     = new MATRIX(Norb,Norb);   *P       = *obj->P;
    C_alp = new MATRIX(Norb,Norb);   *C_alp   = *obj->C_alp;
    C_bet = new MATRIX(Norb,Norb);   *C_bet   = *obj->C_bet;
    Sao   = new MATRIX(Norb,Norb);   *Sao     = *obj->Sao;
    Hao   = new MATRIX(Norb,Norb);   *Hao     = *obj->Hao;
    Fao_alp = new MATRIX(Norb,Norb); *Fao_alp = *obj->Fao_alp;
    Fao_bet = new MATRIX(Norb,Norb); *Fao_bet = *obj->Fao_bet;
    dFao_alp_dP_alp = new MATRIX(Norb,Norb);  *dFao_alp_dP_alp = *obj->dFao_alp_dP_alp;
    dFao_alp_dP_bet = new MATRIX(Norb,Norb);  *dFao_alp_dP_bet = *obj->dFao_alp_dP_bet;
    dFao_bet_dP_alp = new MATRIX(Norb,Norb);  *dFao_bet_dP_alp = *obj->dFao_bet_dP_alp;
    dFao_bet_dP_bet = new MATRIX(Norb,Norb);  *dFao_bet_dP_bet = *obj->dFao_bet_dP_bet;
    E_alp = new MATRIX(Norb,Norb);   *E_alp   = *obj->E_alp;
    E_bet = new MATRIX(Norb,Norb);   *E_bet   = *obj->E_bet;

    Mull_orb_pop_net = obj->Mull_orb_pop_net;
    Mull_orb_pop_gross = obj->Mull_orb_pop_gross;

    bands_alp = obj->bands_alp;
    bands_bet = obj->bands_bet;
    occ_alp = obj->occ_alp;    
    occ_bet = obj->occ_bet;    

}


void Electronic_Structure::check_matrix_dimensionas(MATRIX* A, MATRIX& B, std::string A_name, std::string B_name, std::string func_name){
  if(B.n_cols != A->n_cols){
    cout<<"In "<<func_name<<"\n";
    cout<<"Number of cols of matrix "<<A_name<<" is not equal to the number of cols of matrix "<<B_name<<"\n";
    exit(0);
  }
  if(B.n_rows != A->n_rows){
    cout<<"In "<<func_name<<"\n";
    cout<<"Number of rows of matrix "<<A_name<<" is not equal to the number of rows of matrix "<<B_name<<"\n";
    exit(0);
  }

}

void Electronic_Structure::set_P_alp(MATRIX& x_){
  check_matrix_dimensionas(P_alp, x_, "Electronic::P_alp", "x_", "Electronic::set_P_alp");
  *P_alp = x_;
}
void Electronic_Structure::set_P_bet(MATRIX& x_){
  check_matrix_dimensionas(P_bet, x_, "Electronic::P_bet", "x_", "Electronic::set_P_bet");
  *P_bet = x_;
}
void Electronic_Structure::set_P(MATRIX& x_){
  check_matrix_dimensionas(P, x_, "Electronic::P", "x_", "Electronic::set_P");
  *P = x_;
}
MATRIX Electronic_Structure::get_P_alp(){  return *P_alp; }
MATRIX Electronic_Structure::get_P_bet(){  return *P_bet; }
MATRIX Electronic_Structure::get_P(){  return *P; }


void Electronic_Structure::set_C_alp(MATRIX& x_){
  check_matrix_dimensionas(C_alp, x_, "Electronic::C_alp", "x_", "Electronic::set_C_alp");
  *C_alp = x_;
}
void Electronic_Structure::set_C_bet(MATRIX& x_){
  check_matrix_dimensionas(C_bet, x_, "Electronic::C_bet", "x_", "Electronic::set_C_bet");
  *C_bet = x_;
}
MATRIX Electronic_Structure::get_C_alp(){  return *C_alp; }
MATRIX Electronic_Structure::get_C_bet(){  return *C_bet; }


void Electronic_Structure::set_Sao(MATRIX& x_){
  check_matrix_dimensionas(Sao, x_, "Electronic::Sao", "x_", "Electronic::set_Sao");
  *Sao = x_;
}
void Electronic_Structure::set_Hao(MATRIX& x_){
  check_matrix_dimensionas(Hao, x_, "Electronic::Hao", "x_", "Electronic::set_Hao");
  *Hao = x_;
}
MATRIX Electronic_Structure::get_Sao(){  return *Sao; }
MATRIX Electronic_Structure::get_Hao(){  return *Hao; }


void Electronic_Structure::set_Fao_alp(MATRIX& x_){
  check_matrix_dimensionas(Fao_alp, x_, "Electronic::Fao_alp", "x_", "Electronic::set_Fao_alp");
  *Fao_alp = x_;
}
void Electronic_Structure::set_Fao_bet(MATRIX& x_){
  check_matrix_dimensionas(Fao_bet, x_, "Electronic::Fao_bet", "x_", "Electronic::set_Fao_bet");
  *Fao_bet = x_;
}
MATRIX Electronic_Structure::get_Fao_alp(){  return *Fao_alp; }
MATRIX Electronic_Structure::get_Fao_bet(){  return *Fao_bet; }



void Electronic_Structure::set_dFao_alp_dP_alp(MATRIX& x_){
  check_matrix_dimensionas(dFao_alp_dP_alp, x_, "Electronic::dFao_alp_dP_alp", "x_", "Electronic::set_dFao_alp_dP_alp");
  *dFao_alp_dP_alp = x_;
}
void Electronic_Structure::set_dFao_alp_dP_bet(MATRIX& x_){
  check_matrix_dimensionas(dFao_alp_dP_bet, x_, "Electronic::dFao_alp_dP_bet", "x_", "Electronic::set_dFao_alp_dP_bet");
  *dFao_alp_dP_bet = x_;
}
void Electronic_Structure::set_dFao_bet_dP_alp(MATRIX& x_){
  check_matrix_dimensionas(dFao_bet_dP_alp, x_, "Electronic::dFao_bet_dP_alp", "x_", "Electronic::set_dFao_bet_dP_alp");
  *dFao_bet_dP_alp = x_;
}
void Electronic_Structure::set_dFao_bet_dP_bet(MATRIX& x_){
  check_matrix_dimensionas(dFao_bet_dP_bet, x_, "Electronic::dFao_bet_dP_bet", "x_", "Electronic::set_dFao_bet_dP_bet");
  *dFao_bet_dP_bet = x_;
}
MATRIX Electronic_Structure::get_dFao_alp_dP_alp(){  return *dFao_alp_dP_alp; }
MATRIX Electronic_Structure::get_dFao_alp_dP_bet(){  return *dFao_alp_dP_bet; }
MATRIX Electronic_Structure::get_dFao_bet_dP_alp(){  return *dFao_bet_dP_alp; }
MATRIX Electronic_Structure::get_dFao_bet_dP_bet(){  return *dFao_bet_dP_bet; }



void Electronic_Structure::set_E_alp(MATRIX& x_){
  check_matrix_dimensionas(E_alp, x_, "Electronic::E_alp", "x_", "Electronic::set_E_alp");
  *E_alp = x_;
}
void Electronic_Structure::set_E_bet(MATRIX& x_){
  check_matrix_dimensionas(E_bet, x_, "Electronic::E_bet", "x_", "Electronic::set_E_bet");
  *E_bet = x_;
}
MATRIX Electronic_Structure::get_E_alp(){  return *E_alp; }
MATRIX Electronic_Structure::get_E_bet(){  return *E_bet; }


void Electronic_Structure::excite_alp(int I, int J){

  vector< pair<int,double> > occ_fin;
  excite(I, J, occ_alp, occ_fin);  occ_alp = occ_fin;

  // And recompute density matrix
  compute_density_matrix(occ_alp, C_alp, P_alp);

  // Also update the total density matrix:
  *P = *P_alp + *P_bet;

}

void Electronic_Structure::excite_bet(int I, int J){

  vector< pair<int,double> > occ_fin;
  excite(I, J, occ_bet, occ_fin);  occ_bet = occ_fin;

  // And recompute density matrix
  compute_density_matrix(occ_alp, C_bet, P_bet);

  // Also update the total density matrix:
  *P = *P_alp + *P_bet;

}




int init_numbers(Electronic_Structure& el, System& syst, vector<AO>& basis_ao, Model_Parameters& modprms, double charge){
/** 
  \param[in,out] el The object containing information about the electronic structure of the given system
  \param[in,out] syst The object containing information about the nuclear coordinates and atom types of the given system
  \param[in] basis_ao The vector of AO objects - the AO basis for given calculations
  \param[in,out] modprms The object that contains all the Hamiltonian parameters for given system and method choice
  \param[in] charge The total charge of the system

  Initialization of the central numbers describing electronic variables
  The function computes the number of electrons on a given set of atoms (fragment)
  Then it determines the number of alpha and beta electrons 
  Then it sets the occupations numbers of alpha and beta orbitals to 0 or 1, according to aufbau principle
*/

  int i;

  // Compute total number of valence electrons in the given fragment
  el.Nelec = 0;


  for(i=0;i<syst.Number_of_atoms;i++){

    el.Nelec += modprms.PT[syst.Atoms[i].Atom_element].Nval;  // number of valence electrons for given atom

  }// for i


  // Change the total number of electrons according to the total charge
  // !!! Fractional charge is not supported yet:
  double fractpart, intpart;
  fractpart = modf (charge , &intpart);

  if(fabs(fractpart)>0.01) { cout<<"Error in init_numbers: Fractional charge is not supported yet\n"; exit(0); }  

  el.Nelec -= intpart;


  // Numbers of alpha & beta electrons
  el.Nocc_alp = ((int)el.Nelec%2==0)?el.Nelec/2:(el.Nelec/2+1); 
  el.Nocc_bet = el.Nelec - el.Nocc_alp;


  // Setup occupation numbers for unrestricted shell:
  for(i=0;i<el.Norb;i++){
    if(i<el.Nocc_alp){ el.occ_alp.push_back(pair<int,double>(i,1.0)); }
    else{              el.occ_alp.push_back(pair<int,double>(i,0.0)); }

    if(i<el.Nocc_bet){ el.occ_bet.push_back(pair<int,double>(i,1.0)); }
    else{              el.occ_bet.push_back(pair<int,double>(i,0.0)); }
  }



}

/*
int init_electronic_subsystems(System& syst,
                               vector<vector<int> >& fragments,vector<vector<int> >& atom_to_ao_map,
                               vector<vector<int> >& basis_fo, vector<Electronic>& el,
                               vector<AO>& basis_ao, Model_Parameters& modprms, Nuclear& mol, vector<double>& frag_charge){
/// fragments[n][i] (input) - index of i-th atom belonging to n-th fragment\n
/// atom_to_ao_map[n][i] (input) - index of i-th AO belonging to n-th atom ( global array )\n
/// basis_fo[n][i] (output) - index of i-th AO belonging to n-th fragment/subsystem\n
/// el (output) - will contain Nfrag electronic structure variables\n
/// basis_ao (input) - need for extra info\n
/// modprms (input) - need for extra info\n
/// mol (input) - need for extra info\n
/// frag_charge[n] (input) - is the total charge (minus the number of excess electrons) of n-th fragment \n
/// Example: \n
/
\verbatim
    CH2=CH2 molecule with atoms:
    1  C  with orbitals 2s (0), 2px (1), 2py (2), 2pz(3)
    2  H  with orbitals 1s (4)
    3  H  with orbitals 1s (5)
    4  C  with orbitals 2s (6), 2px (7), 2py (8), 2pz (9)
    5  H  with orbitals 1s (10)
    6  H  with orbitals 1s (11)
    Can be partitioned into 2 fragments: CH2 and CH2
    so:

    fragments = [[0,1,2], [3,4,5]], so e.g. fragments[0][2] = 2 and fragments[1][1] = 4
    atom_to_ao_map = [[0,1,2,3],[4],[5],[6,7,8,9],[10],[11]]
    basis_fo = [[0,1,2,3,4,5],[6,7,8,9,10,11]]
    frag_charge = [0.0, 0.0]   

    In case if we consider cation, e.g. CH2=NH3\^{+} and want to localize the +1 charge on NH3 fragment the frag_charge array can be:
    frag_charge = [0.0, 1.0]

\endverbatim
/

  int Nfrag = fragments.size();
  if(Nfrag<1){ cout<<"Error: in init_electronic_subsystems\n Number of fragments must be more than 0\n"; exit(0);   }// default


  // In case all is ok we are here

  for(int i=0;i<Nfrag;i++){  
    vector<int> x; 

    for(int n=0;n<fragments[i].size();n++){ 
      int a = fragments[i][n]; // n-th atom in i-th fragment

      for(int j=0;j<atom_to_ao_map[a].size();j++){
        x.push_back(atom_to_ao_map[a][j]);  // j-th orbital of a-th atom
      }// for j
    }// for n

    basis_fo.push_back(x); 
    Electronic* ell; 

    ell = new Electronic(x.size());
    init_numbers(fragments[i],ell,basis_ao,modprms,mol,frag_charge[i]);


    el.push_back(ell);

  }// for i

  return Nfrag;
}
*/


}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra


