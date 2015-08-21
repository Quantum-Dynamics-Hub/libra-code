/*********************************************************************************
* Copyright (C) 2014 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
 \file Electronic.cpp
 \brief Implementation of the methods for handling electronic variables

*/
/****************************************************************************
  This file contains following functions:

  int init_numbers(vector<int>& fragment,Electronic* el,vector<AO>& basis_ao, Model_Parameters& modprms, Nuclear& mol,double charge);
  int init_electronic_subsystems(vector<vector<int> >& fragments,vector<vector<int> >& at_orbitals,
                                vector<vector<int> >& basis_fo, vector<Electronic*>& el,
                                vector<AO>& basis_ao, Model_Parameters& modprms,Nuclear&,vector<double>& frag_charge)

****************************************************************************/

#include "Electronic.h"


namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_qm{



/// Implementation of the copy constructor of the class Electronic
Electronic::Electronic(const Electronic& obj){

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

/// Implementation of the constructor from the address pointing to an existing object of class Electronic
Electronic::Electronic(Electronic* obj){

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


/// Initialization of central numbers describing electronic variables\n
/// The function computes the number of electrons on a given set of atoms (fragment)\n
/// Then it determines the number of alpha and beta electrons \n
/// Then it sets the occupations numbers of alpha and beta orbitals to 0 or 1, according to aufbau principle\n
/// charge = total charge on given fragment
int init_numbers(Electronic& el, System& syst, vector<AO>& basis_ao, Model_Parameters& modprms, double charge){

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
int init_electronic_subsystems(vector<vector<int> >& fragments,vector<vector<int> >& at_orbitals,
                               vector<vector<int> >& basis_fo, vector<Electronic*>& el,
                               vector<AO>& basis_ao, Model_Parameters& modprms, Nuclear& mol, vector<double>& frag_charge){
/// fragments[n][i] (input) - index of i-th atom belonging to n-th fragment\n
/// at_orbitals[n][i] (input) - index of i-th AO belonging to n-th atom ( global array )\n
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
    at_orbitals = [[0,1,2,3],[4],[5],[6,7,8,9],[10],[11]]
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

      for(int j=0;j<at_orbitals[a].size();j++){
        x.push_back(at_orbitals[a][j]);  // j-th orbital of a-th atom
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
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian



