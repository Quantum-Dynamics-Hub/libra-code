/*********************************************************************************
* Copyright (C) 2014-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Electronic_Structure.h
  \brief This file describes the Electronic_Structure class - a container for electronic variables.
  Some methods for handling the Electronic_Structure objects are also described here

*/

#ifndef ELECTRONIC_STRUCTURE_H
#define ELECTRONIC_STRUCTURE_H

#include "../../qobjects/libqobjects.h"
#include "../../chemobjects/libchemobjects.h"
#include "../../basis_setups/libbasis_setups.h"
#include "../../control_parameters/libcontrol_parameters.h"
#include "../../model_parameters/libmodel_parameters.h"
#include "../../calculators/libcalculators.h"


/// liblibra namespace
namespace liblibra{

using namespace libqobjects;
using namespace libchemobjects;
using namespace libchemobjects::libchemsys;
using namespace libbasis_setups;
using namespace libcontrol_parameters;
using namespace libmodel_parameters;
using namespace libcalculators;


/// libhamiltonian namespace
namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{



class Electronic_Structure{
/** This class implements the container for storing all kinds of information about electronic structure
  of an atomistic system
*/

  void check_matrix_dimensionas(MATRIX* A, MATRIX& B, std::string A_name, std::string B_name, std::string func_name);

  public:

  //-------------- Constructors --------------
  Electronic_Structure(){
  /** The default constructor - initializes all pointers to NULL-pointer. Sets the numbers to zero.
  */

    Norb  = 0;
    Nocc_alp = Nocc_bet = 0;
    Nelec = 0.0;
    P_alp = NULL; P_bet = NULL; P = NULL; 
    C_alp = NULL; C_bet = NULL;
    Sao   = NULL; Hao = NULL;
    Fao_alp = NULL; Fao_bet = NULL;
    dFao_alp_dP_alp = NULL;  dFao_alp_dP_bet = NULL;
    dFao_bet_dP_alp = NULL;  dFao_bet_dP_bet = NULL;
    E_alp = NULL; E_bet = NULL;
  }

  Electronic_Structure(int n){
  /** \param[in] n The number of orbitals (the dimensionality of the electronic problem) in given system

  The constructor with one argument - it also allocates memory for the variables
  */

    Norb = n;
    Nocc_alp = Nocc_bet = 0;
    Nelec = 0.0;
    P_alp = new MATRIX(Norb,Norb);
    P_bet = new MATRIX(Norb,Norb);
    P     = new MATRIX(Norb,Norb);
    C_alp = new MATRIX(Norb,Norb);
    C_bet = new MATRIX(Norb,Norb);
    Sao   = new MATRIX(Norb,Norb);
    Hao   = new MATRIX(Norb,Norb);
    Fao_alp = new MATRIX(Norb,Norb);
    Fao_bet = new MATRIX(Norb,Norb);
    dFao_alp_dP_alp = new MATRIX(Norb,Norb); *dFao_alp_dP_alp = 0.0;
    dFao_alp_dP_bet = new MATRIX(Norb,Norb); *dFao_alp_dP_bet = 0.0;
    dFao_bet_dP_alp = new MATRIX(Norb,Norb); *dFao_bet_dP_alp = 0.0;
    dFao_bet_dP_bet = new MATRIX(Norb,Norb); *dFao_bet_dP_bet = 0.0;
    E_alp = new MATRIX(Norb,Norb);
    E_bet = new MATRIX(Norb,Norb);   
    Mull_orb_pop_net = vector<double>(Norb,0.0);
    Mull_orb_pop_gross = vector<double>(Norb,0.0);

  }

  Electronic_Structure(const Electronic_Structure&);   ///< Copy constructor;
  Electronic_Structure(Electronic_Structure*);         ///< Constructor from the address of the other existing object

  Electronic_Structure& operator=(const Electronic_Structure&);


  
  ~Electronic_Structure(){
  /** The destructor - frees the memory for all matrices, vectors, sets them to NULL-pointer, sets the numbers to zero
  */
    delete P_alp;   delete P_bet;    delete P;        
    delete C_alp;   delete C_bet;
    delete Sao;     delete Hao;
    delete Fao_alp; delete Fao_bet;
    delete dFao_alp_dP_alp;  delete dFao_alp_dP_bet;
    delete dFao_bet_dP_alp;  delete dFao_bet_dP_bet;
    delete E_alp;   delete E_bet;

    Mull_orb_pop_net.clear();
    Mull_orb_pop_gross.clear();

    if(bands_alp.size()>0){ bands_alp.clear(); }
    if(bands_bet.size()>0){ bands_bet.clear(); }
    if(occ_alp.size()>0){ occ_alp.clear(); }
    if(occ_bet.size()>0){ occ_bet.clear(); }
   
    Norb  = 0;
    Nocc_alp = Nocc_bet = 0;
    Nelec = 0.0;
    P_alp = NULL; P_bet = NULL; P = NULL; 
    C_alp = NULL; C_bet = NULL;
    Sao   = NULL; Hao = NULL;
    Fao_alp = NULL; Fao_bet = NULL;
    dFao_alp_dP_alp = NULL;  dFao_alp_dP_bet = NULL;
    dFao_bet_dP_alp = NULL;  dFao_bet_dP_bet = NULL;
    E_alp = NULL; E_bet = NULL;

  }


 
  //-------------- Data members --------------
  int Norb;        ///< the total number of orbitals in this subsystem
  int Nocc_alp;    ///< the number of occupied alpha orbitals in this subsystem
  int Nocc_bet;    ///< the number of occupied beta orbitals in this subsystem
  double Nelec;    ///< the number of electrons in this subsystem

  vector< pair<int,double> > bands_alp; ///< orbital indices and orbital energies, alpha-channel
  vector< pair<int,double> > bands_bet; ///< orbital indices and orbital energies, beta-channel
  vector< pair<int,double> > occ_alp;   ///< orbital indices and orbital occupation numbers, alpha-channel
  vector< pair<int,double> > occ_bet;   ///< orbital indices and orbital occupation numbers, beta-channel
  double get_bands_alp(int indx){ return bands_alp[indx].second; }
  double get_bands_bet(int indx){ return bands_bet[indx].second; }
  double get_occ_alp(int indx){ return occ_alp[indx].second; }
  double get_occ_bet(int indx){ return occ_bet[indx].second; }
  

  // Density matrices
  MATRIX* P_alp;                        ///< Density matrix, alpha-channel
  MATRIX* P_bet;                        ///< Density matrix, beta-channel
  MATRIX* P;                            ///< Density matrix, total
  void set_P_alp(MATRIX& x_);
  void set_P_bet(MATRIX& x_);
  void set_P(MATRIX& x_);
  MATRIX get_P_alp();
  MATRIX get_P_bet();
  MATRIX get_P();


  // Wfc coefficients
  MATRIX* C_alp;                        ///< MO-LCAO coefficients, alpha-channel
                                        ///< C_alp[k][i] - is the i-th MO (column) k-th AO (row)
  MATRIX* C_bet;                        ///< MO-LCAO coefficients, beta-channel                     
                                        ///< C_bet[k][i] - is the i-th MO (column) k-th AO (row)
  void set_C_alp(MATRIX& x_);
  void set_C_bet(MATRIX& x_);
  MATRIX get_C_alp();
  MATRIX get_C_bet();

  // Overlaps
  MATRIX* Sao;                          ///< Overlap in AO basis
  void set_Sao(MATRIX& x_);
  MATRIX get_Sao();

  // Core Hamiltonian
  MATRIX* Hao;                          ///< Core Hamiltonian in AO basis
  void set_Hao(MATRIX& x_);
  MATRIX get_Hao();

  // Fock matrices
  MATRIX* Fao_alp;                      ///< Fock matrix in AO basis, alpha-channel
  MATRIX* Fao_bet;                      ///< Fock matrix in AO basis, beta-channel
  void set_Fao_alp(MATRIX& x_);
  void set_Fao_bet(MATRIX& x_);
  MATRIX get_Fao_alp();
  MATRIX get_Fao_bet();


  // Corrections:
  MATRIX* dFao_alp_dP_alp;              ///< derivative of the Fock matrix (alpha component) w.r.t. alpha component of density matrix
  MATRIX* dFao_alp_dP_bet;              ///< derivative of the Fock matrix (alpha component) w.r.t. beta component of density matrix
  MATRIX* dFao_bet_dP_alp;              ///< derivative of the Fock matrix (beta component) w.r.t. alpha component of density matrix
  MATRIX* dFao_bet_dP_bet;              ///< derivative of the Fock matrix (beta component) w.r.t. beta component of density matrix
  void set_dFao_alp_dP_alp(MATRIX& x_);
  void set_dFao_alp_dP_bet(MATRIX& x_);
  void set_dFao_bet_dP_alp(MATRIX& x_);
  void set_dFao_bet_dP_bet(MATRIX& x_);
  MATRIX get_dFao_alp_dP_alp();
  MATRIX get_dFao_alp_dP_bet();
  MATRIX get_dFao_bet_dP_alp();
  MATRIX get_dFao_bet_dP_bet();


  // Eigenvalues
  MATRIX* E_alp;                        ///< MO energies, alpha-channel
  MATRIX* E_bet;                        ///< MO energies, beta-channel
  void set_E_alp(MATRIX& x_);
  void set_E_bet(MATRIX& x_);
  MATRIX get_E_alp();
  MATRIX get_E_bet();

  vector<double> Mull_orb_pop_net;      ///< Net Mulliken populations for all (molecular) orbitals
  vector<double> Mull_orb_pop_gross;    ///< Gross Mulliken populations for all (molecular) orbitals



  void excite_alp(int I,int J);
  void excite_bet(int I,int J);

};

// Electronic.cpp
int init_numbers(vector<int>&,Electronic_Structure*, vector<AO>&, Model_Parameters&, System&, double);

//int init_electronic_subsystems(vector<vector<int> >&,vector<vector<int> >&, vector<vector<int> >&, vector<Electronic*>&,
//                               vector<AO>&,Model_Parameters&, Nuclear&,vector<double>&);


/*
// NAC.cpp
void compute_overlap_nac(Electronic*, Electronic*, int, int, vector<int>&,
                         Electronic*, Electronic*, int, int, vector<int>&,
                         vector<AO>&,vector<AO>&,
                         double, vector<double*>&, int,MATRIX*, MATRIX*, MATRIX*, int);

void collect_matrices(int, vector<MATRIX*>&, vector<int>&, vector<int>&,
                      MATRIX*, vector<int>&,vector<int>&);
*/


}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra


#endif // ELECTRONIC_STRUCTURE_H

