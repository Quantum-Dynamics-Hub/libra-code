/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Hamiltonian_QM.cpp
  \brief The file implements functions for quantum-mechanical Hamiltonian calculations as well
  as the classes for organizing such computations in an object-oriented way. 
*/

#include "Hamiltonian_QM.h"
#include "SCF.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{



void Hamiltonian_core(
  System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF
){
/**
  \param[in,out] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the ab initio calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[out] Hao The Atomistic Hamiltonian we compute
  \param[in] Sao The AO overlap matrix
  \param[in] DF The debug flag - controlls the amount of additionally printed information

  The generic function for computing the core Hamiltonian for a given system
*/

  if(prms.hamiltonian=="hf"){

    Hamiltonian_core_hf(syst, basis_ao,  prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, DF);

  }
  else if(prms.hamiltonian=="indo"){

    Hamiltonian_core_indo(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, DF);

  }
  else if(prms.hamiltonian=="eht"){

    Hamiltonian_core_eht(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, DF);

  }


}


void Hamiltonian_Fock(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
                      Control_Parameters& prms,Model_Parameters& modprms,
                      vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                     ){
/**
  \param[in,out] el The object containing all information about electronic structure of the system, including the computed Fock matrix
  \param[in,out] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the ab initio calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized

  The generic function for computing the Fock Hamiltonian for a given system
*/


  if(prms.hamiltonian=="hf"){

    Hamiltonian_Fock_hf(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

  }
  else if(prms.hamiltonian=="indo"){

    Hamiltonian_Fock_indo(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

  }
  else if(prms.hamiltonian=="eht"){

    Hamiltonian_Fock_eht(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

  }



}

void Hamiltonian_Fock(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
                      Control_Parameters& prms,Model_Parameters& modprms,
                      vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                     ){
  Hamiltonian_Fock(&el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

}



void derivative_couplings
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
  int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3,
  MATRIX& Dmo_a_x, MATRIX& Dmo_a_y, MATRIX& Dmo_a_z,
  MATRIX& Dmo_b_x, MATRIX& Dmo_b_y, MATRIX& Dmo_b_z
){
/**
  \param[in] el The object containing all information about electronic structure of the system
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the ab initio calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] Hao The core Hamiltonian matrix
  \param[in] Sao The AO overlap matrix
  \param[in] Norb The number of orbitals in the system
  \param[in] at_indx The index of the atom w.r.t. whose coordinates we are computing the derivative couplings
  \param[in] x_period Then number of periodic shells in X direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] y_period Then number of periodic shells in Y direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] z_period Then number of periodic shells in Z direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] t1 The periodicity vector along a crystal direction ("X")
  \param[in] t2 The periodicity vector along b crystal direction ("Y")
  \param[in] t3 The periodicity vector along c crystal direction ("Z")
  \param[out] Dmo_a_x The derivative coupling of the alpha-orbitals w.r.t. x-projection of the selected atom 
  \param[out] Dmo_a_y The derivative coupling of the alpha-orbitals w.r.t. y-projection of the selected atom 
  \param[out] Dmo_a_z The derivative coupling of the alpha-orbitals w.r.t. z-projection of the selected atom 
  \param[out] Dmo_b_x The derivative coupling of the beta-orbitals w.r.t. x-projection of the selected atom 
  \param[out] Dmo_b_y The derivative coupling of the beta-orbitals w.r.t. y-projection of the selected atom 
  \param[out] Dmo_b_z The derivative coupling of the beta-orbitals w.r.t. z-projection of the selected atom 


  The function for computing the derivative couplings for a given system

  WARNING: Use this version with caution - need more testing!
*/




  int i,j;
    
  MATRIX* dHao_dx; dHao_dx = new MATRIX(Norb, Norb);
  MATRIX* dHao_dy; dHao_dy = new MATRIX(Norb, Norb);
  MATRIX* dHao_dz; dHao_dz = new MATRIX(Norb, Norb);

  MATRIX* dSao_dx; dSao_dx = new MATRIX(Norb, Norb);
  MATRIX* dSao_dy; dSao_dy = new MATRIX(Norb, Norb);
  MATRIX* dSao_dz; dSao_dz = new MATRIX(Norb, Norb);

  MATRIX* dFao_alp_dx; dFao_alp_dx = new MATRIX(Norb, Norb);
  MATRIX* dFao_alp_dy; dFao_alp_dy = new MATRIX(Norb, Norb);
  MATRIX* dFao_alp_dz; dFao_alp_dz = new MATRIX(Norb, Norb);

  MATRIX* dFao_bet_dx; dFao_bet_dx = new MATRIX(Norb, Norb);
  MATRIX* dFao_bet_dy; dFao_bet_dy = new MATRIX(Norb, Norb);
  MATRIX* dFao_bet_dz; dFao_bet_dz = new MATRIX(Norb, Norb);

  MATRIX* Dao_x; Dao_x = new MATRIX(Norb, Norb);
  MATRIX* Dao_y; Dao_y = new MATRIX(Norb, Norb);
  MATRIX* Dao_z; Dao_z = new MATRIX(Norb, Norb);

  int DF = 0;

  if(prms.hamiltonian=="indo"){

    Hamiltonian_core_deriv_indo(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, DF, at_indx, *dHao_dx, *dHao_dy, *dHao_dz, *dSao_dx, *dSao_dy, *dSao_dz );
    Hamiltonian_Fock_derivs_indo(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, at_indx, *dHao_dx, *dHao_dy, *dHao_dz, *dFao_alp_dx, *dFao_alp_dy, *dFao_alp_dz, *dFao_bet_dx, *dFao_bet_dy, *dFao_bet_dz);

  }

  update_derivative_coupling_matrix(x_period, y_period, z_period, t1, t2, t3, atom_to_ao_map, ao_to_atom_map, basis_ao, at_indx, *Dao_x, *Dao_y, *Dao_z);


  // !!! Because in INDO  S = I => Dao = 0.0
  *dSao_dx = 0.0;
  *dSao_dy = 0.0;
  *dSao_dz = 0.0;
  *Dao_x = 0.0;
  *Dao_y = 0.0;
  *Dao_z = 0.0;


  MATRIX* T; T = new MATRIX(Norb,Norb);
  MATRIX* A_ax; A_ax = new MATRIX(Norb,Norb);
  MATRIX* A_ay; A_ay = new MATRIX(Norb,Norb);
  MATRIX* A_az; A_az = new MATRIX(Norb,Norb);
  MATRIX* A_bx; A_bx = new MATRIX(Norb,Norb);
  MATRIX* A_by; A_by = new MATRIX(Norb,Norb);
  MATRIX* A_bz; A_bz = new MATRIX(Norb,Norb);


  *T = (*el.E_alp) * (*el.C_alp) * (*Dao_x) * (*el.C_alp).T();
  *A_ax = (*el.C_alp) * (*dFao_alp_dx) * (*el.C_alp).T() - (*T + (*T).T());

  *T = (*el.E_alp) * (*el.C_alp) * (*Dao_y) * (*el.C_alp).T();
  *A_ay = (*el.C_alp) * (*dFao_alp_dy) * (*el.C_alp).T() - (*T + (*T).T());

  *T = (*el.E_alp) * (*el.C_alp) * (*Dao_z) * (*el.C_alp).T();
  *A_az = (*el.C_alp) * (*dFao_alp_dz) * (*el.C_alp).T() - (*T + (*T).T());

  *T = (*el.E_bet) * (*el.C_bet) * (*Dao_x) * (*el.C_bet).T();
  *A_bx = (*el.C_bet) * (*dFao_bet_dx) * (*el.C_bet).T() - (*T + (*T).T());

  *T = (*el.E_bet) * (*el.C_bet) * (*Dao_y) * (*el.C_bet).T();
  *A_by = (*el.C_bet) * (*dFao_bet_dy) * (*el.C_bet).T() - (*T + (*T).T());

  *T = (*el.E_bet) * (*el.C_bet) * (*Dao_z) * (*el.C_bet).T();
  *A_bz = (*el.C_bet) * (*dFao_bet_dz) * (*el.C_bet).T() - (*T + (*T).T());


  MATRIX* dEa_dx; dEa_dx = new MATRIX(Norb,Norb);
  MATRIX* dEa_dy; dEa_dy = new MATRIX(Norb,Norb);
  MATRIX* dEa_dz; dEa_dz = new MATRIX(Norb,Norb);
  MATRIX* dEb_dx; dEb_dx = new MATRIX(Norb,Norb);
  MATRIX* dEb_dy; dEb_dy = new MATRIX(Norb,Norb);
  MATRIX* dEb_dz; dEb_dz = new MATRIX(Norb,Norb);


  for(i=0;i<Norb;i++){
    dEa_dx->set(i,i,A_ax->get(i,i));
    dEa_dy->set(i,i,A_ay->get(i,i));
    dEa_dz->set(i,i,A_az->get(i,i));
    dEb_dx->set(i,i,A_bx->get(i,i));
    dEb_dy->set(i,i,A_by->get(i,i));
    dEb_dz->set(i,i,A_bz->get(i,i));

    for(j=0;j<Norb;j++){
      if(i!=j){
        double dEa = ( el.E_alp->get(j,j)-el.E_alp->get(i,i) );
        if(fabs(dEa)>1e-12){
          Dmo_a_x.set(i,j,A_ax->get(i,j)/dEa);
          Dmo_a_y.set(i,j,A_ay->get(i,j)/dEa);
          Dmo_a_z.set(i,j,A_az->get(i,j)/dEa);
        }

        double dEb = ( el.E_bet->get(j,j)-el.E_bet->get(i,i) );
        if(fabs(dEb)>1e-12){
          Dmo_b_x.set(i,j,A_bx->get(i,j)/dEb);
          Dmo_b_y.set(i,j,A_by->get(i,j)/dEb);
          Dmo_b_z.set(i,j,A_bz->get(i,j)/dEb);
        }
      }// i!=j
      else{
          Dmo_a_x.set(i,j,0.0);
          Dmo_a_y.set(i,j,0.0);
          Dmo_a_z.set(i,j,0.0);
          Dmo_b_x.set(i,j,0.0);
          Dmo_b_y.set(i,j,0.0);
          Dmo_b_z.set(i,j,0.0);
      }
    }// for j
  }// for i


  //================= Clean up - delete temporary memory blocks ==============
  delete dHao_dx;  delete dHao_dy;  delete dHao_dz;
  delete dSao_dx;  delete dSao_dy;  delete dSao_dz;
  delete dFao_alp_dx;  delete dFao_alp_dy;  delete dFao_alp_dz;
  delete dFao_bet_dx;  delete dFao_bet_dy;  delete dFao_bet_dz;
  delete Dao_x;  delete Dao_y;  delete Dao_z;  delete T;
  delete A_ax;  delete A_ay;  delete A_az;
  delete A_bx;  delete A_by;  delete A_bz;
  delete dEa_dx;  delete dEa_dy;  delete dEa_dz;
  delete dEb_dx;  delete dEb_dy;  delete dEb_dz;


}



void derivative_couplings1
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
  int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3,
  MATRIX& Dmo_a_x, MATRIX& Dmo_a_y, MATRIX& Dmo_a_z,
  MATRIX& Dmo_b_x, MATRIX& Dmo_b_y, MATRIX& Dmo_b_z,
  MATRIX& dEa_dx,  MATRIX& dEa_dy,  MATRIX& dEa_dz,
  MATRIX& dEb_dx,  MATRIX& dEb_dy,  MATRIX& dEb_dz
){
/**
  \param[in] el The object containing all information about electronic structure of the system
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the ab initio calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] Hao The core Hamiltonian matrix
  \param[in] Sao The AO overlap matrix
  \param[in] Norb The number of orbitals in the system
  \param[in] at_indx The index of the atom w.r.t. whose coordinates we are computing the derivative couplings
  \param[in] x_period Then number of periodic shells in X direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] y_period Then number of periodic shells in Y direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] z_period Then number of periodic shells in Z direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] t1 The periodicity vector along a crystal direction ("X")
  \param[in] t2 The periodicity vector along b crystal direction ("Y")
  \param[in] t3 The periodicity vector along c crystal direction ("Z")
  \param[out] Dmo_a_x The derivative coupling of the alpha-orbitals w.r.t. x-projection of the selected atom 
  \param[out] Dmo_a_y The derivative coupling of the alpha-orbitals w.r.t. y-projection of the selected atom 
  \param[out] Dmo_a_z The derivative coupling of the alpha-orbitals w.r.t. z-projection of the selected atom 
  \param[out] Dmo_b_x The derivative coupling of the beta-orbitals w.r.t. x-projection of the selected atom 
  \param[out] Dmo_b_y The derivative coupling of the beta-orbitals w.r.t. y-projection of the selected atom 
  \param[out] Dmo_b_z The derivative coupling of the beta-orbitals w.r.t. z-projection of the selected atom 
  \param[out] dEa_dx The derivative of the alpha eigenvalues w.r.t. x-projection of the selected atom
  \param[out] dEa_dy The derivative of the alpha eigenvalues w.r.t. y-projection of the selected atom
  \param[out] dEa_dz The derivative of the alpha eigenvalues w.r.t. z-projection of the selected atom
  \param[out] dEb_dx The derivative of the beta eigenvalues w.r.t. x-projection of the selected atom
  \param[out] dEb_dy The derivative of the beta eigenvalues w.r.t. y-projection of the selected atom
  \param[out] dEb_dz The derivative of the beta eigenvalues w.r.t. z-projection of the selected atom


  The function computes the derivative couplings and the derivatives of the eigenvalues w.r.t. given nuclear coordinates

  This is simpler version of the derivative couplings:
  From the derivation of Hellman-Feynman theorem:

  <i| dH/dR |j> = dE_i/dR * delta_ij - (E_i - E_j)*D_ij

  where D_ij = <i| dj/dR > - derivative coupling in MO basis

  <i| dH/dR |j> - is also in MO basis, but can be transformed from
  the AO basis

*/


  int i,j;
    
  MATRIX* dHao_dx; dHao_dx = new MATRIX(Norb, Norb);
  MATRIX* dHao_dy; dHao_dy = new MATRIX(Norb, Norb);
  MATRIX* dHao_dz; dHao_dz = new MATRIX(Norb, Norb);

  MATRIX* dSao_dx; dSao_dx = new MATRIX(Norb, Norb);
  MATRIX* dSao_dy; dSao_dy = new MATRIX(Norb, Norb);
  MATRIX* dSao_dz; dSao_dz = new MATRIX(Norb, Norb);

  MATRIX* dFao_alp_dx; dFao_alp_dx = new MATRIX(Norb, Norb);
  MATRIX* dFao_alp_dy; dFao_alp_dy = new MATRIX(Norb, Norb);
  MATRIX* dFao_alp_dz; dFao_alp_dz = new MATRIX(Norb, Norb);

  MATRIX* dFao_bet_dx; dFao_bet_dx = new MATRIX(Norb, Norb);
  MATRIX* dFao_bet_dy; dFao_bet_dy = new MATRIX(Norb, Norb);
  MATRIX* dFao_bet_dz; dFao_bet_dz = new MATRIX(Norb, Norb);


  int DF = 0;

  if(prms.hamiltonian=="indo"){

    Hamiltonian_core_deriv_indo(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, DF, at_indx, *dHao_dx, *dHao_dy, *dHao_dz, *dSao_dx, *dSao_dy, *dSao_dz );
    Hamiltonian_Fock_derivs_indo(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, at_indx, *dHao_dx, *dHao_dy, *dHao_dz, *dFao_alp_dx, *dFao_alp_dy, *dFao_alp_dz, *dFao_bet_dx, *dFao_bet_dy, *dFao_bet_dz);

  }

// Don't need this
//  update_derivative_coupling_matrix(x_period, y_period, z_period, t1, t2, t3, atom_to_ao_map, ao_to_atom_map, basis_ao, at_indx, *Dao_x, *Dao_y, *Dao_z);


  MATRIX* A_ax; A_ax = new MATRIX(Norb,Norb);
  MATRIX* A_ay; A_ay = new MATRIX(Norb,Norb);
  MATRIX* A_az; A_az = new MATRIX(Norb,Norb);
  MATRIX* A_bx; A_bx = new MATRIX(Norb,Norb);
  MATRIX* A_by; A_by = new MATRIX(Norb,Norb);
  MATRIX* A_bz; A_bz = new MATRIX(Norb,Norb);

  // From AO basis to MO basis
  *A_ax = (*el.C_alp).T() * (*dFao_alp_dx) * (*el.C_alp);
  *A_ay = (*el.C_alp).T() * (*dFao_alp_dy) * (*el.C_alp);
  *A_az = (*el.C_alp).T() * (*dFao_alp_dz) * (*el.C_alp);
  *A_bx = (*el.C_bet).T() * (*dFao_bet_dx) * (*el.C_bet);
  *A_by = (*el.C_bet).T() * (*dFao_bet_dy) * (*el.C_bet);
  *A_bz = (*el.C_bet).T() * (*dFao_bet_dz) * (*el.C_bet);


  for(i=0;i<Norb;i++){
    dEa_dx.set(i,i,A_ax->get(i,i));
    dEa_dy.set(i,i,A_ay->get(i,i));
    dEa_dz.set(i,i,A_az->get(i,i));
    dEb_dx.set(i,i,A_bx->get(i,i));
    dEb_dy.set(i,i,A_by->get(i,i));
    dEb_dz.set(i,i,A_bz->get(i,i));

    for(j=0;j<Norb;j++){
      if(i!=j){
        double dEa = ( el.E_alp->get(j,j)-el.E_alp->get(i,i) );
        if(fabs(dEa)>1e-12){
          double nac;
          nac = A_ax->get(i,j)/dEa; if(fabs(nac)>1.0){ nac /= fabs(nac); }     Dmo_a_x.set(i,j,nac);
          nac = A_ay->get(i,j)/dEa; if(fabs(nac)>1.0){ nac /= fabs(nac); }     Dmo_a_y.set(i,j,nac);
          nac = A_az->get(i,j)/dEa; if(fabs(nac)>1.0){ nac /= fabs(nac); }     Dmo_a_z.set(i,j,nac);
        }

        double dEb = ( el.E_bet->get(j,j)-el.E_bet->get(i,i) );
        if(fabs(dEb)>1e-12){
          double nac;
          nac = A_bx->get(i,j)/dEb; if(fabs(nac)>1.0){ nac /= fabs(nac); }     Dmo_b_x.set(i,j,nac);
          nac = A_by->get(i,j)/dEb; if(fabs(nac)>1.0){ nac /= fabs(nac); }     Dmo_b_y.set(i,j,nac);
          nac = A_bz->get(i,j)/dEb; if(fabs(nac)>1.0){ nac /= fabs(nac); }     Dmo_b_z.set(i,j,nac);
        }
      }// i!=j
      else{ // i==j - diagonal terms
          Dmo_a_x.set(i,j,0.0);
          Dmo_a_y.set(i,j,0.0);
          Dmo_a_z.set(i,j,0.0);
          Dmo_b_x.set(i,j,0.0);
          Dmo_b_y.set(i,j,0.0);
          Dmo_b_z.set(i,j,0.0);
      }
    }// for j
  }// for i


  //================= Clean up - delete temporary memory blocks ==============
  delete dHao_dx;  delete dHao_dy;  delete dHao_dz;
  delete dSao_dx;  delete dSao_dy;  delete dSao_dz;
  delete dFao_alp_dx;  delete dFao_alp_dy;  delete dFao_alp_dz;
  delete dFao_bet_dx;  delete dFao_bet_dy;  delete dFao_bet_dz;
  delete A_ax;  delete A_ay;  delete A_az;
  delete A_bx;  delete A_by;  delete A_bz;


}

void derivative_couplings1
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
  int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3,
  MATRIX& Dmo_a_x, MATRIX& Dmo_a_y, MATRIX& Dmo_a_z,
  MATRIX& Dmo_b_x, MATRIX& Dmo_b_y, MATRIX& Dmo_b_z
){
/**
  \param[in] el The object containing all information about electronic structure of the system
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the ab initio calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] Hao The core Hamiltonian matrix
  \param[in] Sao The AO overlap matrix
  \param[in] Norb The number of orbitals in the system
  \param[in] at_indx The index of the atom w.r.t. whose coordinates we are computing the derivative couplings
  \param[in] x_period Then number of periodic shells in X direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] y_period Then number of periodic shells in Y direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] z_period Then number of periodic shells in Z direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] t1 The periodicity vector along a crystal direction ("X")
  \param[in] t2 The periodicity vector along b crystal direction ("Y")
  \param[in] t3 The periodicity vector along c crystal direction ("Z")
  \param[out] Dmo_a_x The derivative coupling of the alpha-orbitals w.r.t. x-projection of the selected atom 
  \param[out] Dmo_a_y The derivative coupling of the alpha-orbitals w.r.t. y-projection of the selected atom 
  \param[out] Dmo_a_z The derivative coupling of the alpha-orbitals w.r.t. z-projection of the selected atom 
  \param[out] Dmo_b_x The derivative coupling of the beta-orbitals w.r.t. x-projection of the selected atom 
  \param[out] Dmo_b_y The derivative coupling of the beta-orbitals w.r.t. y-projection of the selected atom 
  \param[out] Dmo_b_z The derivative coupling of the beta-orbitals w.r.t. z-projection of the selected atom 


  The function computes the derivative couplings 

  This is simpler version of the derivative couplings:
  From the derivation of Hellman-Feynman theorem:

  <i| dH/dR |j> = dE_i/dR * delta_ij - (E_i - E_j)*D_ij

  where D_ij = <i| dj/dR > - derivative coupling in MO basis

  <i| dH/dR |j> - is also in MO basis, but can be transformed from
  the AO basis

*/


  MATRIX* dEa_dx; dEa_dx = new MATRIX(Norb,Norb);
  MATRIX* dEa_dy; dEa_dy = new MATRIX(Norb,Norb);
  MATRIX* dEa_dz; dEa_dz = new MATRIX(Norb,Norb);
  MATRIX* dEb_dx; dEb_dx = new MATRIX(Norb,Norb);
  MATRIX* dEb_dy; dEb_dy = new MATRIX(Norb,Norb);
  MATRIX* dEb_dz; dEb_dz = new MATRIX(Norb,Norb);

  derivative_couplings1(el, syst, basis_ao,
  prms,modprms,atom_to_ao_map,ao_to_atom_map,
  Hao, Sao, Norb, at_indx, x_period, y_period, z_period, t1, t2, t3,
  Dmo_a_x, Dmo_a_y, Dmo_a_z,  Dmo_b_x, Dmo_b_y, Dmo_b_z,
  *dEa_dx,  *dEa_dy,  *dEa_dz,  *dEb_dx,  *dEb_dy,  *dEb_dz);


  delete dEa_dx;  delete dEa_dy;  delete dEa_dz;
  delete dEb_dx;  delete dEb_dy;  delete dEb_dz;


}




VECTOR force
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
  int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3){
/**
  \param[in] el The object containing all information about electronic structure of the system
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the ab initio calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] Hao The core Hamiltonian matrix 
  \param[in] Sao The AO overlap matrix
  \param[in] Norb The number of orbitals in the system
  \param[in] at_indx The index of the atom w.r.t. whose coordinates we are computing the derivative couplings
  \param[in] x_period Then number of periodic shells in X direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] y_period Then number of periodic shells in Y direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] z_period Then number of periodic shells in Z direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] t1 The periodicity vector along a crystal direction ("X")
  \param[in] t2 The periodicity vector along b crystal direction ("Y")
  \param[in] t3 The periodicity vector along c crystal direction ("Z")

  Compute the force (corresponding to the electronic structure given by el object) acting on the selected (with index
  at_indx ) atom. Return the force vector. This is only electronic energy contribution

  This version of the force computation is good for semiempirical methods with identity overlap matrix (e.g. INDO, but not EHT)
*/


//  cout<<"in force (Hamiltonian_QM.cpp)\n";
  Timer tim1;

  int i,j;
    
  MATRIX* dHao_dx; dHao_dx = new MATRIX(Norb, Norb);
  MATRIX* dHao_dy; dHao_dy = new MATRIX(Norb, Norb);
  MATRIX* dHao_dz; dHao_dz = new MATRIX(Norb, Norb);

  MATRIX* dSao_dx; dSao_dx = new MATRIX(Norb, Norb);
  MATRIX* dSao_dy; dSao_dy = new MATRIX(Norb, Norb);
  MATRIX* dSao_dz; dSao_dz = new MATRIX(Norb, Norb);

  MATRIX* dFao_alp_dx; dFao_alp_dx = new MATRIX(Norb, Norb);
  MATRIX* dFao_alp_dy; dFao_alp_dy = new MATRIX(Norb, Norb);
  MATRIX* dFao_alp_dz; dFao_alp_dz = new MATRIX(Norb, Norb);

  MATRIX* dFao_bet_dx; dFao_bet_dx = new MATRIX(Norb, Norb);
  MATRIX* dFao_bet_dy; dFao_bet_dy = new MATRIX(Norb, Norb);
  MATRIX* dFao_bet_dz; dFao_bet_dz = new MATRIX(Norb, Norb);

  int DF = 0;

  if(prms.hamiltonian=="indo"){

    Hamiltonian_core_deriv_indo(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, DF, at_indx, *dHao_dx, *dHao_dy, *dHao_dz, *dSao_dx, *dSao_dy, *dSao_dz );
    Hamiltonian_Fock_derivs_indo(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, at_indx, *dHao_dx, *dHao_dy, *dHao_dz, *dFao_alp_dx, *dFao_alp_dy, *dFao_alp_dz, *dFao_bet_dx, *dFao_bet_dy, *dFao_bet_dz);

  }


  VECTOR F; F = 0.0;
  // But this one seems to work pretty well - at least for INDO!
  F.x  = (  (*el.P_alp) * (*dHao_dx + *dFao_alp_dx)  ).tr();
  F.x += (  (*el.P_bet) * (*dHao_dx + *dFao_bet_dx)  ).tr();
  F.x  = -0.5 * F.x;

  F.y  = (  (*el.P_alp) * (*dHao_dy + *dFao_alp_dy)  ).tr();
  F.y += (  (*el.P_bet) * (*dHao_dy + *dFao_bet_dy)  ).tr();
  F.y  = -0.5 * F.y;

  F.z  = (  (*el.P_alp) * (*dHao_dz + *dFao_alp_dz)  ).tr();
  F.z += (  (*el.P_bet) * (*dHao_dz + *dFao_bet_dz)  ).tr();
  F.z  = -0.5 * F.z;



  //================= Clean up - delete temporary memory blocks ==============
  delete dHao_dx;  delete dHao_dy;  delete dHao_dz;
  delete dSao_dx;  delete dSao_dy;  delete dSao_dz;
  delete dFao_alp_dx;  delete dFao_alp_dy;  delete dFao_alp_dz;
  delete dFao_bet_dx;  delete dFao_bet_dy;  delete dFao_bet_dz;


  return F;

}



VECTOR force_extended
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
  int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3){
/**
  \param[in] el The object containing all information about electronic structure of the system
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the ab initio calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] Hao The core Hamiltonian matrix 
  \param[in] Sao The AO overlap matrix
  \param[in] Norb The number of orbitals in the system
  \param[in] at_indx The index of the atom w.r.t. whose coordinates we are computing the derivative couplings
  \param[in] x_period Then number of periodic shells in X direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] y_period Then number of periodic shells in Y direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] z_period Then number of periodic shells in Z direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] t1 The periodicity vector along a crystal direction ("X")
  \param[in] t2 The periodicity vector along b crystal direction ("Y")
  \param[in] t3 The periodicity vector along c crystal direction ("Z")

  Compute the force (corresponding to the electronic structure given by el object) acting on the selected (with index
  at_indx ) atom. Return the force vector. This is only electronic energy contribution

  This is the extended version of forces computations - just in case we need 
  more general formulation - but at present (for INDO) it is not used

*/




  int i,j;
    
  MATRIX* dHao_dx; dHao_dx = new MATRIX(Norb, Norb);
  MATRIX* dHao_dy; dHao_dy = new MATRIX(Norb, Norb);
  MATRIX* dHao_dz; dHao_dz = new MATRIX(Norb, Norb);

  MATRIX* dSao_dx; dSao_dx = new MATRIX(Norb, Norb);
  MATRIX* dSao_dy; dSao_dy = new MATRIX(Norb, Norb);
  MATRIX* dSao_dz; dSao_dz = new MATRIX(Norb, Norb);

  MATRIX* dFao_alp_dx; dFao_alp_dx = new MATRIX(Norb, Norb);
  MATRIX* dFao_alp_dy; dFao_alp_dy = new MATRIX(Norb, Norb);
  MATRIX* dFao_alp_dz; dFao_alp_dz = new MATRIX(Norb, Norb);

  MATRIX* dFao_bet_dx; dFao_bet_dx = new MATRIX(Norb, Norb);
  MATRIX* dFao_bet_dy; dFao_bet_dy = new MATRIX(Norb, Norb);
  MATRIX* dFao_bet_dz; dFao_bet_dz = new MATRIX(Norb, Norb);

  MATRIX* Dao_x; Dao_x = new MATRIX(Norb, Norb);
  MATRIX* Dao_y; Dao_y = new MATRIX(Norb, Norb);
  MATRIX* Dao_z; Dao_z = new MATRIX(Norb, Norb);

  int DF = 0;


  if(prms.hamiltonian=="indo"){

    Hamiltonian_core_deriv_indo(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, DF, at_indx, *dHao_dx, *dHao_dy, *dHao_dz, *dSao_dx, *dSao_dy, *dSao_dz );
    Hamiltonian_Fock_derivs_indo(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, at_indx, *dHao_dx, *dHao_dy, *dHao_dz, *dFao_alp_dx, *dFao_alp_dy, *dFao_alp_dz, *dFao_bet_dx, *dFao_bet_dy, *dFao_bet_dz);

  }

  update_derivative_coupling_matrix(x_period, y_period, z_period, t1, t2, t3, atom_to_ao_map, ao_to_atom_map, basis_ao, at_indx, *Dao_x, *Dao_y, *Dao_z);

  // For debug
  /*
  cout<<"dHao_dx = \n";    dHao_dx->show_matrix();
  cout<<"dHao_dy = \n";    dHao_dy->show_matrix();
  cout<<"dHao_dz = \n";    dHao_dz->show_matrix();

  cout<<"dFao_alp_dx = \n";    dFao_alp_dx->show_matrix();
  cout<<"dFao_alp_dy = \n";    dFao_alp_dy->show_matrix();
  cout<<"dFao_alp_dz = \n";    dFao_alp_dz->show_matrix();

  cout<<"Dao_x = \n";    Dao_x.show_matrix();
  cout<<"Dao_y = \n";    Dao_y.show_matrix();
  cout<<"Dao_z = \n";    Dao_z.show_matrix();
  cout<<"Checking properties of Dao_x vs. dSao_dx\n";
  (( Dao_x + Dao_x.T() ) - dSao_dx).show_matrix();
  cout<<"Checking properties of Dao_y vs. dSao_dy\n";
  (( Dao_y + Dao_y.T() ) - dSao_dy).show_matrix();
  cout<<"Checking properties of Dao_z vs. dSao_dz\n";
  (( Dao_z + Dao_z.T() ) - dSao_dz).show_matrix();
  */

  // !!! Because in INDO  S = I => Dao = 0.0
  *dSao_dx = 0.0;
  *dSao_dy = 0.0;
  *dSao_dz = 0.0;
  *Dao_x = 0.0;
  *Dao_y = 0.0;
  *Dao_z = 0.0;


  MATRIX* T; T = new MATRIX(Norb,Norb);
  MATRIX* A_ax; A_ax = new MATRIX(Norb,Norb);
  MATRIX* A_ay; A_ay = new MATRIX(Norb,Norb);
  MATRIX* A_az; A_az = new MATRIX(Norb,Norb);
  MATRIX* A_bx; A_bx = new MATRIX(Norb,Norb);
  MATRIX* A_by; A_by = new MATRIX(Norb,Norb);
  MATRIX* A_bz; A_bz = new MATRIX(Norb,Norb);


  *T = (*el.E_alp) * (*el.C_alp) * (*Dao_x) * (*el.C_alp).T();
  *A_ax = (*el.C_alp) * (*dFao_alp_dx) * (*el.C_alp).T() - (*T + (*T).T());

  *T = (*el.E_alp) * (*el.C_alp) * (*Dao_y) * (*el.C_alp).T();
  *A_ay = (*el.C_alp) * (*dFao_alp_dy) * (*el.C_alp).T() - (*T + (*T).T());

  *T = (*el.E_alp) * (*el.C_alp) * (*Dao_z) * (*el.C_alp).T();
  *A_az = (*el.C_alp) * (*dFao_alp_dz) * (*el.C_alp).T() - (*T + (*T).T());

  *T = (*el.E_bet) * (*el.C_bet) * (*Dao_x) * (*el.C_bet).T();
  *A_bx = (*el.C_bet) * (*dFao_bet_dx) * (*el.C_bet).T() - (*T + (*T).T());

  *T = (*el.E_bet) * (*el.C_bet) * (*Dao_y) * (*el.C_bet).T();
  *A_by = (*el.C_bet) * (*dFao_bet_dy) * (*el.C_bet).T() - (*T + (*T).T());

  *T = (*el.E_bet) * (*el.C_bet) * (*Dao_z) * (*el.C_bet).T();
  *A_bz = (*el.C_bet) * (*dFao_bet_dz) * (*el.C_bet).T() - (*T + (*T).T());


  MATRIX* dEa_dx; dEa_dx = new MATRIX(Norb,Norb);
  MATRIX* dEa_dy; dEa_dy = new MATRIX(Norb,Norb);
  MATRIX* dEa_dz; dEa_dz = new MATRIX(Norb,Norb);
  MATRIX* dEb_dx; dEb_dx = new MATRIX(Norb,Norb);
  MATRIX* dEb_dy; dEb_dy = new MATRIX(Norb,Norb);
  MATRIX* dEb_dz; dEb_dz = new MATRIX(Norb,Norb);

  MATRIX* Dmo_a_x; Dmo_a_x = new MATRIX(Norb,Norb);
  MATRIX* Dmo_a_y; Dmo_a_y = new MATRIX(Norb,Norb);
  MATRIX* Dmo_a_z; Dmo_a_z = new MATRIX(Norb,Norb);
  MATRIX* Dmo_b_x; Dmo_b_x = new MATRIX(Norb,Norb);
  MATRIX* Dmo_b_y; Dmo_b_y = new MATRIX(Norb,Norb);
  MATRIX* Dmo_b_z; Dmo_b_z = new MATRIX(Norb,Norb);


  for(i=0;i<Norb;i++){
    dEa_dx->set(i,i,A_ax->get(i,i));
    dEa_dy->set(i,i,A_ay->get(i,i));
    dEa_dz->set(i,i,A_az->get(i,i));
    dEb_dx->set(i,i,A_bx->get(i,i));
    dEb_dy->set(i,i,A_by->get(i,i));
    dEb_dz->set(i,i,A_bz->get(i,i));

    for(j=0;j<Norb;j++){
      if(i!=j){
        double dEa = ( el.E_alp->get(j,j)-el.E_alp->get(i,i) );
        if(fabs(dEa)>1e-12){
          Dmo_a_x->set(i,j,A_ax->get(i,j)/dEa);
          Dmo_a_y->set(i,j,A_ay->get(i,j)/dEa);
          Dmo_a_z->set(i,j,A_az->get(i,j)/dEa);
        }

        double dEb = ( el.E_bet->get(j,j)-el.E_bet->get(i,i) );
        if(fabs(dEb)>1e-12){
          Dmo_b_x->set(i,j,A_bx->get(i,j)/dEb);
          Dmo_b_y->set(i,j,A_by->get(i,j)/dEb);
          Dmo_b_z->set(i,j,A_bz->get(i,j)/dEb);
        }
      }// i!=j
    }// for j
  }// for i

  // For debug:
  /*
  cout<<"dEa_dx = \n";  dEa_dx.show_matrix();
  cout<<"dEa_dy = \n";  dEa_dy.show_matrix();
  cout<<"dEa_dz = \n";  dEa_dz.show_matrix();
  cout<<"dEb_dx = \n";  dEb_dx.show_matrix();
  cout<<"dEb_dy = \n";  dEb_dy.show_matrix();
  cout<<"dEb_dz = \n";  dEb_dz.show_matrix();
  cout<<"Dmo_a_x = \n"; Dmo_a_x.show_matrix();
  cout<<"Dmo_a_y = \n"; Dmo_a_y.show_matrix();
  cout<<"Dmo_a_z = \n"; Dmo_a_z.show_matrix();
  cout<<"Dmo_b_x = \n"; Dmo_b_x.show_matrix();
  cout<<"Dmo_b_y = \n"; Dmo_b_y.show_matrix();
  cout<<"Dmo_b_z = \n"; Dmo_b_z.show_matrix();
  */

  // P = C * O * C.T(), and C^T * S * C = O  =>  O = C^T * S * P * S * C
  // The two below definitions would be good, but more expensive
  // Oa = C_alp.T() * Sao * el.get_P_alp() * Sao * C_alp
  // Ob = C_bet.T() * Sao * el.get_P_bet() * Sao * C_bet

  MATRIX* Oa;  Oa = new MATRIX(Norb,Norb);
  MATRIX* Ob;  Ob = new MATRIX(Norb,Norb);

  for(i=0;i<Norb;i++){
    Oa->set(el.occ_alp[i].first, el.occ_alp[i].first, el.occ_alp[i].second);
    Ob->set(el.occ_bet[i].first, el.occ_bet[i].first, el.occ_bet[i].second);
  }
//  *Oa = (*el.C_alp).T() * (Sao) * (*el.P_alp) * (Sao) * (*el.C_alp);
//  *Ob = (*el.C_bet).T() * (Sao) * (*el.P_bet) * (Sao) * (*el.C_bet);



  MATRIX* dP_alp_dx; dP_alp_dx = new MATRIX(Norb, Norb);
  MATRIX* dP_alp_dy; dP_alp_dy = new MATRIX(Norb, Norb);
  MATRIX* dP_alp_dz; dP_alp_dz = new MATRIX(Norb, Norb);

  MATRIX* dP_bet_dx; dP_bet_dx = new MATRIX(Norb, Norb);
  MATRIX* dP_bet_dy; dP_bet_dy = new MATRIX(Norb, Norb);
  MATRIX* dP_bet_dz; dP_bet_dz = new MATRIX(Norb, Norb);


  // Compute the derivatives of the density matrix
  *T = (*el.C_alp) * (*el.C_alp).T() * (*Dao_x) * (*el.C_alp);
  *dP_alp_dx = *(el.C_alp) * ((*Dmo_a_x) * (*Oa) - (*Oa) * (*Dmo_a_x)) * (*el.C_alp).T() - (*T + (*T).T());

  *T = (*el.C_alp) * (*el.C_alp).T() * (*Dao_y) * (*el.C_alp);
  *dP_alp_dy = *(el.C_alp) * ((*Dmo_a_y) * (*Oa) - (*Oa) * (*Dmo_a_y)) * (*el.C_alp).T() - (*T + (*T).T());

  *T = (*el.C_alp) * (*el.C_alp).T() * (*Dao_z) * (*el.C_alp);
  *dP_alp_dz = *(el.C_alp) * ((*Dmo_a_z) * (*Oa) - (*Oa) * (*Dmo_a_z)) * (*el.C_alp).T() - (*T + (*T).T());


  *T = (*el.C_bet) * (*el.C_bet).T() * (*Dao_x) * (*el.C_bet);
  *dP_bet_dx = *(el.C_bet) * ((*Dmo_a_x) * (*Oa) - (*Oa) * (*Dmo_a_x)) * (*el.C_bet).T() - (*T + (*T).T());

  *T = (*el.C_bet) * (*el.C_bet).T() * (*Dao_y) * (*el.C_bet);
  *dP_bet_dy = *(el.C_bet) * ((*Dmo_a_y) * (*Oa) - (*Oa) * (*Dmo_a_y)) * (*el.C_bet).T() - (*T + (*T).T());

  *T = (*el.C_bet) * (*el.C_bet).T() * (*Dao_z) * (*el.C_bet);
  *dP_bet_dz = *(el.C_bet) * ((*Dmo_a_z) * (*Oa) - (*Oa) * (*Dmo_a_z)) * (*el.C_bet).T() - (*T + (*T).T());

 

  VECTOR F; F = 0.0;

  // These are the forces from my derivation, but they do not work for some reason - it seems like there is some
  // sort of sign uncertainty in dP_alp_dx or some other fluctuations accumulate on top of it. Yet, theoretically
  // it looks more accurate than the version below - in tests the agreement with numerical forces is better on average,
  // but the fluctuations are terrible.
  /*
  F.x  = (  (*el.P_alp) * (*dHao_dx + *dFao_alp_dx)  ).tr()  +  ((*dP_alp_dx) * (Hao + (*el.Fao_alp)) ).tr();
  F.x += (  (*el.P_bet) * (*dHao_dx + *dFao_bet_dx)  ).tr()  +  ((*dP_bet_dx) * (Hao + (*el.Fao_bet)) ).tr();
  F.x  = -0.5 * F.x;                                                                                         
                                                                                                             
  F.y  = (  (*el.P_alp) * (*dHao_dy + *dFao_alp_dy)  ).tr()  +  ((*dP_alp_dy) * (Hao + (*el.Fao_alp)) ).tr();
  F.y += (  (*el.P_bet) * (*dHao_dy + *dFao_bet_dy)  ).tr()  +  ((*dP_bet_dy) * (Hao + (*el.Fao_bet)) ).tr();
  F.y  = -0.5 * F.y;                                                                                         
                                                                                                             
  F.z  = (  (*el.P_alp) * (*dHao_dz + *dFao_alp_dz)  ).tr()  +  ((*dP_alp_dz) * (Hao + (*el.Fao_alp)) ).tr();
  F.z += (  (*el.P_bet) * (*dHao_dz + *dFao_bet_dz)  ).tr()  +  ((*dP_bet_dz) * (Hao + (*el.Fao_bet)) ).tr();
  F.z  = -0.5 * F.z;
  */


  // But this one seems to work pretty well - at least for INDO!
  F.x  = (  (*el.P_alp) * (*dHao_dx + *dFao_alp_dx)  ).tr();
  F.x += (  (*el.P_bet) * (*dHao_dx + *dFao_bet_dx)  ).tr();
  F.x  = -0.5 * F.x;

  F.y  = (  (*el.P_alp) * (*dHao_dy + *dFao_alp_dy)  ).tr();
  F.y += (  (*el.P_bet) * (*dHao_dy + *dFao_bet_dy)  ).tr();
  F.y  = -0.5 * F.y;

  F.z  = (  (*el.P_alp) * (*dHao_dz + *dFao_alp_dz)  ).tr();
  F.z += (  (*el.P_bet) * (*dHao_dz + *dFao_bet_dz)  ).tr();
  F.z  = -0.5 * F.z;





  //================= Clean up - delete temporary memory blocks ==============
  delete dHao_dx;  delete dHao_dy;  delete dHao_dz;
  delete dSao_dx;  delete dSao_dy;  delete dSao_dz;
  delete dFao_alp_dx;  delete dFao_alp_dy;  delete dFao_alp_dz;
  delete dFao_bet_dx;  delete dFao_bet_dy;  delete dFao_bet_dz;
  delete Dao_x;  delete Dao_y;  delete Dao_z;  delete T;
  delete A_ax;  delete A_ay;  delete A_az;
  delete A_bx;  delete A_by;  delete A_bz;
  delete dEa_dx;  delete dEa_dy;  delete dEa_dz;
  delete dEb_dx;  delete dEb_dy;  delete dEb_dz;
  delete Dmo_a_x;  delete Dmo_a_y;  delete Dmo_a_z;
  delete Dmo_b_x;  delete Dmo_b_y;  delete Dmo_b_z;
  delete Oa;  delete Ob;
  delete dP_alp_dx;  delete dP_alp_dy;  delete dP_alp_dz;
  delete dP_bet_dx;  delete dP_bet_dy;  delete dP_bet_dz;



  return F;

}




double energy_and_forces
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms,Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
){
/**
  \param[in,out] el The object containing all information about electronic structure of the system
  \param[in,out] syst The object defining molecular structure of the chemical system
  \param[in,out] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the ab initio calculations
  \param[in,out] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized

  Compute the energy and forces (corresponding to the electronic structure given by el object) acting on all atoms
  of the system. Returns the total electronic + nuclear energy. The forces are returned into the corresponding variables 
  of the syst object.

  This function will perform all calculations (including the parameters) from the scratch, so that the 
  function can be used with MD or optimization context

*/

  /// Presently, configured for INDO calculations 

  int opt = 1;  // 1 - for INDO, 0 - for CNDO/CNDO2

  int i,j;

  
  //=========== STEP 1: Update positions of the basis AOs ================
  /// Update positions of the basis AOs to those given by the syst variables (atomic centers)

  for(i=0;i<el.Norb;i++){
    int at_indx = ao_to_atom_map[i];
    VECTOR at_r; at_r = syst.Atoms[at_indx].Atom_RB.rb_cm;
    at_r = at_r - basis_ao[i].primitives[0].R; // now this is the displacement

    basis_ao[i].shift_position(at_r);

  }// for i
        

  //=========== STEP 2: Depending on hamiltonian to use, set internal parameters ================
  /// Depending on hamiltonian to use, set internal parameters for faster calculations

  if(prms.hamiltonian=="eht"){ 

    vector<string> mol_at_types;
    for(i=0;i<syst.Number_of_atoms;i++){  mol_at_types.push_back(syst.Atoms[i].Atom_element); }

    set_parameters_eht_mapping(modprms, basis_ao);
    set_parameters_eht_mapping1(modprms,syst.Number_of_atoms,mol_at_types);
  }

  //=========== STEP 3: Overlap matrix ================
  /// Update the AO overlap matrix

  MATRIX* Sao; Sao = new MATRIX(el.Norb, el.Norb);

  int x_period = 0;    int y_period = 0;    int z_period = 0;
  VECTOR t1,t2,t3;

  update_overlap_matrix(x_period, y_period, z_period, t1, t2, t3, basis_ao, *Sao);


  //=========== STEP 4: Parameters ================
  /// Precompute model parameters

  vector<double> eri;  vector<double> V_AB;
     
  if(prms.hamiltonian=="indo"){
    Sao->Init_Unit_Matrix(1.0);  
    indo_core_parameters(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, opt, 0);
//    compute_all_indo_core_parameters_derivs(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, opt);

  }


  //=========== STEP 5: Core and Fock matrices ================
  /// Re-compute core and guess Fock matrices

  int debug = 0;
  double degen = 1.0;
  double kT = 0.025;
  double etol = 0.0001;
  int pop_opt = prms.pop_opt; //  #  0 -  integer populations,  1 - Fermi distribution              

  MATRIX* Hao; Hao = new MATRIX(el.Norb, el.Norb);

  Hamiltonian_core(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, *Hao,  *Sao, debug);  

  // Compute guess density
  Fock_to_P(Hao, Sao, el.Nocc_alp, degen, kT, etol, pop_opt, el.E_alp, el.C_alp, el.P_alp, el.bands_alp, el.occ_alp); 
  Fock_to_P(Hao, Sao, el.Nocc_bet, degen, kT, etol, pop_opt, el.E_bet, el.C_bet, el.P_bet, el.bands_bet, el.occ_bet);

  *el.Hao = *Hao;
  *el.Sao = *Sao;


  // Compute guess Fock matrix
  Hamiltonian_Fock(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);


  //==============  STEP 6: SCF iterations =======================
  /// Perform SCF iterations to get molecular orbitals and electronic energy

  double E = scf(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, 0);

  // Update Fock matrix at the converged density:
  Fock_to_P(el.Fao_alp, el.Sao, el.Nocc_alp, degen, kT, etol, pop_opt, el.E_alp, el.C_alp, el.P_alp, el.bands_alp, el.occ_alp); 
  Fock_to_P(el.Fao_bet, el.Sao, el.Nocc_bet, degen, kT, etol, pop_opt, el.E_bet, el.C_bet, el.P_bet, el.bands_bet, el.occ_bet);


  /*  
  cout<<"Fao_alp = \n"; el.Fao_alp->show_matrix();
  cout<<"Sao = \n"; el.Sao->show_matrix();
  cout<<"C_alp = \n"; el.C_alp->show_matrix();
  cout<<"E_alp = \n"; el.E_alp->show_matrix();
  */  



  //==============  STEP 7: Now compute forces for all atoms =====================
  // - electronic contributions
  /// Compute electronic contributions to forces for all atoms

  for(int n=0;n<syst.Number_of_atoms;n++){

    syst.Atoms[n].Atom_RB.rb_force = 
 
    force(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, *Hao, *Sao,
          el.Norb, n, x_period, y_period, z_period, t1, t2, t3);

    //cout<<"Force "<<n<<" = "<<syst.Atoms[n].Atom_RB.rb_force<<endl;

  }// for n


  // - nuclear-nuclear repulsion
  /// and then include nuclear-nuclear repulsion contributions (to energy and forces)
  vector<double> Zeff;
  vector<VECTOR> G;
  vector<VECTOR> R;
  for(i=0;i<syst.Number_of_atoms;i++){
    R.push_back(syst.Atoms[i].Atom_RB.rb_cm);
    Zeff.push_back(modprms.PT[syst.Atoms[i].Atom_element].Zeff);
  }

  double Enucl = energy_nucl(R, Zeff, G);
  for(int n=0;n<syst.Number_of_atoms;n++){   syst.Atoms[n].Atom_RB.rb_force -= G[n];   }


  E = E + Enucl;


  delete Sao;  delete Hao;
 
  /// Return the total energy (and the total forces are returned in the syst object)

  return E;

}





listHamiltonian_QM::listHamiltonian_QM(const listHamiltonian_QM& ob){   ///< Copy constructor;

  Nelec = ob.Nelec;     ///< The number of electrons in this sub-system
  Norb = ob.Norb;       ///< The total number of orbitals in this sub-system
  prms = ob.prms;       ///< Control parameters defining how to perform calculations for this sub-system
                        ///< Note: different levels of treatment are possible for different sub-systems
  modprms = ob.modprms; ///< Model parameters for this sub-system

  //============= Ground state =================
  basis_ao = ob.basis_ao;             ///< Basis for this sub-system
  atom_to_ao_map = ob.atom_to_ao_map; ///< Nuclei --> AO mapping for this sub-system
  ao_to_atom_map = ob.ao_to_atom_map; ///< AO --> Nuclei mapping for this sub-system
  el = new Electronic_Structure(ob.el);

  //============ Excited states ================
  basis_ex = ob.basis_ex;  ///< Excitations for this sub-system -  may be the same as in prms, but may be different  

}

void listHamiltonian_QM::operator=(const listHamiltonian_QM& ob){   ///< Copying one listHamiltonian_QM into the other one

  Nelec = ob.Nelec;     ///< The number of electrons in this sub-system
  Norb = ob.Norb;       ///< The total number of orbitals in this sub-system
  prms = ob.prms;       ///< Control parameters defining how to perform calculations for this sub-system
                        ///< Note: different levels of treatment are possible for different sub-systems
  modprms = ob.modprms; ///< Model parameters for this sub-system

  //============= Ground state =================
  basis_ao = ob.basis_ao;             ///< Basis for this sub-system
  atom_to_ao_map = ob.atom_to_ao_map; ///< Nuclei --> AO mapping for this sub-system
  ao_to_atom_map = ob.ao_to_atom_map; ///< AO --> Nuclei mapping for this sub-system
  *el = *ob.el;

  //============ Excited states ================
  basis_ex = ob.basis_ex;  ///< Excitations for this sub-system -  may be the same as in prms, but may be different  

}


listHamiltonian_QM::~listHamiltonian_QM(){

  if(basis_ao.size()>0){  basis_ao.clear(); }
  if(atom_to_ao_map.size()){  atom_to_ao_map.clear(); }
  if(ao_to_atom_map.size()){  ao_to_atom_map.clear(); }
  if(basis_ex.size()){ basis_ex.size(); }

  delete el;

}


listHamiltonian_QM::listHamiltonian_QM(std::string ctrl_filename,System& syst){
/**
  \param[in] ctrl_filename The name of the file that defines all the settings for the electronic structure calculations, including
  the information on which atomistic Hamiltonian parameters to use
  \param[in] syst The object containing structural information about system. 

  Performes a sequanece of the initialization steps to form the QM Hamiltonian (eventually), but doesn't form the Hamiltonian yet
  Also, adds a "ground state excitation" - the reference configuration, so that the excitonic basis of lenght 1 is defined
  Optionally: Add more excitonic states (by using add_excitation) for excited-states calculations
  The function listHamiltonian_QM::compute_scf() can be safely called after this constructor
*/


  init(ctrl_filename, syst);
  add_excitation(0,1,0,1);

}

void listHamiltonian_QM::init(std::string ctrl_filename,System& syst){
/**
  \param[in] ctrl_filename The name of the file that defines all the settings for the electronic structure calculations, including
  the information on which atomistic Hamiltonian parameters to use
  \param[in] syst The object containing structural information about system. 

  Performes a sequanece of the initialization steps to form the QM Hamiltonian (eventually), but doesn't form the Hamiltonian yet
  Usually should follow by a set of add_excitation() function calls, to create electronic (excitonic) basis 
*/


  //=========== STEP 1: Create control parameters (setting computation options) ================
  /// Create the Control_Paramters object from the ctrl_filename file

  libcontrol_parameters::get_parameters_from_file(ctrl_filename, prms);


  //=========== STEP 2:  Create model parameters and load them from file (using control parameters options) ================
  /// Initialize/read model parameters (need basis info)

  if(prms.hamiltonian=="eht"/* or prms.hamiltonian=="geht"*/){
    set_parameters_eht(prms, modprms);
  }
  else if(prms.hamiltonian=="indo"){
    set_parameters_indo(prms, modprms);
  }
//  else if(prms.hamiltonian=="geht1"){
//    set_parameters_geht1(prms, modprms); 
//  }
//  else if(prms.hamiltonian=="geht2"){
//    set_parameters_geht2(prms, modprms); 
//  }

  //=========== STEP 3: Set basis (STO-3G_DZ) ================
  /// Set STO-3G_DZ basis

  //------- Input --------------
  vector<std::string> mol_at_types;
  vector<VECTOR> R;

  for(int i=0; i<syst.Number_of_atoms;i++){
    mol_at_types.push_back(syst.Atoms[i].Atom_element);
    R.push_back(syst.Atoms[i].Atom_RB.rb_cm);
  }

  //-------- Output -----------
  int verb = 0;
  set_basis_STO_3G_DZ(mol_at_types, R, modprms, verb, basis_ao, Nelec, Norb, atom_to_ao_map, ao_to_atom_map);


   //=========== STEP 4: Electronic structure ================
   /// Ctreate the Electronic_Structure object and initialize its numbers

   el = new Electronic_Structure(Norb);
   el->Nelec = Nelec;
   el->Nocc_alp = Nelec/2;
   el->Nocc_bet = Nelec - el->Nocc_alp;


  //========== STEP 5: Depending on hamiltonian to use, set internal parameters ================
  /// Depending on hamiltonian to use, set internal parameters for faster calculations
  /// this step runs after AO basis is set!!!

  if(prms.hamiltonian=="eht" /*|| prms.hamiltonian=="geht" ||
     prms.hamiltonian=="geht1" || prms.hamiltonian=="geht2"*/){
      set_parameters_eht_mapping(modprms, basis_ao);
      set_parameters_eht_mapping1(modprms,syst.Number_of_atoms,mol_at_types);
  }

  //=========== STEP 6: Overlap matrix ================
  /// Precompute AO overlap matrix and store it in the created above Electronic_Structure object

  int x_period = 0;
  int y_period = 0;
  int z_period = 0;
  VECTOR t1, t2, t3;

  update_overlap_matrix(x_period, y_period, z_period, t1, t2, t3, basis_ao, *el->Sao); 

  //=========== STEP 7: Method-specific Parameters ================
  /// Set up Hamiltonian-type-specific parameters

  if(prms.hamiltonian=="indo"){
    /// For INDO: a) set overlap matrix to the identiy matrix; b) precompute INDO core parameters

    int opt = 1; 
    el->Sao->Init_Unit_Matrix(1.0);  
    indo_core_parameters(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, opt,1);
  }
  else if(prms.hamiltonian=="cndo"||prms.hamiltonian=="cndo2"){
    /// For CNDO and CNDO/2: a) set overlap matrix to the identiy matrix; b) precompute CNDO/2 core parameters

    int opt = 0; 
    el->Sao->Init_Unit_Matrix(1.0);  
    indo_core_parameters(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, opt,1);
  }

/*
  //=========== STEP 8: Core Hamiltonian ================
  int debug = 1;
  Hamiltonian_core(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, *el->Hao,  *el->Sao, debug);


  //=========== STEP 9: Guess density matrix ================
  std::string eigen_method="generalized";
  vector<Timer> bench_t2(4);

  Fock_to_P(el->Norb, el->Nocc_alp, 1, el->Nocc_alp, eigen_method, prms.pop_opt,
            el->Fao_alp, el->Sao, el->C_alp, el->E_alp, el->bands_alp, el->occ_alp, el->P_alp,
            bench_t2);

  Fock_to_P(el->Norb, el->Nocc_bet, 1, el->Nocc_bet, eigen_method, prms.pop_opt,
            el->Fao_bet, el->Sao, el->C_bet, el->E_bet, el->bands_bet, el->occ_bet, el->P_bet,
            bench_t2);

  //=========== STEP 10: SCF solution ================
  scf(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, 0); 
*/
    

}


void listHamiltonian_QM::compute_overlap(System& syst){
/**
  \param[in,out] syst The object containing structural information about system. 

  Computes AO basis overlap matrix

  Use after listHamiltonian_QM::init() and add_excitation() or non-default constructor !!!
*/

  int x_period = 0;
  int y_period = 0;
  int z_period = 0;
  VECTOR t1, t2, t3;

  update_overlap_matrix(x_period, y_period, z_period, t1, t2, t3, basis_ao, *el->Sao); 

}


void listHamiltonian_QM::compute_core_Hamiltonian(System& syst){
/**
  \param[in,out] syst The object containing structural information about system. 

  Computes core Hamiltonian of the system (initialization must be performed first)

  Use after listHamiltonian_QM::init() and add_excitation() or non-default constructor !!!
*/


  int debug = 0;
  Hamiltonian_core(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, *el->Hao,  *el->Sao, debug);


}



double listHamiltonian_QM::compute_scf(System& syst){
/**
  \param[in,out] syst The object containing structural information about system. It also will contain the forces on active state
  in the end of this calculation.

  Performes a sequanece of the preparatory steps for the SCF as well as the SCF calculations for the quantum system

  Use after listHamiltonian_QM::init() and add_excitation()  !!!
*/


/*
  //=========== STEP 1: Create control parameters (setting computation options) ================
  libcontrol_parameters::get_parameters_from_file(ctrl_filename, prms);


  //=========== STEP 2:  Create model parameters and load them from file (using control parameters options) ================
  // Initialize/read model parameters (need basis info)
  if(prms.hamiltonian=="eht" or prms.hamiltonian=="geht"){
    set_parameters_eht(prms, modprms);
  }
  else if(prms.hamiltonian=="indo"){
    set_parameters_indo(prms, modprms);
  }
  else if(prms.hamiltonian=="geht1"){
    set_parameters_geht1(prms, modprms); 
  }
  else if(prms.hamiltonian=="geht2"){
    set_parameters_geht2(prms, modprms); 
  }

  //=========== STEP 3: Set basis (STO-3G_DZ) ================
  //------- Input --------------
  vector<std::string> mol_at_types;
  vector<VECTOR> R;

  for(int i=0; i<syst.Number_of_atoms;i++){
    mol_at_types.push_back(syst.Atoms[i].Atom_element);
    R.push_back(syst.Atoms[i].Atom_RB.rb_cm);
  }

  //-------- Output -----------
  int verb = 0;
  set_basis_STO_3G_DZ(mol_at_types, R, modprms, verb, basis_ao, Nelec, Norb, atom_to_ao_map, ao_to_atom_map);


   //=========== STEP 4: Electronic structure ================
   el = new Electronic_Structure(Norb);
   el->Nelec = Nelec;
   el->Nocc_alp = Nelec/2;
   el->Nocc_bet = Nelec - el->Nocc_alp;


  //=========== STEP 5: Depending on hamiltonian to use, set internal parameters ================
  // this step runs after AO basis is set!!!

  if(prms.hamiltonian=="eht" || prms.hamiltonian=="geht" ||
     prms.hamiltonian=="geht1" || prms.hamiltonian=="geht2"){
      set_parameters_eht_mapping(modprms, basis_ao);
      set_parameters_eht_mapping1(modprms,syst.Number_of_atoms,mol_at_types);
  }
  //=========== STEP 6: Overlap matrix ================
  int x_period = 0;
  int y_period = 0;
  int z_period = 0;
  VECTOR t1, t2, t3;

  update_overlap_matrix(x_period, y_period, z_period, t1, t2, t3, basis_ao, *el->Sao); 

  //=========== STEP 7: Method-specific Parameters ================
  if(prms.hamiltonian=="indo"){
    int opt = 1; 
    el->Sao->Init_Unit_Matrix(1.0);  
    indo_core_parameters(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, opt,1);
  }
  else if(prms.hamiltonian=="cndo"||prms.hamiltonian=="cndo2"){
    int opt = 0; 
    el->Sao->Init_Unit_Matrix(1.0);  
    indo_core_parameters(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, opt,1);
  }

*/
  //=========== STEP 8: Core Hamiltonian ================
  /// Form core Hamiltonian

  int debug = 0;
  Hamiltonian_core(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, *el->Hao,  *el->Sao, debug);


  //=========== STEP 9: Guess density matrix ================
  /// Compute guss density matrix (from the core Hamiltonian)

  std::string eigen_method="generalized";
  vector<Timer> bench_t2(4);

  Fock_to_P(el->Norb, el->Nocc_alp, 1, el->Nocc_alp, eigen_method, prms.pop_opt,
            el->Fao_alp, el->Sao, el->C_alp, el->E_alp, el->bands_alp, el->occ_alp, el->P_alp,
            bench_t2);

  Fock_to_P(el->Norb, el->Nocc_bet, 1, el->Nocc_bet, eigen_method, prms.pop_opt,
            el->Fao_bet, el->Sao, el->C_bet, el->E_bet, el->bands_bet, el->occ_bet, el->P_bet,
            bench_t2);

  //=========== STEP 10: SCF solution ================
  /// Perform the SCF calculations

  double E = scf(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, 0); 

  /// Return the total electronic energy of the system

  return E;    

}


void listHamiltonian_QM::set_electronic_structure(Electronic_Structure& el_){
/**
  \param[in,out] el_ The object containing all the information about electronic structure of the quantum subsystem.
  Make the internal pointer to Electronic_Structure object point to an external object
  This way also the changes made externally will affect the internal state, and vice versa - the 
  internal modifications and updates will be available in the external object.
*/

  el = &el_;  // by reference

}

double listHamiltonian_QM::energy_and_forces(System& syst){ 
/**
  \param[in,out] syst The object containing structural information about system. It also will contain the forces on active state
  in the end of this calculation.

  Computes energy and forces (inclding excited state) for quantum part of the chemical system (given by syst)
*/

  return libhamiltonian_qm::energy_and_forces(*el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map );

}


void listHamiltonian_QM::add_excitation(int f_o, int f_s, int t_o, int t_s){
/**
  \param[in] f_o "from orbital" The index of the orbital from which electron is excited
  \param[in] f_s "from spin" The index of the spin from which electron is excited (1 - alpha; -1 - beta)
  \param[in] t_o "to orbital" The index of the orbital to which electron is excited
  \param[in] t_s "to spin" The index of the spin to which electron is excited (1 - alpha; -1 - beta)

  Creates a new excitation and adds it to quantum Hamiltonian - as one of the basis states for 
  excited states calculations

  Note: It is important to create such excitations before attempting NA-MD calculations 
  The number of such excitation (calls of this + excite_bet function) must be the number of electronic 
  states in the Electronic object (not the Electronic_Structure !!!) - 1 (ground state is already included)
*/

  int sz = basis_ex.size();

  if(sz==0){ basis_ex.push_back(excitation(f_o, f_s, t_o, t_s)); }
  else{

    int is_new = 1;
    for(int i=0;i<sz;i++){ 
      if(basis_ex[i].from_orbit[0]==f_o){
        if(basis_ex[i].from_spin[0]==f_s){
          if(basis_ex[i].to_orbit[0]==t_o){
            if(basis_ex[i].to_spin[0]==t_s){
              is_new = 0;
            }
          }
        }
      }
    }// for i

    if(is_new==1){  basis_ex.push_back(excitation(f_o, f_s, t_o, t_s)); }

  }

}

//void listHamiltonian_QM::set_excitonic_basis(boost::python::list basis_ex){
//}

void listHamiltonian_QM::excite_alp(int I, int J){
/**
  \param[in] I the index of the source orbital involved in the alpha-excitation
  \param[in] J the index of the target orbital involved in the alpha-excitation

  This function creates an excitation - the electronic excited state - the basis state for NA-MD
  This excitation is for the alpha electron going from orbital with index I to that with index J (I-->J)  
*/


  el->excite_alp(I,J);

}

void listHamiltonian_QM::excite_bet(int I, int J){
/**
  \param[in] I the index of the source orbital involved in the beta-excitation
  \param[in] J the index of the target orbital involved in the beta-excitation

  This function creates an excitation - the electronic excited state - the basis state for NA-MD
  This excitation is for the beta electron going from orbital with index I to that with index J (I-->J)  
*/


  el->excite_bet(I,J);

}




}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra


