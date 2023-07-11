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
  \file Hamiltonian_EHT.cpp
  \brief The file implements functions for extended Huckel theory (EHT) calculations
*/

#include "Hamiltonian_EHT.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{



void Hamiltonian_core_eht
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[out] Hao The pointer to the matrix object in which the core Hamiltonian will be stored
  \param[in] Sao The pointer to the AO overlap matrix 
  \param[in] DF Debug flag - controlls how much of extra info is printed out
  
  Compute the core EHT (extended Huckel theory) Hamiltonian

  Options:
  prms.eht_formula == 0 - unweighted 
  prms.eht_formula == 1 - weighted
  prms.eht_formula == 2 - Calzaferi
  prms.eht_formula == 3 - for the developments

*/

           
  int i,j,n, a, b, I,J;
  double delt, delt2, delt4;

  int sz = syst.Number_of_atoms; // number of atoms in this fragment  

  int Norb = basis_ao.size(); // how many AOs are included in this fragment
  if(Norb!=Hao->n_cols){  
    cout<<"In Hamiltonian_core_eht: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
 
  //========================= Core diagonal elements ==========================  
  // Diagonal elements = set to IPs
  for(i=0;i<Norb;i++){                                              //    old
    Hao->M[i*Norb+i] = Hao->M[i*Norb+i] = modprms.orb_params[i].IP; //modprms.PT[basis_ao[i].element].IP[basis_ao[i].ao_shell];
  }// for i


  // Off-diagonal elements
  double K_const = 0.0;// default value

  // Off-diagonal elements
  for(i=0;i<Norb;i++){
    a = ao_to_atom_map[i];
    for(j=i+1;j<Norb;j++){          
        b = ao_to_atom_map[j];
        // This is old and slow version
        //K_const = modprms.eht_k.get_K_value(basis_ao[i].element,basis_ao[i].ao_shell,basis_ao[j].element,basis_ao[j].ao_shell);
        K_const = modprms.meht_k.get_K_value(0, i,j);
        //cout<<"i= "<<i<<" j= "<<j<<" K_const= "<<K_const<<endl;

        if(prms.eht_formula==0){  // Unweighted formula

          Hao->M[i*Norb+j] = 0.5*K_const*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j])*Sao->M[i*Norb+j]; 
          Hao->M[j*Norb+i] = Hao->M[i*Norb+j];
        }

        else if(prms.eht_formula==1){  // Weighted formula:        

          delt = (Hao->M[i*Norb+i]-Hao->M[j*Norb+j])/(Hao->M[i*Norb+i]+Hao->M[j*Norb+j]);
          delt2 = delt*delt;
          delt4 = delt2*delt2;
        
          Hao->M[i*Norb+j] = 0.5*(K_const + delt2 + (1.0 - K_const)*delt4)*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j])*Sao->M[i*Norb+j];
          Hao->M[j*Norb+i] = Hao->M[i*Norb+j];

        }

        else if(prms.eht_formula==2){  // Calzaferi formula:        

          delt = (Hao->M[i*Norb+i]-Hao->M[j*Norb+j])/(Hao->M[i*Norb+i]+Hao->M[j*Norb+j]);
          delt2 = delt*delt;
          delt4 = delt2*delt2;
          
          double rab = (syst.Atoms[a].Atom_RB.rb_cm - syst.Atoms[b].Atom_RB.rb_cm).length();
          double delta = 0.13;


          std::string elt_i = basis_ao[i].element;
          std::string elt_j = basis_ao[j].element;
          std::string sh_i  = basis_ao[i].ao_shell;
          std::string sh_j  = basis_ao[j].ao_shell;


          float n_i = modprms.PT[elt_i].Nquant[sh_i];
          float n_j = modprms.PT[elt_j].Nquant[sh_j];
          int nz_i  = modprms.PT[elt_i].Nzeta[sh_i];
          int nz_j  = modprms.PT[elt_j].Nzeta[sh_j];

          double ri, rj;  // radii

          // for atom i 
          if(nz_i==1){  ri = (n_i/modprms.PT[elt_i].zetas[sh_i][0]); }
          else if(nz_i==2){
            double z1 = modprms.PT[elt_i].zetas[sh_i][0];
            double z2 = modprms.PT[elt_i].zetas[sh_i][1];
            double c1 = modprms.PT[elt_i].coeffs[sh_i][0];
            double c2 = modprms.PT[elt_i].coeffs[sh_i][1];

            ri = n_i/(c1*c1*z1 + c2*c2*z2 + 
                       ( pow(2.0,2.0*n_i)*pow(z1*z2, n_i+0.5)/pow((z1+z2),2.0*n_i)
                       ) 
                     ); 

          }

          // for atom j
          if(nz_j==1){  rj = (n_j/modprms.PT[elt_j].zetas[sh_j][0]); }
          else if(nz_j==2){
            double z1 = modprms.PT[elt_j].zetas[sh_j][0];
            double z2 = modprms.PT[elt_j].zetas[sh_j][1];
            double c1 = modprms.PT[elt_j].coeffs[sh_j][0];
            double c2 = modprms.PT[elt_j].coeffs[sh_j][1];

            rj = n_j/(c1*c1*z1 + c2*c2*z2 + 
                       ( pow(2.0,2.0*n_j)*pow(z1*z2, n_j+0.5)/pow((z1+z2),2.0*n_j)
                       ) 
                     ); 

          }

                                             
          double d0 = ri + rj;    

          K_const = 1.0 + (0.75 + delt2 - 0.75*delt4)*exp(-delta*(rab - d0));
        
          Hao->M[i*Norb+j] = 0.5*K_const*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j])*Sao->M[i*Norb+j];
          Hao->M[j*Norb+i] = Hao->M[i*Norb+j];

        }

        else if(prms.eht_formula==3){       
         // Add your options here and next

        }// ==3
     
    }// for bj
  }// for i


}


void Hamiltonian_core_eht
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[out] Hao The matrix object in which the core Hamiltonian will be stored
  \param[in] Sao The AO overlap matrix 
  \param[in] DF Debug flag - controlls how much of extra info is printed out
  
  Compute the core EHT (extended Huckel theory) Hamiltonian - Python-friendly version

  Options:
  prms.eht_formula == 0 - unweighted 
  prms.eht_formula == 1 - weighted
  prms.eht_formula == 2 - Calzaferi
  prms.eht_formula == 3 - for the developments

*/


  Hamiltonian_core_eht( syst, basis_ao, prms, modprms,  atom_to_ao_map, ao_to_atom_map, &Hao, &Sao, DF);

}


void Hamiltonian_core_deriv_eht
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF,
  int c,
  MATRIX* dHao_dx, MATRIX* dHao_dy, MATRIX* dHao_dz, 
  MATRIX* dSao_dx, MATRIX* dSao_dy, MATRIX* dSao_dz
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[out] Hao The pointer to the matrix object in which the core Hamiltonian will be stored
  \param[out] Sao The pointer to the AO overlap matrix computed here
  \param[in] DF Debug flag - controlls how much of extra info is printed out
  \param[in] c The index of the atom w.r.t. which coordinates we take the derivatives
  \param[out] dHao_dx The derivative of the Hamiltonian w.r.t. the x-coordinate of the selected atom
  \param[out] dHao_dy The derivative of the Hamiltonian w.r.t. the y-coordinate of the selected atom
  \param[out] dHao_dz The derivative of the Hamiltonian w.r.t. the z-coordinate of the selected atom
  \param[out] dSao_dx The derivative of the AO overlap matrix w.r.t. the x-coordinate of the selected atom
  \param[out] dSao_dy The derivative of the AO overlap matrix w.r.t. the y-coordinate of the selected atom
  \param[out] dSao_dz The derivative of the AO overlap matrix w.r.t. the z-coordinate of the selected atom

  Compute the core EHT (extended Huckel theory) Hamiltonian and its derivatives w.r.t. specific nuclear DOFs

  Options:
  prms.eht_formula == 0 - unweighted 
  prms.eht_formula == 1 - weighted
  prms.eht_formula == 2 - Calzaferi
  prms.eht_formula == 3 - for the developments
  
*/


  //================ Basically, here we compute derivatives of the core Hamiltonian ========================

  int i,j,k,a,b,I,J,A,B;
  double delt, delt2, delt4;
  VECTOR dIdA,dIdB;

  int Norb = basis_ao.size(); // how many AOs are included in this fragment
  if(Norb!=Hao->n_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_eht is called\n";
    cout<<"In Hamiltonian_core_deriv_eht: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  if(Norb!=dHao_dx->n_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_eht is called\n";
    cout<<"In Hamiltonian_core_deriv_eht: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  if(Norb!=dHao_dy->n_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_eht is called\n";
    cout<<"In Hamiltonian_core_deriv_eht: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  if(Norb!=dHao_dz->n_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_eht is called\n";
    cout<<"In Hamiltonian_core_deriv_eht: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }

  int sz = syst.Number_of_atoms; // number of atoms in this fragment



  stringstream ss(stringstream::in | stringstream::out);
  std::string out;
  (ss << c);  ss >> out;



  //----------- Compute core terms of the Hamiltonian ------------
  *Hao = 0.0;
  *dHao_dx = 0.0;
  *dHao_dy = 0.0;
  *dHao_dz = 0.0;

  *dSao_dx = 0.0;
  *dSao_dy = 0.0;
  *dSao_dz = 0.0;


  for(i=0;i<Norb;i++){  // global orbital indices

    // Diagonal elements = set to IPs
    // No contributions to diagonal elements of gradients
    a = ao_to_atom_map[i];
                                                                    //            old
    Hao->M[i*Norb+i] = Hao->M[i*Norb+i] = modprms.orb_params[i].IP; // modprms.PT[basis_ao[i].element].IP[basis_ao[i].ao_shell];
    
    
    //-------------- Off-diagonal terms of the core matrix ---------
    for(j=0;j<Norb;j++){

      if(j!=i){
        b = ao_to_atom_map[j];

        if(b==a){ ;; }  // different orbitals centered on the same atom - give zero (not true for hybrid orbitals)
        else{           // centered on different atoms - use overlap formula


          VECTOR dSda, dSdb, dSdc;
          double sao_ij = gaussian_overlap(&basis_ao[i],&basis_ao[j], 1, 1, dSda, dSdb);
          // i - on atom a
          // j - on atom b

          dSdc = 0.0;
          if(c==a){ dSdc += dSda; }
          if(c==b){ dSdc += dSdb; }
          

          double K_const = modprms.meht_k.get_K_value(0, i,j);

          if(prms.eht_formula==0){  // Unweighted formula

            double beta_ij = 0.5*K_const*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j]);

            Hao->M[i*Norb+j] = beta_ij*Sao->M[i*Norb+j]; 

            dHao_dx->M[i*Norb+j] += beta_ij * dSdc.x;
            dHao_dy->M[i*Norb+j] += beta_ij * dSdc.y;
            dHao_dz->M[i*Norb+j] += beta_ij * dSdc.z;

            dSao_dx->M[i*Norb+j] += dSdc.x;
            dSao_dy->M[i*Norb+j] += dSdc.y;
            dSao_dz->M[i*Norb+j] += dSdc.z;


          }// eht_formula == 0

          else if(prms.eht_formula==1){  // Weighted formula:        

            delt = (Hao->M[i*Norb+i]-Hao->M[j*Norb+j])/(Hao->M[i*Norb+i]+Hao->M[j*Norb+j]);
            delt2 = delt*delt;
            delt4 = delt2*delt2;
        
            double beta_ij = 0.5*(K_const + delt2 + (1.0 - K_const)*delt4)*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j]);

            Hao->M[i*Norb+j] = beta_ij * Sao->M[i*Norb+j];

            dHao_dx->M[i*Norb+j] += beta_ij * dSdc.x;
            dHao_dy->M[i*Norb+j] += beta_ij * dSdc.y;
            dHao_dz->M[i*Norb+j] += beta_ij * dSdc.z;

            dSao_dx->M[i*Norb+j] += dSdc.x;
            dSao_dy->M[i*Norb+j] += dSdc.y;
            dSao_dz->M[i*Norb+j] += dSdc.z;


          }
          else if(prms.eht_formula==2){  // Calzaferi formula:        

            delt = (Hao->M[i*Norb+i]-Hao->M[j*Norb+j])/(Hao->M[i*Norb+i]+Hao->M[j*Norb+j]);
            delt2 = delt*delt;
            delt4 = delt2*delt2;
          
            double rab = (syst.Atoms[a].Atom_RB.rb_cm - syst.Atoms[b].Atom_RB.rb_cm).length();
            double delta = 0.13;


            std::string elt_i = basis_ao[i].element;
            std::string elt_j = basis_ao[j].element;
            std::string sh_i  = basis_ao[i].ao_shell;
            std::string sh_j  = basis_ao[j].ao_shell;


            float n_i = modprms.PT[elt_i].Nquant[sh_i];
            float n_j = modprms.PT[elt_j].Nquant[sh_j];
            int nz_i  = modprms.PT[elt_i].Nzeta[sh_i];
            int nz_j  = modprms.PT[elt_j].Nzeta[sh_j];

            double ri, rj;  // radii

            // for atom i 
            if(nz_i==1){  ri = (n_i/modprms.PT[elt_i].zetas[sh_i][0]); }
            else if(nz_i==2){
              double z1 = modprms.PT[elt_i].zetas[sh_i][0];
              double z2 = modprms.PT[elt_i].zetas[sh_i][1];
              double c1 = modprms.PT[elt_i].coeffs[sh_i][0];
              double c2 = modprms.PT[elt_i].coeffs[sh_i][1];

              ri = n_i/(c1*c1*z1 + c2*c2*z2 + 
                         ( pow(2.0,2.0*n_i)*pow(z1*z2, n_i+0.5)/pow((z1+z2),2.0*n_i)
                         ) 
                       ); 

            }

            // for atom j
            if(nz_j==1){  rj = (n_j/modprms.PT[elt_j].zetas[sh_j][0]); }
            else if(nz_j==2){
              double z1 = modprms.PT[elt_j].zetas[sh_j][0];
              double z2 = modprms.PT[elt_j].zetas[sh_j][1];
              double c1 = modprms.PT[elt_j].coeffs[sh_j][0];
              double c2 = modprms.PT[elt_j].coeffs[sh_j][1];

              rj = n_j/(c1*c1*z1 + c2*c2*z2 + 
                         ( pow(2.0,2.0*n_j)*pow(z1*z2, n_j+0.5)/pow((z1+z2),2.0*n_j)
                         ) 
                       ); 
            }

            double d0 = ri + rj;    

            K_const = 1.0 + (0.75 + delt2 - 0.75*delt4)*exp(-delta*(rab - d0));
        
            double beta_ij = 0.5*K_const*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j]);
             

            VECTOR dbeta_ij_dr; dbeta_ij_dr = 0.0;
            if(c==a){
              dbeta_ij_dr = 0.5*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j])*
                       (-delta*(syst.Atoms[a].Atom_RB.rb_cm - syst.Atoms[b].Atom_RB.rb_cm).unit())*(K_const - 1.0);
            }
            if(c==b){
              dbeta_ij_dr -= 0.5*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j])*
                       (-delta*(syst.Atoms[a].Atom_RB.rb_cm - syst.Atoms[b].Atom_RB.rb_cm).unit())*(K_const - 1.0);
            }


            Hao->M[i*Norb+j] = beta_ij * Sao->M[i*Norb+j];

            dHao_dx->M[i*Norb+j] += (beta_ij * dSdc.x + dbeta_ij_dr.x * Sao->M[i*Norb+j]);
            dHao_dy->M[i*Norb+j] += (beta_ij * dSdc.y + dbeta_ij_dr.y * Sao->M[i*Norb+j]);
            dHao_dz->M[i*Norb+j] += (beta_ij * dSdc.z + dbeta_ij_dr.z * Sao->M[i*Norb+j]);

            dSao_dx->M[i*Norb+j] += dSdc.x;
            dSao_dy->M[i*Norb+j] += dSdc.y;
            dSao_dz->M[i*Norb+j] += dSdc.z;


        }

        else if(prms.eht_formula==3){       
         // Add your options here and next

        }// ==3




        }// else a!=b
      }// j!=i
    }// for j    
  }// for i


}


void Hamiltonian_core_deriv_eht
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF,
  int c,
  MATRIX& dHao_dx, MATRIX& dHao_dy, MATRIX& dHao_dz, 
  MATRIX& dSao_dx, MATRIX& dSao_dy, MATRIX& dSao_dz
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[out] Hao The matrix object in which the core Hamiltonian will be stored
  \param[out] Sao The AO overlap matrix computed here
  \param[in] DF Debug flag - controlls how much of extra info is printed out
  \param[in] c The index of the atom w.r.t. which coordinates we take the derivatives
  \param[out] dHao_dx The derivative of the Hamiltonian w.r.t. the x-coordinate of the selected atom
  \param[out] dHao_dy The derivative of the Hamiltonian w.r.t. the y-coordinate of the selected atom
  \param[out] dHao_dz The derivative of the Hamiltonian w.r.t. the z-coordinate of the selected atom
  \param[out] dSao_dx The derivative of the AO overlap matrix w.r.t. the x-coordinate of the selected atom
  \param[out] dSao_dy The derivative of the AO overlap matrix w.r.t. the y-coordinate of the selected atom
  \param[out] dSao_dz The derivative of the AO overlap matrix w.r.t. the z-coordinate of the selected atom

  Compute the core EHT (extended Huckel theory) Hamiltonian and its derivatives w.r.t. specific nuclear DOFs - Python-friendly version

  Options:
  prms.eht_formula == 0 - unweighted 
  prms.eht_formula == 1 - weighted
  prms.eht_formula == 2 - Calzaferi
  prms.eht_formula == 3 - for the developments
  
*/


  Hamiltonian_core_deriv_eht
  ( syst, basis_ao, prms, modprms,  atom_to_ao_map, ao_to_atom_map,  &Hao, &Sao, DF, c,
    &dHao_dx, &dHao_dy, &dHao_dz,   &dSao_dx, &dSao_dy, &dSao_dz);

}



void Hamiltonian_Fock_eht(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
                          Control_Parameters& prms, Model_Parameters& modprms,
                          vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                         ){
/**
  \param[in,out] el The electronic structre properties of the system
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  
  Compute the EHT Fock Hamiltonian. Well, we are now talking about even further generalized EHT - the one that includes
  self-consistent charges and so. That is why we get density-matrix-dependent Fock matrix.
  Our options:

  prms.eht_sce_formula==0  - charge-independent matrix, so F = H (core)
  prms.eht_sce_formula==1  - we first add orbital energy correction that is linearly-dependent on atomic Mulliken charges, then
                             the corrected orbital energies are used in one of the mixing formulas:
                             
    prms.eht_formula == 0 - unweighted 
    prms.eht_formula == 1 - weighted
    prms.eht_formula == 2 - Calzaferi
    prms.eht_formula == 3 - for the developments

  On top of this all, we also have a QE-type of correction:

  prms.eht_electrostatics>=1

    a) this gives on-site Fock matrix correction of the INDO type (breaking spin symmetry, in general) - this is actually controlled by the
    model parameters

    b) also, this correction gives the Fock matrix terms that originate from the energy term: xi_a * Q_a + 1/2*J_aa * Q_a^2 -on-site energies

  prms.eht_electrostatics>=2

    same as   prms.eht_electrostatics>=1 but also the terms originating from the energy terms 1/2 * J_ab * Q_a * Q_b are added to the 
    Fock matrix

*/


  int i,j,k,n,I,J,K,a,b,A,B;
  double delt, delt2, delt4;

  int Norb = basis_ao.size(); // how many AOs are included in this fragment
  if(Norb!=el->Hao->n_cols){  
    cout<<"In Hamiltonian_Fock_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }


  // Formation of the Fock matrix: Core part
  *el->Fao_alp = 0.0; 
  *el->Fao_bet = 0.0; 


  //============ Update charges and populations ==========================
  *el->P = *el->P_alp + *el->P_bet;


  update_Mull_orb_pop(el->P, el->Sao, el->Mull_orb_pop_gross, el->Mull_orb_pop_net);

  vector<double> Zeff(syst.Number_of_atoms, 0.0);
  vector<double> Mull_charges_gross(syst.Number_of_atoms, 0.0);
  vector<double> Mull_charges_net(syst.Number_of_atoms, 0.0);

  for(a=0;a<syst.Number_of_atoms;a++){ Zeff[a] = modprms.PT[syst.Atoms[a].Atom_element].Zeff; } // e.g. 4 for STO-3G C

  update_Mull_charges(ao_to_atom_map, Zeff, el->Mull_orb_pop_gross, el->Mull_orb_pop_net, Mull_charges_gross, Mull_charges_net);

  for(a=0;a<syst.Number_of_atoms;a++){ 
    syst.Atoms[a].Atom_mull_charge_gross = Mull_charges_gross[a]; 
    syst.Atoms[a].Atom_mull_charge_net = Mull_charges_net[a]; 
  }



  //============== Charge-corrected Fock matrix ==============================
  // Depending of SCE, do different type of correction to diagonal elements
  if(prms.eht_sce_formula==0){ 
    *el->Fao_alp = *el->Hao;
    *el->Fao_bet = *el->Hao;

  }
  else if(prms.eht_sce_formula==1){

    // Now modify IPs - diagonal elements
    for(int i=0;i<Norb;i++){
      int a = ao_to_atom_map[i];
      double Ai = 0.0;  // typically a positive number
      double Q = syst.Atoms[a].Atom_mull_charge_gross; // charge of the atom, on which i-th AO is sitting

      if(Q>0){ Ai = modprms.PT[basis_ao[a].element].J_param1[basis_ao[a].ao_shell]; }
      else{    Ai = modprms.PT[basis_ao[a].element].J_param2[basis_ao[a].ao_shell]; }

      el->Fao_alp->M[i*Norb+i] = el->Hao->M[i*Norb+i] - (Ai * Q);

    }// for a


    // Off-diagonal elements
    double K_const = 0.0;// default value
  
    for(i=0;i<Norb;i++){
      a = ao_to_atom_map[i];
      for(j=i+1;j<Norb;j++){          
          b = ao_to_atom_map[j];
          // This is old and slow version
          //K_const = modprms.eht_k.get_K_value(basis_ao[i].element,basis_ao[i].ao_shell,basis_ao[j].element,basis_ao[j].ao_shell);
          K_const = modprms.meht_k.get_K_value(0,i,j);
  
          if(prms.eht_formula==0){  // Unweighted formula
  
            el->Fao_alp->M[i*Norb+j] = 0.5*K_const*(el->Fao_alp->M[i*Norb+i] + el->Fao_alp->M[j*Norb+j]) * el->Sao->M[i*Norb+j]; 
            el->Fao_alp->M[j*Norb+i] = el->Fao_alp->M[i*Norb+j];
          }
  
          else if(prms.eht_formula==1){  // Weighted formula:        
  
            delt = (el->Fao_alp->M[i*Norb+i] - el->Fao_alp->M[j*Norb+j])/(el->Fao_alp->M[i*Norb+i] + el->Fao_alp->M[j*Norb+j]);
            delt2 = delt*delt;
            delt4 = delt2*delt2;
          
            el->Fao_alp->M[i*Norb+j] = 0.5*(K_const + delt2 + (1.0 - K_const)*delt4)*(el->Fao_alp->M[i*Norb+i]+el->Fao_alp->M[j*Norb+j])*el->Sao->M[i*Norb+j];
            el->Fao_alp->M[j*Norb+i] = el->Fao_alp->M[i*Norb+j];
  
          }
  
          else if(prms.eht_formula==2){  // Calzaferi formula:        
  
            delt = (el->Fao_alp->M[i*Norb+i] - el->Fao_alp->M[j*Norb+j])/(el->Fao_alp->M[i*Norb+i] + el->Fao_alp->M[j*Norb+j]);
            delt2 = delt*delt;
            delt4 = delt2*delt2;
            
            double rab = (syst.Atoms[a].Atom_RB.rb_cm - syst.Atoms[b].Atom_RB.rb_cm).length();
            double delta = 0.13;
  
  
            std::string elt_i = basis_ao[i].element;
            std::string elt_j = basis_ao[j].element;
            std::string sh_i  = basis_ao[i].ao_shell;
            std::string sh_j  = basis_ao[j].ao_shell;
  
  
            float n_i = modprms.PT[elt_i].Nquant[sh_i];
            float n_j = modprms.PT[elt_j].Nquant[sh_j];
            int nz_i  = modprms.PT[elt_i].Nzeta[sh_i];
            int nz_j  = modprms.PT[elt_j].Nzeta[sh_j];
  
            double ri, rj;  // radii
  
            // for atom i 
            if(nz_i==1){  ri = (n_i/modprms.PT[elt_i].zetas[sh_i][0]); }
            else if(nz_i==2){
              double z1 = modprms.PT[elt_i].zetas[sh_i][0];
              double z2 = modprms.PT[elt_i].zetas[sh_i][1];
              double c1 = modprms.PT[elt_i].coeffs[sh_i][0];
              double c2 = modprms.PT[elt_i].coeffs[sh_i][1];
  
              ri = n_i/(c1*c1*z1 + c2*c2*z2 + 
                         ( pow(2.0,2.0*n_i)*pow(z1*z2, n_i+0.5)/pow((z1+z2),2.0*n_i)
                         ) 
                       ); 
  
            }
  
            // for atom j
            if(nz_j==1){  rj = (n_j/modprms.PT[elt_j].zetas[sh_j][0]); }
            else if(nz_j==2){
              double z1 = modprms.PT[elt_j].zetas[sh_j][0];
              double z2 = modprms.PT[elt_j].zetas[sh_j][1];
              double c1 = modprms.PT[elt_j].coeffs[sh_j][0];
              double c2 = modprms.PT[elt_j].coeffs[sh_j][1];
  
              rj = n_j/(c1*c1*z1 + c2*c2*z2 + 
                         ( pow(2.0,2.0*n_j)*pow(z1*z2, n_j+0.5)/pow((z1+z2),2.0*n_j)
                         ) 
                       ); 
  
            }
  
                                               
            double d0 = ri + rj;    
  
            K_const = 1.0 + (0.75 + delt2 - 0.75*delt4)*exp(-delta*(rab - d0));
          
            el->Fao_alp->M[i*Norb+j] = 0.5*K_const*(el->Fao_alp->M[i*Norb+i]+el->Fao_alp->M[j*Norb+j])*el->Sao->M[i*Norb+j];
            el->Fao_alp->M[j*Norb+i] = el->Fao_alp->M[i*Norb+j];
  
          }
  
          else if(prms.eht_formula==3){       
           // Add your options here and next
  
          }// ==3
       
      }// for bj
    }// for i


    // So far we have computed only Fao_alp
    *el->Fao_bet = *el->Fao_alp;

  }// eht_sce_formula == 1

  else{
    cout<<"Warning (in Hamiltonian_Fock_eht):  no eht_sce_formula="<<prms.eht_sce_formula<<" is known. Skipping...\n";
  }




  //========================= Additional QEq-like terms ==========================  
  if(prms.eht_electrostatics>=1){

    for(a=0;a<syst.Number_of_atoms;a++){  // over all atoms - effects of electronegativities

      // Now add this to Fock matrix:
      int orb_a = atom_to_ao_map[a][0];  // first orbital of the atom a - assume that coefficients do not depend on orbital, only
                                         // on atom type
      double xi_a = modprms.orb_params[orb_a].J_param1; // PT[basis_ao[orb_A].element].J_param1[basis_ao[orb_A].ao_shell]; 
      double J_aa = modprms.orb_params[orb_a].J_param2; // PT[basis_ao[orb_A].element].J_param2[basis_ao[orb_A].ao_shell]; 
      double Qa = syst.Atoms[a].Atom_mull_charge_gross;

    
      for(i=0;i<atom_to_ao_map[a].size();i++){  // over all orbitals on a
        for(j=0;j<atom_to_ao_map[a].size();j++){  // over all orbitals on a
    
  
          //************ INDO-type on-site exchange *************   
          // Exchange effects - only local, to preserve rotational invariance - this is similar to INDO
          double eri_ab = modprms.meht_k.get_K_value(1,i,j); // a.u. of energy - exchange integral
    
          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*eri_ab; 
            el->Fao_bet->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*eri_ab; 
          }
          else{
            el->Fao_alp->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*eri_ab; 
            el->Fao_bet->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*eri_ab; 
          }
    
    

          //************ This corresponds to contribution to energy:  xi_a * Q_a + 1/2*J_aa * Q_a^2 *************
                
          el->Fao_alp->M[i*Norb+j] += (xi_a + J_aa * Qa) * (-el->Sao->M[i*el->Norb+j]);
          el->Fao_bet->M[i*Norb+j] += (xi_a + J_aa * Qa) * (-el->Sao->M[i*el->Norb+j]);

          //*****************************************************************************************************
    
        }// for j
      }// for i
    }// for a

  }// prms.eht_electrostatics>=1
 
  if(prms.eht_electrostatics>=2){ 

    //***************** This corresponds to contribution to energy: 1/2 * J_ab * Q_a * Q_b ****************


    // Contribution from dQ_a/dP_ij

    for(a=0;a<syst.Number_of_atoms;a++){  // over all atoms - effects of electronegativities    
      double Qa = syst.Atoms[a].Atom_mull_charge_gross;
      int orb_a = atom_to_ao_map[a][0];  // first orbital of the atom a - assume that coefficients do not depend on orbital, only
                                         // on atom type

      for(i=0;i<atom_to_ao_map[a].size();i++){  // over all orbitals on a
        for(j=0;j<atom_to_ao_map[a].size();j++){  // over all orbitals on a


          for(b=0;b<syst.Number_of_atoms;b++){  // over all atoms - effects of electronegativities

            if(b!=a){

              int orb_b = atom_to_ao_map[b][0];  // first orbital of the atom b - assume that coefficients do not depend on orbital, only
                                                 // on atom type

              double Qb = syst.Atoms[b].Atom_mull_charge_gross;
              double dist = (syst.Atoms[a].Atom_RB.rb_cm - syst.Atoms[b].Atom_RB.rb_cm).length();

              double K2_const = modprms.meht_k.get_K_value(2,orb_a,orb_b); // a.u. of energy
              double K3_const = modprms.meht_k.get_K_value(3,orb_a,orb_b); // a.u. of length 
              double K4_const = modprms.meht_k.get_K_value(4,orb_a,orb_b); // a.u. of length
              if(K4_const<0.0){  K4_const = 0.0; }

              double J_ab = 1.0/sqrt(dist*dist + K3_const*K3_const);
      
       
              double f = ERF(K4_const * dist);
              J_ab = J_ab * f +  (1.0 - f) * K2_const;  // add long-range electrostatics smoothly
       
          
              el->Fao_alp->M[i*Norb+j] += (J_ab * Qb) * (-el->Sao->M[i*el->Norb+j]);
              el->Fao_bet->M[i*Norb+j] += (J_ab * Qb) * (-el->Sao->M[i*el->Norb+j]);

            }// b!=a          
          }// for b

        }// for j     
      }// for i
    }// for a


    // Contribution from dQ_b/dP_ij

    for(b=0;b<syst.Number_of_atoms;b++){  // over all atoms - effects of electronegativities    
      double Qb = syst.Atoms[b].Atom_mull_charge_gross;
      int orb_b = atom_to_ao_map[b][0];  // first orbital of the atom b - assume that coefficients do not depend on orbital, only
                                         // on atom type

      for(i=0;i<atom_to_ao_map[b].size();i++){  // over all orbitals on b
        for(j=0;j<atom_to_ao_map[b].size();j++){  // over all orbitals on b


          for(a=0;a<syst.Number_of_atoms;a++){  // over all atoms - effects of electronegativities

            if(b!=a){

              int orb_a = atom_to_ao_map[a][0];  // first orbital of the atom a - assume that coefficients do not depend on orbital, only
                                                 // on atom type

              double Qa = syst.Atoms[a].Atom_mull_charge_gross;

              double dist = (syst.Atoms[a].Atom_RB.rb_cm - syst.Atoms[b].Atom_RB.rb_cm).length();

              double K2_const = modprms.meht_k.get_K_value(2,orb_a,orb_b); // a.u. of energy
              double K3_const = modprms.meht_k.get_K_value(3,orb_a,orb_b); // a.u. of length  
              double K4_const = modprms.meht_k.get_K_value(4,orb_a,orb_b); // a.u. of length
              if(K4_const<0.0){  K4_const = 0.0; }

              double J_ab = 1.0/sqrt(dist*dist + K3_const*K3_const);
      
       
              double f = ERF(K4_const * dist);
              J_ab = J_ab * f +  (1.0 - f) * K2_const;  // add long-range electrostatics smoothly
       
          
              el->Fao_alp->M[i*Norb+j] += (J_ab * Qa) * (-el->Sao->M[i*el->Norb+j]);
              el->Fao_bet->M[i*Norb+j] += (J_ab * Qa) * (-el->Sao->M[i*el->Norb+j]);

            }// b!=a          
          }// for a

        }// for j     
      }// for i
    }// for b



  }// if eht_electrostatics == 2


    
}



void Hamiltonian_Fock_eht(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
                           Control_Parameters& prms, Model_Parameters& modprms,
                           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                          ){
/**
  \param[in,out] el The electronic structre properties of the system
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized

  Just the Python-friendly version  
  Compute the EHT Fock Hamiltonian. Well, we are now talking about even further generalized EHT - the one that includes
  self-consistent charges and so. That is why we get density-matrix-dependent Fock matrix.
  Our options:

  prms.eht_sce_formula==0  - charge-independent matrix, so F = H (core)
  prms.eht_sce_formula==1  - we first add orbital energy correction that is linearly-dependent on atomic Mulliken charges, then
                             the corrected orbital energies are used in one of the mixing formulas:
                             
    prms.eht_formula == 0 - unweighted 
    prms.eht_formula == 1 - weighted
    prms.eht_formula == 2 - Calzaferi
    prms.eht_formula == 3 - for the developments

  On top of this all, we also have a QE-type of correction:

  prms.eht_electrostatics>=1

    a) this gives on-site Fock matrix correction of the INDO type (breaking spin symmetry, in general) - this is actually controlled by the
    model parameters

    b) also, this correction gives the Fock matrix terms that originate from the energy term: xi_a * Q_a + 1/2*J_aa * Q_a^2 -on-site energies

  prms.eht_electrostatics>=2

    same as   prms.eht_electrostatics>=1 but also the terms originating from the energy terms 1/2 * J_ab * Q_a * Q_b are added to the 
    Fock matrix

*/


  Hamiltonian_Fock_eht(&el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

}





}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra
