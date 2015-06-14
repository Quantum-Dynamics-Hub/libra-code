#include "Hamiltonian.h"

/****************************************************************************

  This file contains following functions:

  void Hamiltonian_core_hf(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                           vector<int>& fragment,vector<int>& basis_fo,vector<AO>& basis_ao, vector<vector<int> >& at_orbitals,
                           MATRIX* Hao,MATRIX* Sao,Memory* mem)

  void Hamiltonian_Fock_hf(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                           vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                           Electronic* el,Electronic* el0, Memory* mem)


****************************************************************************/

void Hamiltonian_core_hf(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                         vector<int>& fragment,vector<int>& basis_fo,vector<AO>& basis_ao, vector<vector<int> >& at_orbitals,
                         MATRIX* Hao,MATRIX* Sao,Memory* mem){
// prms - control parameters
// modprms - model parameters
// mol - nuclear geometry and charges
// basis_fo - basis of fragment obitals - basis_fo[i] - is index of basis_ao, which is i-th in the given group
// basis_ao - the full/complete pool of AOs (entire system)
// Hao, Sao - Hamiltonian and overlap matrices

  int i,j,n, I,J;
  VECTOR da,db,dc;
  
  int Norb = basis_fo.size(); // how many AOs are included in this fragment
  if(Norb!=Hao->num_of_cols){  
    cout<<"In Hamiltonian_core_hf: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }


  for(i=0;i<Norb;i++){
    I = basis_fo[i];  
    for(j=0;j<Norb;j++){
      J = basis_fo[j];  

      Hao->M[i*Norb+j] = KINETIC_INTEGRAL(basis_ao[I],basis_ao[J],da,db);

      for(n=0;n<mol.Nnucl;n++){
        Hao->M[i*Norb+j] -= mol.Zeff[n]*NUCLEAR_ATTRACTION_INTEGRAL(basis_ao[I],basis_ao[J],mol.R[n],n,da,db,dc);
      }// for n

    }// for j
  }// for i


}

void Hamiltonian_Fock_hf(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                          vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                          Electronic* el,Electronic* el0, Memory* mem){
// This functions constructs HF Fock matrix

  int a,b,c,d,A,B,C,D;

  int Norb = basis_fo.size(); // how many AOs are included in this fragment
  if(Norb!=el->Hao->num_of_cols){  
    cout<<"In Hamiltonian_Fock_hf: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }


  // Update total density matrix
  *el->P = *el->P_alp + *el->P_bet;


  update_Mull_orb_pop(el->P, el->Sao, el->Mull_orb_pop_gross, el->Mull_orb_pop_net);

  // Here we need a local version of at_orbitals - only those atoms that carry orbitals spanning basis_fo
  // we assume that there is no overlaps of atoms between two subsystems

  update_Mull_charges(fragment, basis_fo, at_orbitals, mol.Zeff, el->Mull_orb_pop_gross, el->Mull_orb_pop_net,mol.Mull_charges_gross, mol.Mull_charges_net);


  // Compute Fock matrices: Core part
  *el->Fao_alp = *el->Hao;
  *el->Fao_bet = *el->Hao;



  // Formation of the Fock matrix: add Coulomb and Exchange parts
  for(a=0;a<Norb;a++){
    A = basis_fo[a];  
    for(b=0;b<Norb;b++){
      B = basis_fo[b];  
      for(c=0;c<Norb;c++){
        C = basis_fo[c];  
        for(d=0;d<Norb;d++){
          D = basis_fo[d];  

          //  (P_cd * (ab|cd) - P_alp_cd*(ad|cb))
//          double J_abcd = ELECTRON_REPULSION_INTEGRAL(basis_ao[A],basis_ao[B],basis_ao[C],basis_ao[D]);
//          double K_adcb = ELECTRON_REPULSION_INTEGRAL(basis_ao[A],basis_ao[D],basis_ao[C],basis_ao[B]);
          double J_abcd,K_adcb;
          modprms.hf_int.get_JK_values(A,B,C,D,J_abcd,K_adcb);

          if(prms.use_rosh){
            el->Fao_alp->M[a*Norb+b] += (el->P->M[c*Norb+d]*J_abcd - 0.5*el->P->M[c*Norb+d]*K_adcb);
            el->Fao_bet->M[a*Norb+b] += (el->P->M[c*Norb+d]*J_abcd - 0.5*el->P->M[c*Norb+d]*K_adcb);
          }
          else{
            el->Fao_alp->M[a*Norb+b] += (el->P->M[c*Norb+d]*J_abcd - el->P_alp->M[c*Norb+d]*K_adcb);
            el->Fao_bet->M[a*Norb+b] += (el->P->M[c*Norb+d]*J_abcd - el->P_bet->M[c*Norb+d]*K_adcb);
          }


        }// for d
      }// for c
    }// for b
  }// for a

}

