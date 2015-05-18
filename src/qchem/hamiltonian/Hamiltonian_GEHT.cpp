#include "Hamiltonian.h"

/****************************************************************************

  This file contains following functions:

  void Hamiltonian_core_geht(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                             vector<int>& fragment,vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                             MATRIX* Hao,MATRIX* Sao,Memory* mem)

  void Hamiltonian_Fock_geht(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                             vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                             Electronic* el,Electronic* el0, Memory* mem)




****************************************************************************/

void Hamiltonian_core_geht(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                          vector<int>& fragment,vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                          MATRIX* Hao,MATRIX* Sao,Memory* mem){

// Generalized EHT

// prms - control parameters
// modprms - model parameters
// mol - nuclear geometry and charges
// basis_fo - basis of fragment obitals - basis_fo[i] - is index of basis_ao, which is i-th in the given group
// basis_ao - the full/complete pool of AOs (entire system)
// Hao, Sao - Hamiltonian and overlap matrices


           
  int i,j,n, I,J;
  double delt, delt2, delt4;
  

  int Norb = basis_fo.size(); // how many AOs are included in this fragment
  if(Norb!=Hao->num_of_cols){  
    cout<<"In Hamiltonian_core_eht: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
 
  //========================= Core diagonal elements ==========================  
  // Diagonal elements = set to IPs
  for(i=0;i<Norb;i++){  
    I = basis_fo[i];  
    Hao->M[i*Norb+i] = modprms.PT[basis_ao[I].element].IP[basis_ao[I].ao_shell];
  }// for i


  //========================= Off-diagonal elements ==========================    
  double K_const = 0.0;// = 1.75;

  // Off-diagonal elements
  for(i=0;i<Norb;i++){
    I = basis_fo[i];
    for(j=i+1;j<Norb;j++){
        J = basis_fo[j];

//        K_const = modprms.eht_k.get_K_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);
        K_const = modprms.meht_k.get_K_value(I,J);


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
          
          double rab = (mol.R[basis_ao[i].at_indx] - mol.R[basis_ao[j].at_indx]).length();
          double delta = 0.13;
          double d0 =  (modprms.PT[basis_ao[i].element].PQN/modprms.PT[basis_ao[i].element].Zeta + 
                       modprms.PT[basis_ao[j].element].PQN/modprms.PT[basis_ao[j].element].Zeta
                      );    


          K_const = 1.0 + (0.75 + delt2 - 0.75*delt4)*exp(-delta*(rab - d0));
        
          Hao->M[i*Norb+j] = 0.5*K_const*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j])*Sao->M[i*Norb+j];
          Hao->M[j*Norb+i] = Hao->M[i*Norb+j];

        }

        else if(prms.eht_formula==3){       // This is slightly generalized eht_formula == 1

//          double K1_const = modprms.eht_k.get_K1_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);
//          double K2_const = modprms.eht_k.get_K2_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);

          double K1_const = modprms.meht_k.get_K1_value(I,J);
          double K2_const = modprms.meht_k.get_K2_value(I,J);


          Hao->M[i*Norb+j]  = (0.5*K_const*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j]) + K1_const)*Sao->M[i*Norb+j];        


          Hao->M[j*Norb+i] = Hao->M[i*Norb+j];

        }// == 3
     
    }// for bj
  }// for i



}


void Hamiltonian_Fock_geht(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                           vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                           Electronic* el,Electronic* el0, Memory* mem){


  int Norb = basis_fo.size(); // how many AOs are included in this fragment
  if(Norb!=el->Hao->num_of_cols){  
    cout<<"In Hamiltonian_Fock_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }

//  cout<<"In Hamiltonian_Fock_geht\n";


  // Update charges
  *el->P = *el->P_alp + *el->P_bet;  // just in case - this is supposed to be already performed outside

  update_Mull_orb_pop(el->P, el->Sao, el->Mull_orb_pop_gross, el->Mull_orb_pop_net);

  // Here we need a local version of at_orbitals - only those atoms that carry orbitals spanning basis_fo
  // we assume that there is no overlaps of atoms between two subsystems

  update_Mull_charges(fragment, basis_fo, at_orbitals, mol.Zeff, el->Mull_orb_pop_gross, el->Mull_orb_pop_net,mol.Mull_charges_gross, mol.Mull_charges_net);


  // Recompute (core) EHT Hamiltonian - it may be charge-dependent (well, this is not strictly speaking a "core")
  // Now we don't need it, because it does not depend on charge
  //Hamiltonian_core_eht(prms,modprms,mol,fragment,basis_fo,basis_ao, at_orbitals, el->Hao,el->Sao,mem);


  // Compute Fock matrices
  // Core contribution
  *el->Fao_alp = *el->Hao;
  *el->Fao_bet = *el->Hao;


  // Diagonal electrostatics
  vector<double> Wii(Norb,0.0);

//  cout<<"Testing Fock_geht\n";

  // Precompute electrostatic correction factor for all orbitals
  for(int i=0;i<Norb;i++){
    int I = basis_fo[i];

    Wii[i] = 0.0;

    for(int j=0;j<Norb;j++){
      int J = basis_fo[j];

//      double K3_const_old = modprms.eht_k.get_K3_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);
//      double K4_const_old = modprms.eht_k.get_K4_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);

      double K3_const = modprms.meht_k.get_K3_value(I,J);
      double K4_const = modprms.meht_k.get_K4_value(I,J);

//      cout<<"I= "<<I<<" J= "<<J<<" K3_old= "<<K3_const_old<<" K3_new= "<<K3_const<<" K4_old= "<<K4_const_old<<" K4_new= "<<K4_const<<endl;
//      cout<<"I= "<<I<<" J= "<<J<<" K3_new= "<<K3_const<<" K4_new= "<<K4_const<<endl;


      if(basis_ao[J].at_indx!=basis_ao[I].at_indx){  // orbitals on all other atoms

     
        double dist2 = (mol.R[basis_ao[I].at_indx]-mol.R[basis_ao[J].at_indx]).length2();

        // This approximates Q[n] * ERI[n,I]
        Wii[i] -= K3_const*K3_const * mol.Mull_charges_gross[basis_ao[J].at_indx]/sqrt((K4_const*K4_const) + dist2);

      }// if - different atoms


      // Contribution from periodic images - well, this is not really converging, but at least first environment will be
      // taken into account
      if(prms.x_period>=1 || prms.y_period>=1 || prms.z_period>=1){

        for(int nx=-prms.x_period;nx<=prms.x_period;nx++){
          for(int ny=-prms.y_period;ny<=prms.y_period;ny++){
            for(int nz=-prms.z_period;nz<=prms.z_period;nz++){
        
              if(nx==0 && ny==0 && nz==0){ }
              else{
                
                VECTOR TV = nx*prms.t1 + ny*prms.t2 + nz*prms.t3;
     
                double dist2 = (mol.R[basis_ao[I].at_indx]-mol.R[basis_ao[J].at_indx]-TV).length();
                // This summation corresponds to k = 0 (Gamma-point)
     
                // This approximates - Q[n] * eri[]
                Wii[i] -= K3_const*K3_const*mol.Mull_charges_gross[basis_ao[J].at_indx]/sqrt((K4_const*K4_const) + dist2);
     
              }
        
            }// for nz
          }// for ny
        }// for nx 

      }// if - periodic



    }// for j
  }// for i


  // Now use diagonal correction factors to update all elements of Fock matrix
  for(int i=0;i<Norb;i++){

    el->Fao_alp->M[i*Norb+i] += Wii[i]; //0.5*(Wii[i]+Wii[j])*el->Sao->M[i*el->Norb+j];
    el->Fao_bet->M[i*Norb+i] += Wii[i]; //0.5*(Wii[i]+Wii[j])*el->Sao->M[i*el->Norb+j];

  }// i



  // Compute mapping from globar orbital index K to local orbital index k (local for this fragment)
  vector<int> loc_indx(el->Norb, -1);  // dimension of this array is the full dimension 
  
  for(int n=0;n<Norb;n++){ // loop over fragment orbitals
    loc_indx[ basis_fo[n] ] = n;  
  }



  // This is generalization of INDO Fock
  // Formation of the Fock matrix: add Coulomb and Exchange parts    
  for(int i=0;i<Norb;i++){
    int I = basis_fo[i];
    int A = basis_ao[I].at_indx;  // basis_fo[i] on atom a - global atom index

//    int a = -1;                   // local index of atom A
//    for(int n=0;n<fragment.size();n++){ if(fragment[n]==A){ a = n; }  } 


    for(int j=0;j<Norb;j++){
      int J = basis_fo[j];
      int B = basis_ao[J].at_indx;

//      int b = -1;                   // local index of atom B
//      for(int n=0;n<fragment.size();n++){ if(fragment[n]==B){ b = n; }  } 

 
      // On-site Coulomb and exchange
      if(I==J){  // Diagonal terms  - same orbital

        for(int kk=0;kk<at_orbitals[A].size();kk++){                               // for all orbitals on atom A
          int K = at_orbitals[A][kk];                                              // global orbital index
          int k = loc_indx[K];  //for(int n=0;n<Norb;n++){ if(basis_fo[n] == K) { k = n; } }  // local orbital index
       

          double ii_kk, ik_ik; ii_kk = ik_ik = 0.0;

          // Integrals should be in groups: [s-s],  [s-p, p-p], [s-p,s-d, p-d,p-p,d-d]          
//          double ii_kk_old = modprms.eht_k.get_K1_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[K].element,basis_ao[K].ao_shell);
//          double ik_ik_old = modprms.eht_k.get_K2_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[K].element,basis_ao[K].ao_shell);

          ii_kk = modprms.meht_k.get_K1_value(I,K);
          ik_ik = modprms.meht_k.get_K2_value(I,K);

//      cout<<"I= "<<I<<" K= "<<K<<" ii_kk_old= "<<ii_kk_old<<" ii_kk_new= "<<ii_kk<<" ik_ik_old= "<<ik_ik_old<<" ik_ik_new= "<<ik_ik<<endl;



          if(prms.use_rosh){ // Restricted open-shell
            el->Fao_alp->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - 0.5*el->P->M[k*Norb+k]*ik_ik);
            el->Fao_bet->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - 0.5*el->P->M[k*Norb+k]*ik_ik);

          }
          else{ // unrestricted
            el->Fao_alp->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - el->P_alp->M[k*Norb+k]*ik_ik);
            el->Fao_bet->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - el->P_bet->M[k*Norb+k]*ik_ik);
          }
  
        }// for kk - all orbitals on atom A
   
  
      }// I==J
      else{      // Off-diagonal terms

        // On-site Coulomb and exchange  
        if(A==B){ // different orbitals are on the same atom
          double ij_ij,ii_jj; ij_ij = ii_jj = 0.0;


          // Integrals should be in groups: [s-s],  [s-p, p-p], [s-p,s-d, p-d,p-p,d-d]
//          ii_jj = modprms.eht_k.get_K1_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);
//          ij_ij = modprms.eht_k.get_K2_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);

          ii_jj = modprms.meht_k.get_K1_value(I,J);
          ij_ij = modprms.meht_k.get_K2_value(I,J);



          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - 0.5*el->P->M[i*Norb+j])*ij_ij - 0.5*el->P->M[i*Norb+j]*ii_jj );
            el->Fao_bet->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - 0.5*el->P->M[i*Norb+j])*ij_ij - 0.5*el->P->M[i*Norb+j]*ii_jj );
          }
          else{  
            el->Fao_alp->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - el->P_alp->M[i*Norb+j])*ij_ij - el->P_alp->M[i*Norb+j]*ii_jj );
            el->Fao_bet->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - el->P_bet->M[i*Norb+j])*ij_ij - el->P_bet->M[i*Norb+j]*ii_jj );
          }
  
        }// A==B - different orbitals are on the same atom

  
        else{ // different orbitals are on different atoms

          // This approximates ERI[I,J] - this should depend only on atom types, not on orbitals
//          double K3_const = modprms.eht_k.get_K3_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);
//          double K4_const = modprms.eht_k.get_K4_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);

          double K3_const = modprms.meht_k.get_K3_value(I,J);
          double K4_const = modprms.meht_k.get_K4_value(I,J);


          double dist = (mol.R[A]-mol.R[B]).length();

          double eri_ab = K3_const*K3_const/sqrt((K4_const*K4_const) + dist*dist);


          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*eri_ab; 
            el->Fao_bet->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*eri_ab; 
          }
          else{
            el->Fao_alp->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*eri_ab; 
            el->Fao_bet->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*eri_ab; 
          }

//          cout<<"Different orbital on different atoms part: \n";
//          cout<<basis_ao[I].element<<" "<<basis_ao[J].element<<" K3 = "<<K3_const/eV<<" K4= "<<K4_const<<"  dist= "<<dist<<" eri_ab= "<<eri_ab/eV<<endl;

        }// different orbitals on different atoms
    
      }// I!=J


  
  
    }// for j
  }// for i

//  cout<<"end of Hamiltonian_Fock_geht\n";

    
}

