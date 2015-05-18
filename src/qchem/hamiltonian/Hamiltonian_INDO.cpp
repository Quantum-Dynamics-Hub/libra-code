#include "Hamiltonian.h"

/****************************************************************************

  This file contains following functions:

  void indo_core_parameters(vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao, Nuclear& mol,
                            vector<double>& eri, vector<double>& V_AB, Memory* mem, int opt)


  void Hamiltonian_core_indo(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                             vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao, vector<vector<int> >& at_orbitals, 
                             MATRIX* Hao,MATRIX* Sao,Memory* mem, vector<double>& eri, vector<double>& V_AB)



  void get_integrals(int i,int j,vector<AO>& basis_ao, double eri_aa, double G1, double F2, double& ii_jj,double& ij_ij)

  void Hamiltonian_Fock_indo(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                             vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                             Electronic* el, Memory* mem,vector<double>& eri, vector<double>& V_AB)



****************************************************************************/

void indo_core_parameters(vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao, Nuclear& mol,
                          vector<double>& eri, vector<double>& V_AB, Memory* mem, int opt){
// opt == 0 - cndo
// opt == 1 - indo

  int DF = 0;

  int i,j,a,b,i1,a1,I,J,A;
  VECTOR da,db,dc;


  if(DF){ cout<<"in indo_core_parameters\n"; }

  // Allocate memory if needed
  int sz = fragment.size(); // number of atoms in the group <fragment>

  if(eri.size()!=sz*sz){  cout<<"In indo_core_parameters: eri array is not allocated\nDo allocation...\n"; 
    eri.clear(); eri = vector<double>(sz*sz,0.0); 
  }
  if(V_AB.size()!=sz*sz){  cout<<"In indo_core_parameters: V_AB array is not allocated\nDo allocation...\n"; 
    V_AB.clear(); V_AB = vector<double>(sz*sz,0.0); 
  }


  // Compute global indices of s-type orbitals on each atom
  vector<int> sorb_indx(sz,0); // global index of s-type orbital on i-th atom

  for(i=0;i<sz;i++){
    a = fragment[i];  // global index of i-th atom in the group <fragment>

    for(i1=0;i1<basis_fo.size();i1++){ // all orbitals in the group <fragment>
      I = basis_fo[i1]; // global index of i1-th orbital in the fragment

      if(basis_ao[I].at_indx==i && basis_ao[I].ao_shell_type=="s" ){  sorb_indx[i] = I; }
          
    }// for i1
  }// for i

  // Printing mapping
  if(DF){ 
    cout<<"i - runs over local indices of atoms in given fragment\n";
    cout<<"sorb_indx[i] - is the global index of s-type orbital (assuming only one) centered on atom i\n";
    for(i=0;i<sorb_indx.size();i++){
      cout<<"i= "<<i<<" sorb_indx[i]= "<<sorb_indx[i]<<endl;
    }
  }
    
  
  // Compute ERIs and V_AB
  for(a=0;a<sz;a++){
    for(b=0;b<sz;b++){

      I = sorb_indx[a];
      J = sorb_indx[b];

     
      // ERI
      eri[a*sz+b] = ELECTRON_REPULSION_INTEGRAL(basis_ao[I],basis_ao[I],basis_ao[J],basis_ao[J]);


      // V_AB
      int B = fragment[b]; // global index of atom b
      if(opt==0){
        V_AB[a*sz+b] = mol.Zeff[B]*NUCLEAR_ATTRACTION_INTEGRAL(basis_ao[J],basis_ao[J],mol.R[B],B,da,db,dc,mem->aux,mem->n_aux,mem->auxv,mem->n_auxv);// V_AB[a][b]
      }
      else if(opt==1){
        V_AB[a*sz+b] = mol.Zeff[B]*eri[a*sz+b];
      }

      if(DF){ cout<<"a= "<<a<<" b= "<<b<<" I= "<<I<<" J= "<<J<<" eri= "<<eri[a*sz+b]<<" V_AB= "<<V_AB[a*sz+b]<<endl; }

    }// for j
  }// for i


}



void Hamiltonian_core_indo(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                           vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao, vector<vector<int> >& at_orbitals, 
                           MATRIX* Hao,MATRIX* Sao,Memory* mem, vector<double>& eri, vector<double>& V_AB){

  int DF = 0;

  int i,j,k,a,b,I,J,A,B;
  VECTOR dIdA,dIdB;
  //cout<<"in Hamiltonian_core_indo\n";

  int Norb = basis_fo.size(); // how many AOs are included in this fragment
  if(Norb!=Hao->num_of_cols){  
    cout<<"In Hamiltonian_core_eht: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  int sz = fragment.size(); // number of atoms in this fragment

  if(eri.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_indo: size of auxiliary eri array is not right\n"; exit(0);}
  if(V_AB.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_indo: size of auxiliary V_AB array is not right\n"; exit(0);}



  // Compute global indices of s-type orbitals on each atom
  vector<int> sorb_indx(sz,0); // global index of s-type orbital on i-th atom

  for(i=0;i<sz;i++){
    a = fragment[i];  // global index of i-th atom in the group <fragment>

    for(int i1=0;i1<basis_fo.size();i1++){ // all orbitals in the group <fragment>
      I = basis_fo[i1]; // global index of i1-th orbital in the fragment

      if(basis_ao[I].at_indx==i && basis_ao[I].ao_shell_type=="s" ){  sorb_indx[i] = I; }
          
    }// for i1
  }// for i




  //----------- Compute core terms of the Hamiltonian ------------
  *Hao = 0.0;



  for(i=0;i<Norb;i++){
    I = basis_fo[i];         ///< global index of the i-th AO in this fragment
    A = basis_ao[I].at_indx; ///< global index of the atom on which I-th AO is localized
    a = -1;                  ///< local index of atom A in the fragment
    for(j=0;j<sz;j++){   if(fragment[j] == A) { a = j; }    }


    // values of IP are different for cndo (just IPs) and cndo2 and indo ( 0.5*(IP + EA) )
    Hao->M[i*Norb+i] += modprms.PT[basis_ao[I].element].IP[basis_ao[I].ao_shell];

    if(DF){
      cout<<"Setting diagonal element i = "<<i<<" I= "<<I<<endl;
      cout<<"Contribution from IP = "<<Hao->M[i*Norb+i]<<endl;
    }


    /// The code below is the same for CNDO, CNDO2 and INDO - but the difference comes in use of different G1 and F2 parameters
    /// for CNDO and CNDO2 they are zero
    double G1 = modprms.PT[basis_ao[I].element].G1[basis_ao[I].ao_shell];
    double F2 = modprms.PT[basis_ao[I].element].F2[basis_ao[I].ao_shell];
    
    if(DF){
      cout<<"A= "<<A<<" G1= "<<G1<<" F2= "<<F2<<" mol.Z[A]= "<<mol.Z[A]<<" basis_ao[I].ao_shell_type= "<<basis_ao[I].ao_shell_type<<endl;
    }
   
    /// Eqs. 3.17 - 3.23 from Pople, Beveridge, Dobosh, JCP 47, 2026 (1967)
    ///
         if(mol.Z[A]==1){ Hao->M[i*Norb+i] -= 0.5*eri[a*sz+a]; }  // H
    else if(mol.Z[A]==3 || mol.Z[A]==11){ 

      if(basis_ao[I].ao_shell_type=="s"){
        Hao->M[i*Norb+i] -= 0.5*eri[a*sz+a]; // s
      }
      else if(basis_ao[I].ao_shell_type=="p"){
        Hao->M[i*Norb+i] -= (0.5*eri[a*sz+a] - G1/12.0); // p
      }
    }  // Li or Na

    else if(mol.Z[A]==4 || mol.Z[A]==12){ 

      if(basis_ao[I].ao_shell_type=="s"){
        Hao->M[i*Norb+i] -= (1.5*eri[a*sz+a] - 0.5*G1); // s
      }
      else if(basis_ao[I].ao_shell_type=="p"){
        Hao->M[i*Norb+i] -= (1.5*eri[a*sz+a] - 0.25*G1); // p
      }
    }  // Be or Mg

    else if( (mol.Z[A]>=5 && mol.Z[A]<=9) || (mol.Z[A]>=13 && mol.Z[A]<=18)){ // B - F or Al - Cl

      double Z = modprms.PT[basis_ao[I].element].Nval; // core charge of atom 

      if(basis_ao[I].ao_shell_type=="s"){      
        Hao->M[i*Norb+i] -= ( ( Z - 0.5 )*eri[a*sz+a] - (1.0/6.0)*( Z - 1.5 )*G1); // s
      }
      else if(basis_ao[I].ao_shell_type=="p"){
        Hao->M[i*Norb+i] -= ( ( Z - 0.5 )*eri[a*sz+a] - (1.0/3.0)*G1 - 0.08*( Z - 2.5 )*F2); // p
      }

    }
    else{  cout<<"Error: INDO is not implemented for elements beyond Cl\n"; exit(0);  }

    if(DF){
      cout<<" + Contribution from Frank-Condon factors = "<<Hao->M[i*Norb+i]<<endl;
    }
    
    //----------------- Coulombic terms --------------

    for(b=0;b<sz;b++){
      if(b!=a){
        Hao->M[i*Norb+i] -= V_AB[a*sz+b];  //  = V_AB = Z_B * eri[A][B] -in INDO
      }
    }// for b

    if(DF){
      cout<<" + Contribution from Coulombic terms = "<<Hao->M[i*Norb+i]<<endl;
    }


    //-------------- Off-diagonal terms of the core matrix ---------
    for(j=0;j<Norb;j++){

      if(j!=i){
        J = basis_fo[j];         ///< global index of the j-th AO in this fragment
        B = basis_ao[J].at_indx; ///< global index of the atom on which J-th AO is localized
        b = -1;                  ///< local index of atom B in the fragment
        for(k=0;k<sz;k++){  if(fragment[k] == B) { b = k; }     }


        if(B==A){ ;; }  // different orbitals centered on the same atom - give zero (not true for hybrid orbitals)
        else{           // centered on different atoms - use overlap formula

          // Overlap is set to identity in INDO, so need to recompute it explicitly
          double sao_ij = OVERLAP_INTEGRAL(basis_ao[I],basis_ao[J],0,dIdA,dIdB,mem->aux,mem->n_aux);

          Hao->M[i*Norb+j] += 0.5*(modprms.PT[basis_ao[I].element].beta0[basis_ao[I].ao_shell] 
                                 + modprms.PT[basis_ao[J].element].beta0[basis_ao[J].ao_shell]) * sao_ij;
        }
      }// j!=i

    }// for j    

  }// for i



}



void get_integrals(int i,int j,vector<AO>& basis_ao, double eri_aa, double G1, double F2, double& ii_jj,double& ij_ij){

  // Compute Coulomb and exchange integrals for the orbitals i and j (global indices), both of which 
  // are centered on the same atom a (global index)
  // eri[a][a] is taken as input argument
  // parameters G1 and F2 (Slater-Condon) are taken as input

  //=====================================================================================================
  // Integrals:
  ij_ij = ii_jj = 0.0; 

  if( basis_ao[i].ao_shell_type=="s" && basis_ao[j].ao_shell_type=="s"){ 
    ij_ij = ii_jj = eri_aa;                  // ss_ss 
  }

  else if( basis_ao[i].ao_shell_type=="s" && basis_ao[j].ao_shell_type=="p"){
    ij_ij = G1/3.0;                          // sx_sx = sy_sy = sz_sz
    ii_jj = eri_aa;                          // ss_xx = ss_yy = ss_zz
  }
  else if( basis_ao[j].ao_shell_type=="s" && basis_ao[i].ao_shell_type=="p"){
    ij_ij = G1/3.0;                          // xs_xs = ys_ys = zs_zs
    ii_jj = eri_aa;                          // xx_ss = yy_ss = zz_ss
  }
  else if( basis_ao[i].ao_shell_type=="p" && basis_ao[j].ao_shell_type=="p" ){ 
    if( (basis_ao[i].x_exp == basis_ao[j].x_exp) && 
        (basis_ao[i].y_exp == basis_ao[j].y_exp) && 
        (basis_ao[i].x_exp == basis_ao[j].x_exp)){
      ij_ij = ii_jj = eri_aa + 0.16*F2;      // xx_xx = yy_yy = zz_zz
    }
    else{
      ij_ij = 0.12*F2;                         // xy_xy = xz_xz = ...
      ii_jj = eri_aa - 0.08*F2;                // xx_yy = xx_zz = ...
    }
  }

  //=====================================================================================================

}



void Hamiltonian_Fock_indo(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                           vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                           Electronic* el, Memory* mem,vector<double>& eri, vector<double>& V_AB){

  
// This function constructs INDO (CNDO/2) Fock matrix
// Unrestricted formulation

  int i,j,k,n,I,J,K,a,b,A,B;

  int Norb = basis_fo.size(); // how many AOs are included in this fragment
  if(Norb!=el->Hao->num_of_cols){  
    cout<<"In Hamiltonian_Fock_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }


  // Formation of the Fock matrix: Core part
  *el->Fao_alp = *el->Hao;
  *el->Fao_bet = *el->Hao;


  // Update charges
  *el->P = *el->P_alp + *el->P_bet;

  update_Mull_orb_pop(el->P, el->Sao, el->Mull_orb_pop_gross, el->Mull_orb_pop_net);


  // Here we need a local version of at_orbitals - only those atoms that carry orbitals spanning basis_fo
  // we assume that there is no overlaps of atoms between two subsystems
  update_Mull_charges(fragment, basis_fo, at_orbitals, mol.Zeff, el->Mull_orb_pop_gross, el->Mull_orb_pop_net,mol.Mull_charges_gross, mol.Mull_charges_net);


    
  // Formation of the Fock matrix: add Coulomb and Exchange parts    
  for(i=0;i<Norb;i++){
    I = basis_fo[i];
    A = basis_ao[I].at_indx;  // basis_fo[i] on atom a - global atom index
    a = -1;                   // local index of atom A
    for(n=0;n<fragment.size();n++){ if(fragment[n]==A){ a = n; }  } 


    for(j=0;j<Norb;j++){
      J = basis_fo[j];
      B = basis_ao[J].at_indx;
      b = -1;                   // local index of atom B
      for(n=0;n<fragment.size();n++){ if(fragment[n]==B){ b = n; }  } 


      if(I==J){  // Diagonal terms


        for(int kk=0;kk<at_orbitals[A].size();kk++){                       // for all orbitals on atom A
          K = at_orbitals[A][kk];                                          // global orbital index
          k = -1;  for(n=0;n<Norb;n++){ if(basis_fo[n] == K) { k = n; } }  // local orbital index


          double ii_kk, ik_ik; ii_kk = ik_ik = 0.0;
          double G1 = modprms.PT[basis_ao[I].element].G1[basis_ao[I].ao_shell];
          double F2 = modprms.PT[basis_ao[I].element].F2[basis_ao[I].ao_shell];
          double eri_aa = eri[a*fragment.size()+a];

          get_integrals(I,K,basis_ao,eri_aa,G1,F2,ii_kk,ik_ik);

          if(prms.use_rosh){ // Restricted open-shell
            el->Fao_alp->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - 0.5*el->P->M[k*Norb+k]*ik_ik);
            el->Fao_bet->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - 0.5*el->P->M[k*Norb+k]*ik_ik);

          }
          else{ // unrestricted
            el->Fao_alp->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - el->P_alp->M[k*Norb+k]*ik_ik);
            el->Fao_bet->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - el->P_bet->M[k*Norb+k]*ik_ik);
          }
  
        }// for kk - all orbitals on atom A
  

        // Contributions from all other atoms to the diagonal terms
        // don't worry that b determined above will be rewritten - this is ok
        for(b=0;b<fragment.size();b++){
          B = fragment[b];

          if(B!=A){

            // Compute density matrix due to all orbitals on atom b 
            el->Fao_alp->M[i*Norb+i] += (mol.Zeff[B] - mol.Mull_charges_net[B])*eri[a*fragment.size()+b]; 
            el->Fao_bet->M[i*Norb+i] += (mol.Zeff[B] - mol.Mull_charges_net[B])*eri[a*fragment.size()+b]; 
          }
        }// for b
  
  
      }// i==j
      else{      // Off-diagonal terms
  
        if(A==B){ // different orbitals are on the same atom
          double ij_ij,ii_jj; ij_ij = ii_jj = 0.0;

          double G1 = modprms.PT[basis_ao[I].element].G1[basis_ao[I].ao_shell];
          double F2 = modprms.PT[basis_ao[I].element].F2[basis_ao[I].ao_shell];
          double eri_aa = eri[a*fragment.size()+a];
          get_integrals(I,J,basis_ao,eri_aa,G1,F2,ii_jj,ij_ij);

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

          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*eri[a*fragment.size()+b]; 
            el->Fao_bet->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*eri[a*fragment.size()+b]; 
          }
          else{
            el->Fao_alp->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*eri[a*fragment.size()+b]; 
            el->Fao_bet->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*eri[a*fragment.size()+b]; 
          }

        }
  
  
      }// i!=j
  
  
    }// for j
  }// for i


}
