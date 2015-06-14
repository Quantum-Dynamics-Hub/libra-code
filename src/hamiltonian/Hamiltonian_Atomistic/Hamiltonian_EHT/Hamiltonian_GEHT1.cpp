#include "Hamiltonian.h"

/****************************************************************************

  This file contains following functions:

  void geht1_core_parameters(vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao, Nuclear& mol,
                             vector<double>& eri, vector<double>& V_AB, Memory* mem, int opt)


  void Hamiltonian_core_geht1(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                              vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao, vector<vector<int> >& at_orbitals, 
                              MATRIX* Hao,MATRIX* Sao,Memory* mem, vector<double>& eri, vector<double>& V_AB)



  void get_geht1_integrals(int i,int j,vector<AO>& basis_ao, double eri_aa, double G1, double F2, double& ii_jj,double& ij_ij)

  void Hamiltonian_Fock_geht1(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                              vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                              Electronic* el, Memory* mem,vector<double>& eri, vector<double>& V_AB)



****************************************************************************/

void geht1_core_parameters(vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao, Nuclear& mol,
                           vector<double>& eri, vector<double>& V_AB, Memory* mem, int opt){
// opt == 0 - cndo
// opt == 1 - indo

  int DF = 0;

  int i,j,a,b,i1,a1,I,J,A;
  VECTOR da,db,dc;


  if(DF){ cout<<"in geht1_core_parameters\n"; }

  // Allocate memory if needed
  int sz = fragment.size(); // number of atoms in the group <fragment>

  if(eri.size()!=sz*sz){  cout<<"In geht1_core_parameters: eri array is not allocated\nDo allocation...\n"; 
    eri.clear(); eri = vector<double>(sz*sz,0.0); 
  }
  if(V_AB.size()!=sz*sz){  cout<<"In geht1_core_parameters: V_AB array is not allocated\nDo allocation...\n"; 
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


/*
void Hamiltonian_core_geht1(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                            vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao, vector<vector<int> >& at_orbitals, 
                            MATRIX* Hao,MATRIX* Sao,Memory* mem, vector<double>& eri, vector<double>& V_AB){

  int DF = 0;

  int i,j,k,a,b,I,J,A,B;
  VECTOR dIdA,dIdB;


  int Norb = basis_fo.size(); // how many AOs are included in this fragment
  if(Norb!=Hao->num_of_cols){  
    cout<<"In Hamiltonian_core_geht1: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  int sz = fragment.size(); // number of atoms in this fragment

  if(eri.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_geht1: size of auxiliary eri array is not right\n"; exit(0);}
  if(V_AB.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_geht1: size of auxiliary V_AB array is not right\n"; exit(0);}



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

*/

void get_geht1_integrals(int i,int j,vector<AO>& basis_ao, double eri_aa, double G1, double F2, double& ii_jj,double& ij_ij){

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


/*
void Hamiltonian_Fock_geht1(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
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

//          get_geht1_integrals(I,K,basis_ao,eri_aa,G1,F2,ii_kk,ik_ik);
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
//          get_geht1_integrals(I,J,basis_ao,eri_aa,G1,F2,ii_jj,ij_ij);
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

*/



void Hamiltonian_core_geht1(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
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
    Hao->M[i*Norb+i] = modprms.orb_params[i].IP; //modprms.PT[basis_ao[I].element].IP[basis_ao[I].ao_shell];
  }// for i


  //========================= Off-diagonal elements ==========================    
  double K_const = 0.0;// = 1.75;

  // Off-diagonal elements
  for(i=0;i<Norb;i++){
    I = basis_fo[i];
    for(j=i+1;j<Norb;j++){
        J = basis_fo[j];

        K_const = modprms.meht_k.get_K_value(I,J);

//        cout<<"I= "<<I<<" J= "<<J<<" K= "<<K_const<<endl;


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


        }// eht_formula==2 (Calzaferi)


        else if(prms.eht_formula==3){       // This is slightly generalized eht_formula == 0

          double K1_const = modprms.meht_k.get_K1_value(I,J); //!!!! Careful when using with self-consistent electrostatics !!!!
          // K1 value will also be used in different formula

//          double K2_const = modprms.meht_k.get_K2_value(I,J);


          Hao->M[i*Norb+j]  = (0.5*K_const*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j]) )*Sao->M[i*Norb+j] + K1_const;        


          Hao->M[j*Norb+i] = Hao->M[i*Norb+j];

        }// == 3
     
    }// for bj
  }// for i



}


// Here we actually start with energy (taken from charge equilibration) and then derive the correct Fock matrix
// no additional Fock matrix corrections are needed
// 
void Hamiltonian_Fock_geht1_v1(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                               vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                               Electronic* el,Electronic* el0, Memory* mem){


  int Norb = basis_fo.size(); // how many AOs are included in this fragment
  if(Norb!=el->Hao->num_of_cols){  
    cout<<"In Hamiltonian_Fock_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }


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


  // These matrices will be just zero
  *el->dFao_alp_dP_alp = 0.0;
  *el->dFao_alp_dP_bet = 0.0;
  *el->dFao_bet_dP_alp = 0.0;
  *el->dFao_bet_dP_bet = 0.0;



  //========================= Charge-corrected diagonal elements ==========================  
  if(prms.eht_electrostatics>=1){

  // Now use diagonal correction factors to update all elements of Fock matrix
  // SC-EHT
  for(int i=0;i<Norb;i++){
    int I = basis_fo[i];          // global index of orbital i
    int a = basis_ao[I].at_indx;  // global index of the atom on which orbital i is centered



    for(int j=0;j<Norb;j++){
      int J = basis_fo[j];          // global index of orbital j
      int b = basis_ao[J].at_indx;  // global index of the atom on which orbital j is centered



      //***************************************************
      if(a==b){  // orbitals are on the same atom - add fixed exchange to all orbitals

        // Exchange effects - only local, to preserve rotational invariance
        double K1_const = modprms.meht_k.get_K1_value(I,J); // a.u. of energy - exchange integral
//        double K4_const = modprms.meht_k.get_K4_value(I,J); // a.u. of length

//      double dist = (mol.R[a]-mol.R[b]).length();    
      double eri_ab = K1_const; // * ERFC(K4_const * dist);  


        if(prms.use_rosh){
          el->Fao_alp->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*eri_ab; 
          el->Fao_bet->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*eri_ab; 
        }
        else{
          el->Fao_alp->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*eri_ab; 
          el->Fao_bet->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*eri_ab; 
        }

      }// a==b
      //***************************************************      


 


      for(int A=0;A<mol.Nnucl;A++){  // over all atoms - effects of electronegativities
        double dQA_dPij_alp = 0.0;
        double dQA_dPij_bet = 0.0;
        
        double ca = (a==A)?0.5:0.0;
        double cb = (b==A)?0.5:0.0;
        
        dQA_dPij_alp = dQA_dPij_bet = -(ca+cb)*el->Sao->M[i*el->Norb+j];



        // Now add this to Fock matrix:
        int orb_A = at_orbitals[A][0];  // first orbital of the atom A - assume that coefficients do not depend on orbital, only
                                        // on atom type
        double xi_A = modprms.orb_params[orb_A].J_param1; //  PT[basis_ao[orb_A].element].J_param1[basis_ao[orb_A].ao_shell]; 
        double J_AA = modprms.orb_params[orb_A].J_param2; //PT[basis_ao[orb_A].element].J_param2[basis_ao[orb_A].ao_shell]; 
        double QA = mol.Mull_charges_gross[A];


        //************ This crrespond to contribution to energy:  xi_A * Q_A ***************
        el->Fao_alp->M[i*Norb+j] += xi_A * dQA_dPij_alp;
        el->Fao_bet->M[i*Norb+j] += xi_A * dQA_dPij_bet;
        //**********************************************************************************


        //************ This correspond to contribution to energy:  1/2*J_AA * Q_A^2 ********
        el->Fao_alp->M[i*Norb+j] += J_AA * QA * dQA_dPij_alp;
        el->Fao_bet->M[i*Norb+j] += J_AA * QA * dQA_dPij_bet;
        //**********************************************************************************

 
        if(prms.eht_electrostatics>=2){

          for(int B=0;B<mol.Nnucl;B++){  // over all atoms - effects of Coulomb terms
            double dQB_dPij_alp = 0.0;
            double dQB_dPij_bet = 0.0;
          
            ca = (a==B)?0.5:0.0;
            cb = (b==B)?0.5:0.0;
          
            dQB_dPij_alp = dQB_dPij_bet = -(ca+cb)*el->Sao->M[i*el->Norb+j];
          
            // Now add this to Fock matrix:
            int orb_B = at_orbitals[B][0];  // first orbital of the atom B - assume that coefficients do not depend on orbital, only
                                            // on atom type

            double J_BB = modprms.orb_params[orb_B].J_param2; //PT[basis_ao[orb_B].element].J_param2[basis_ao[orb_B].ao_shell]; 
            double QB = mol.Mull_charges_gross[B];
          
          
            // Nishimoto-Mataga-Weiss formula:
            double dist = (mol.R[A]-mol.R[B]).length();
            double J_AB = 1.2/(dist + 2.4/(J_AA+J_BB));
          
          
            // These are additional degrees of freedom to use in more elaborate 
            // functional J_AB
            //double K1_const = modprms.meht_k.get_K1_value(orb_A,orb_B); // a.u. of energy
            //double K2_const = modprms.meht_k.get_K2_value(orb_A,orb_B); // a.u. of energy
            //double K3_const = modprms.meht_k.get_K3_value(orb_A,orb_B); // a.u. of energy
            //double corr = ( K1_const + K2_const * dist + K3_const * dist*dist) * el->Sao->M[i*el->Norb+j];

            double K2_const = modprms.meht_k.get_K2_value(orb_A,orb_B); // a.u. of energy
            double K3_const = modprms.meht_k.get_K3_value(orb_A,orb_B); // a.u. of energy

            double K4_const = modprms.meht_k.get_K4_value(orb_A,orb_B); // a.u. of energy

            if(K4_const<0.0){  K4_const = 0.0; }

            double f = ERF(K4_const * dist);
            J_AB = J_AB * f +  (1.0 - f)*(K2_const + K3_const*dist*dist);  // add long-range electrostatics smoothly

//            cout<<"K3_const = "<<K3_const<<" energy contrib= "<<(1.0 -f)*K3_const*dist*dist<<endl;

          
            //***************** This corresponds to contribution to energy:  J_AB * Q_A * Q_B ****************
            el->Fao_alp->M[i*Norb+j] += 0.5*J_AB * (dQA_dPij_alp * QB + QA * dQB_dPij_alp);
            el->Fao_bet->M[i*Norb+j] += 0.5*J_AB * (dQA_dPij_bet * QB + QA * dQB_dPij_bet);
            //************************************************************************************************
          
          
          
          }// for B

        }// electrostatics >=2
     

      }// for A


    }// for j
  }// i

  }// if eht_electrostatics == 1
    
}





// Here we include charge effects at the level of Fock matrix, and intorduce correction (to get energy right)
void Hamiltonian_Fock_geht1(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                            vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                            Electronic* el,Electronic* el0, Memory* mem){


  int Norb = basis_fo.size(); // how many AOs are included in this fragment
  if(Norb!=el->Hao->num_of_cols){  
    cout<<"In Hamiltonian_Fock_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }

  // Compute mapping from globar orbital index K to local orbital index k (local for this fragment)
  vector<int> loc_indx(el->Norb, -1);  // dimension of this array is the full dimension 
  
  for(int n=0;n<Norb;n++){ // loop over fragment orbitals
    loc_indx[ basis_fo[n] ] = n;  
  }



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


  *el->dFao_alp_dP_alp = 0.0;
  *el->dFao_alp_dP_bet = 0.0;
  *el->dFao_bet_dP_alp = 0.0;
  *el->dFao_bet_dP_bet = 0.0;


  //========================= Charge-corrected diagonal elements ==========================  
  // Diagonal electrostatics
  vector<double> Wii(Norb,0.0);
  vector<double> Aii(Norb,0.0);
  vector<int> aindx(Norb,0.0);

  // Depending of SCE, do different type of correction to diagonal elements
  // Precompute electrostatic correction factor for all orbitals
  for(int i=0;i<Norb;i++){

    int I = basis_fo[i];
    int A = basis_ao[I].at_indx;  // basis_fo[i] on atom a - global atom index

    double A1i = modprms.PT[basis_ao[I].element].J_param1[basis_ao[I].ao_shell]; 
    double A2i = modprms.PT[basis_ao[I].element].J_param2[basis_ao[I].ao_shell]; 
    double Qi = mol.Mull_charges_gross[basis_ao[I].at_indx];


    // On-size terms
    // First term always pushes energy up for negative charge
    Aii[i] = A1i + 2.0*A2i*Qi;
    aindx[i] = A; 

    Wii[i] = -A1i*Qi + A2i*Qi*Qi;



    // This term describes self-energy of the bond chanrge density! 
    // It does not depend on density matrix, but is similar to atomic on-site corrections 
    for(int a=0;a<mol.Nnucl;a++){
      if(a!=basis_ao[I].at_indx){  // orbitals on all other atoms
        for(int kk=0;kk<at_orbitals[a].size();kk++){  // over all orbitals on atom a


          int J = at_orbitals[a][kk];  // kk-th orbital of set on atom a, J - is the global index of the orbital
          int j = loc_indx[J];         // local orbital index

          // To preserve rotational invariance, the constants K1 and K2 should depend only on atom type
          double K1_const = modprms.meht_k.get_K1_value(I,J); // a.u. of energy
          double K2_const = modprms.meht_k.get_K2_value(I,J); // a.u. of energy
          double K3_const = modprms.meht_k.get_K3_value(I,J); // a.u. of energy

          double dist = (mol.R[basis_ao[I].at_indx]-mol.R[a]).length();

          // dist = 0.0 => corr = 0,     so no change for short-range
          // dist = inf => corr = K1 * S so there is an additional splitting for long range (yet, decays as overlap)
          double corr = ( K1_const + K2_const * dist + K3_const * dist*dist) * el->Sao->M[i*el->Norb+j];

          el->dFao_alp_dP_alp->M[i*Norb+j] += corr;
          el->dFao_bet_dP_bet->M[i*Norb+j] += corr;


          el->Fao_alp->M[i*Norb+j] += corr * el->P_alp->M[i*Norb+j];
          el->Fao_bet->M[i*Norb+j] += corr * el->P_bet->M[i*Norb+j];
            

        }// for k
      }// if - different atoms



/*
      // Contribution from periodic images - well, this is not really converging, but at least first environment will be
      // taken into account
      if(prms.x_period>=1 || prms.y_period>=1 || prms.z_period>=1){

        for(int nx=-prms.x_period;nx<=prms.x_period;nx++){
          for(int ny=-prms.y_period;ny<=prms.y_period;ny++){
            for(int nz=-prms.z_period;nz<=prms.z_period;nz++){
        
              if(nx==0 && ny==0 && nz==0){ }  // already included
              else{
                
                VECTOR TV = nx*prms.t1 + ny*prms.t2 + nz*prms.t3;

                double K3_const = modprms.meht_k.get_K3_value(I,J);
                double K4_const = modprms.meht_k.get_K4_value(I,J);     
                double dist2 = (mol.R[basis_ao[I].at_indx]-mol.R[basis_ao[J].at_indx]-TV).length();
                // This summation corresponds to k = 0 (Gamma-point)
     
                // This approximates - Q[n] * eri[]
                Wii[i] -= K3_const*K3_const*mol.Mull_charges_gross[basis_ao[J].at_indx]/sqrt((K4_const*K4_const) + dist2);
     
              }
        
            }// for nz
          }// for ny
        }// for nx 

      }// if - periodic


*/
    }// for j

  }// for i



  // Now use diagonal correction factors to update all elements of Fock matrix
  // SC-EHT
  for(int i=0;i<Norb;i++){
    for(int j=0;j<Norb;j++){

      // This Fock matrix is used to compute energies
      el->Fao_alp->M[i*Norb+j] += 0.5*(Wii[i]+Wii[j])*el->Sao->M[i*el->Norb+j];
      el->Fao_bet->M[i*Norb+j] += 0.5*(Wii[i]+Wii[j])*el->Sao->M[i*el->Norb+j];



      // There should be a correction to Fock matrix used in SCF calculations:
      // compute this part:  dF/dP:

      int I = basis_fo[i];          // global index of orbital i
      int a = basis_ao[I].at_indx;  // global index of the atom on which orbital i is centered

      int J = basis_fo[j];          // global index of orbital j
      int b = basis_ao[J].at_indx;  // global index of the atom on which orbital j is centered


      double dQi_dPij_alp = 0.0;
      double dQi_dPij_bet = 0.0;
      if(a==aindx[i] || b==aindx[i]){ dQi_dPij_alp = dQi_dPij_bet = -0.5*el->Sao->M[i*el->Norb+j]; }

      double dQj_dPij_alp = 0.0;
      double dQj_dPij_bet = 0.0;
      if(a==aindx[j] || b==aindx[j]){ dQj_dPij_alp = dQj_dPij_bet  = -0.5*el->Sao->M[i*el->Norb+j]; }


      el->dFao_alp_dP_alp->M[i*Norb+j] -= 0.5*(Aii[i]*dQi_dPij_alp + Aii[j]*dQj_dPij_alp)*el->Sao->M[i*el->Norb+j];
      el->dFao_alp_dP_bet->M[i*Norb+j] -= 0.5*(Aii[i]*dQi_dPij_bet + Aii[j]*dQj_dPij_bet)*el->Sao->M[i*el->Norb+j];
      el->dFao_bet_dP_alp->M[i*Norb+j] -= 0.5*(Aii[i]*dQi_dPij_alp + Aii[j]*dQj_dPij_alp)*el->Sao->M[i*el->Norb+j];
      el->dFao_bet_dP_bet->M[i*Norb+j] -= 0.5*(Aii[i]*dQi_dPij_bet + Aii[j]*dQj_dPij_bet)*el->Sao->M[i*el->Norb+j];


      // Now, yet another correction - the off-diagonal elements shouldn't go to zero on infinity - otherwise
     //
//      double K1_const = modprms.meht_k.get_K1_value(I,J);
//      double K2_const = modprms.meht_k.get_K2_value(I,J);


    
    }// for j
  }// i



/*
  // This is generalization of INDO Fock
  // Formation of the Fock matrix: add Coulomb and Exchange parts    
  for(int i=0;i<Norb;i++){
    int I = basis_fo[i];
    int A = basis_ao[I].at_indx;  // basis_fo[i] on atom a - global atom index


    for(int j=0;j<Norb;j++){
      int J = basis_fo[j];
      int B = basis_ao[J].at_indx;

 
      // On-site Coulomb and exchange
      if(I==J){  // Diagonal terms  - same orbital

        for(int kk=0;kk<at_orbitals[A].size();kk++){       // for all orbitals on atom A
          int K = at_orbitals[A][kk];                      // global orbital index
          int k = loc_indx[K];                             // local orbital index
       

          double ii_kk, ik_ik; ii_kk = ik_ik = 0.0;
          ii_kk = modprms.meht_k.get_K1_value(I,K);

          if(I==K){  ik_ik = ii_kk; }
          else{
                     ik_ik = modprms.meht_k.get_K2_value(I,K);
          }



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


          ii_jj = modprms.meht_k.get_K1_value(I,J);
          if(I==J){ ij_ij = ii_jj; } // this will never happen, because we have I != J  by definition
          else{
                    ij_ij = modprms.meht_k.get_K2_value(I,J); // for H-only compounds this is zero
          }


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

          double dist = (mol.R[A]-mol.R[B]).length();

          double A2i = modprms.PT[basis_ao[I].element].J_param2[basis_ao[I].ao_shell]; 
          double A2j = modprms.PT[basis_ao[J].element].J_param2[basis_ao[J].ao_shell]; 
          double a2ij = 0.5*(A2i + A2j)/eV;  // unitless

          double eri_ab = modprms.meht_k.get_K1_value(I,J) * ERFC(a2ij*a2ij * dist);  




          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*eri_ab; 
            el->Fao_bet->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*eri_ab; 
          }
          else{
            el->Fao_alp->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*eri_ab; 
            el->Fao_bet->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*eri_ab; 
          }


        }// different orbitals on different atoms
    
      }// I!=J


  
  
    }// for j
  }// for i
*/
//  cout<<"end of Hamiltonian_Fock_geht\n";

    
}




