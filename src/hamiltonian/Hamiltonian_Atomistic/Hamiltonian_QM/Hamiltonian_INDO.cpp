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

#include "Hamiltonian_INDO.h"

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




namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_qm{


vector<int> compute_sorb_indices
( int sz, vector<AO>& basis_ao, vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
){

  // Compute global indices of s-type orbitals on each atom
  vector<int> sorb_indx(sz,0); // global index of s-type orbital on i-th atom

  for(int a=0;a<sz;a++){  // for all atoms
    for(int i=0;i<atom_to_ao_map[a].size();i++){  // all orbitals on given atom (i-dummy)

      int I = atom_to_ao_map[a][i];  // i-th AO on atom a, I - is the global index of this AO in the given basis

      if(basis_ao[I].ao_shell_type=="s" ){  sorb_indx[a] = I; }
          
    }// for j
  }// for i

  return sorb_indx;
}


void compute_indo_core_parameters
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  vector<int>& sorb_indx,
  int opt, int a, int b, double& eri, double& V_AB){

  // Compute ERIs and V_AB

  int I = sorb_indx[a];
  int J = sorb_indx[b];
     
  // ERI
  eri = electron_repulsion_integral(basis_ao[I],basis_ao[I],basis_ao[J],basis_ao[J]); // eri[a][b]

  // V_AB
  int B = b; //ao_to_atom_map[b]; // global index of atom on orbital b
  double Zeff = modprms.PT[syst.Atoms[B].Atom_element].Zeff; 

  if(opt==0){
    V_AB = Zeff*nuclear_attraction_integral(basis_ao[J],basis_ao[J], syst.Atoms[B].Atom_RB.rb_cm);// V_AB[a][b]
  }
  else if(opt==1){
    V_AB = Zeff*eri; //[a*sz+b];
  }

}

void compute_indo_core_parameters_derivs
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  vector<int>& sorb_indx,
  int opt, int a, int b, int c, VECTOR& deri, VECTOR& dV_AB){

  // Compute dERIs and dV_AB:  d ERI[a][b] / dR[c]  and d V[a][b] / dR[c]

  int I = sorb_indx[a];
  int J = sorb_indx[b];
  int K = sorb_indx[c]; 
     
  // ERI
  //eri = electron_repulsion_integral(basis_ao[I],basis_ao[I],basis_ao[J],basis_ao[J]); // eri[a][b]


  VECTOR DA,DB,DC,DD;
  double eri = electron_repulsion_integral(&basis_ao[I],&basis_ao[I],&basis_ao[J],&basis_ao[J],1,1,DA,DB,DC,DD);

  VECTOR deri_dc; deri_dc = 0.0;
  if(a==b){ ;; }
  else{
    if(c==a){ deri += (DA + DB); }
    if(c==b){ deri += (DC + DD); }
  }


  // V_AB
  double V_AB = 0.0;
  int B = b; //ao_to_atom_map[b]; // global index of atom on orbital b
  double Zeff = modprms.PT[syst.Atoms[B].Atom_element].Zeff; 

  if(opt==0){
    V_AB = Zeff*nuclear_attraction_integral(basis_ao[J],basis_ao[J], syst.Atoms[B].Atom_RB.rb_cm);// V_AB[a][b]

    double nai = nuclear_attraction_integral(basis_ao[J],basis_ao[J], syst.Atoms[B].Atom_RB.rb_cm, 1, 1, DA, DB, DC);// V_AB[a][b]

    // Alright - i'm not really sure about this, so skip this part for now
    //if(K==J){ dV_AB += (DA + DB); }
    //if(K==J){ dV_AB += (DA + DB); }

  }
  else if(opt==1){
    V_AB = Zeff*eri; //[a*sz+b];
    dV_AB = Zeff * deri;
  }

}





void indo_core_parameters
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int opt, int DF){
// opt == 0 - cndo
// opt == 1 - indo
  modprms.indo_opt = opt;

//  int DF = 0;
  int i,j,a,b,i1,a1,I,J,A;
  VECTOR da,db,dc;

  if(DF){ cout<<"in indo_core_parameters\n"; }

  // Allocate memory if needed
  int sz = syst.Number_of_atoms; // number of atoms in the system

  if(modprms.eri.size()!=sz*sz){  cout<<"In indo_core_parameters: eri array is not allocated\nDo allocation...\n"; 
    modprms.eri.clear(); modprms.eri = vector<double>(sz*sz,0.0); 
  }
  if(modprms.V_AB.size()!=sz*sz){  cout<<"In indo_core_parameters: V_AB array is not allocated\nDo allocation...\n"; 
    modprms.V_AB.clear(); modprms.V_AB = vector<double>(sz*sz,0.0); 
  }

  vector<int> sorb_indx;
  sorb_indx = compute_sorb_indices(sz,basis_ao,atom_to_ao_map,ao_to_atom_map);


  // Printing mapping
  if(DF){ 
    cout<<"i - runs over indices of atoms in given system\n";
    cout<<"sorb_indx[i] - is the global index of s-type orbital (assuming only one) centered on atom i\n";
    for(i=0;i<sorb_indx.size();i++){
      cout<<"i= "<<i<<" sorb_indx[i]= "<<sorb_indx[i]<<endl;
    }
  }
    
/*
  // Compute ERIs and V_AB
  for(a=0;a<sz;a++){
    for(b=0;b<sz;b++){


      I = sorb_indx[a];
      J = sorb_indx[b];
     
      // ERI
      modprms.eri[a*sz+b] = electron_repulsion_integral(basis_ao[I],basis_ao[I],basis_ao[J],basis_ao[J]);

      // V_AB
      int B = b; //ao_to_atom_map[b]; // global index of atom on orbital b
      if(opt==0){
        modprms.V_AB[a*sz+b] = modprms.PT[syst.Atoms[B].Atom_element].Zeff*nuclear_attraction_integral(basis_ao[J],basis_ao[J], syst.Atoms[B].Atom_RB.rb_cm);// V_AB[a][b]
      }
      else if(opt==1){
        modprms.V_AB[a*sz+b] = modprms.PT[syst.Atoms[B].Atom_element].Zeff*modprms.eri[a*sz+b];
      }

      if(DF){ cout<<"a= "<<a<<" b= "<<b<<" I= "<<I<<" J= "<<J<<" eri= "<<modprms.eri[a*sz+b]<<" V_AB= "<<modprms.V_AB[a*sz+b]<<endl; }

    }// for j
  }// for i
*/

  // Compute ERIs and V_AB
  for(a=0;a<sz;a++){
    for(b=0;b<sz;b++){

    compute_indo_core_parameters
    (syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, sorb_indx, opt, a, b, modprms.eri[a*sz+b], modprms.V_AB[a*sz+b] ); 

    }// for b
  }// for a

}




void Hamiltonian_core_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF
){

  int i,j,k,a,b,I,J,A,B;
  VECTOR dIdA,dIdB;
  //cout<<"in Hamiltonian_core_indo\n";

  int Norb = basis_ao.size(); // how many AOs are included in this fragment
  if(Norb!=Hao->num_of_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_indo is called\n";
    cout<<"In Hamiltonian_core_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  int sz = syst.Number_of_atoms; // number of atoms in this fragment

  if(modprms.eri.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_indo: size of auxiliary eri array is not right\n"; exit(0);}
  if(modprms.V_AB.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_indo: size of auxiliary V_AB array is not right\n"; exit(0);}



  //----------- Compute core terms of the Hamiltonian ------------
  *Hao = 0.0;

  for(i=0;i<Norb;i++){  // global orbital indices
    // values of IP are different for cndo (just IPs) and cndo2 and indo ( 0.5*(IP + EA) )
    a = ao_to_atom_map[i];
    Hao->M[i*Norb+i] += modprms.PT[basis_ao[i].element].IP[basis_ao[i].ao_shell];

    if(DF){
      cout<<"Setting diagonal element i = "<<i<<endl;
      cout<<"Contribution from IP = "<<Hao->M[i*Norb+i]<<endl;
    }

    /// The code below is the same for CNDO, CNDO2 and INDO - but the difference comes in use of different G1 and F2 parameters
    /// for CNDO and CNDO2 they are zero
    double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
    double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
    
    if(DF){
      cout<<"a= "<<a<<" G1= "<<G1<<" F2= "<<F2<<" Atom[a].Atom_Z= "<<syst.Atoms[a].Atom_Z
          <<" basis_ao[i].ao_shell_type= "<<basis_ao[i].ao_shell_type<<endl;
    }
   
    /// Eqs. 3.17 - 3.23 from Pople, Beveridge, Dobosh, JCP 47, 2026 (1967)
    ///
    int Z = syst.Atoms[a].Atom_Z;  // modprms.PT[elt].Z - e.g. 6 for C

         if(Z==1){ Hao->M[i*Norb+i] -= 0.5*modprms.eri[a*sz+a]; }  // H
    else if(Z==3 || Z==11){ 

      if(basis_ao[i].ao_shell_type=="s"){   Hao->M[i*Norb+i] -= 0.5*modprms.eri[a*sz+a]; } // s
      else if(basis_ao[i].ao_shell_type=="p"){  Hao->M[i*Norb+i] -= (0.5*modprms.eri[a*sz+a] - G1/12.0); } // p

    }  // Li or Na

    else if(Z==4 || Z==12){ 

      if(basis_ao[i].ao_shell_type=="s"){  Hao->M[i*Norb+i] -= (1.5*modprms.eri[a*sz+a] - 0.5*G1); } // s
      else if(basis_ao[i].ao_shell_type=="p"){  Hao->M[i*Norb+i] -= (1.5*modprms.eri[a*sz+a] - 0.25*G1); } // p

    }  // Be or Mg

    else if( (Z>=5 && Z<=9) || (Z>=13 && Z<=18)){ // B - F or Al - Cl

      double Z_core = modprms.PT[basis_ao[i].element].Nval; // core charge of atom 

      if(basis_ao[i].ao_shell_type=="s"){      
        Hao->M[i*Norb+i] -= ( ( Z_core - 0.5 )*modprms.eri[a*sz+a] - (1.0/6.0)*( Z_core - 1.5 )*G1); // s
      }
      else if(basis_ao[i].ao_shell_type=="p"){
        Hao->M[i*Norb+i] -= ( ( Z_core - 0.5 )*modprms.eri[a*sz+a] - (1.0/3.0)*G1 - 0.08*( Z_core - 2.5 )*F2); // p
      }

    }
    else{  cout<<"Error: INDO is not implemented for elements beyond Cl\n"; exit(0);  }
    if(DF){ cout<<" + Contribution from Frank-Condon factors = "<<Hao->M[i*Norb+i]<<endl;   }

    
    //----------------- Coulombic terms --------------

    for(b=0;b<sz;b++){
      if(b!=a){
        Hao->M[i*Norb+i] -= modprms.V_AB[a*sz+b];  //  = V_AB = Z_B * eri[A][B] -in INDO
      }
    }// for b
    if(DF){  cout<<" + Contribution from Coulombic terms = "<<Hao->M[i*Norb+i]<<endl;   }


    //-------------- Off-diagonal terms of the core matrix ---------
    for(j=0;j<Norb;j++){

      if(j!=i){
        b = ao_to_atom_map[j];

        if(b==a){ ;; }  // different orbitals centered on the same atom - give zero (not true for hybrid orbitals)
        else{           // centered on different atoms - use overlap formula

          // Overlap is set to identity in INDO, so need to recompute it explicitly
          double sao_ij = gaussian_overlap(basis_ao[i],basis_ao[j]); // 0, dIdA,dIdB), mem->aux, mem->n_aux);

          Hao->M[i*Norb+j] += 0.5*(modprms.PT[basis_ao[i].element].beta0[basis_ao[i].ao_shell] 
                                 + modprms.PT[basis_ao[j].element].beta0[basis_ao[j].ao_shell]) * sao_ij;
        }
      }// j!=i

    }// for j    

  }// for i

}


void Hamiltonian_core_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF
){

  Hamiltonian_core_indo( syst, basis_ao, prms, modprms,  atom_to_ao_map, ao_to_atom_map, &Hao, &Sao, DF);

}


void Hamiltonian_core_deriv_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF,
  int c,
  MATRIX* dHao_dx, MATRIX* dHao_dy, MATRIX* dHao_dz, 
  MATRIX* dSao_dx, MATRIX* dSao_dy, MATRIX* dSao_dz
){

  //================ Basically, here we compute derivatives of the core Hamiltonian ========================

  int i,j,k,a,b,I,J,A,B;
  VECTOR dIdA,dIdB;
//  cout<<"in Hamiltonian_core_deriv_indo\n";

  int Norb = basis_ao.size(); // how many AOs are included in this fragment
  if(Norb!=Hao->num_of_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_indo is called\n";
    cout<<"In Hamiltonian_core_deriv_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  if(Norb!=dHao_dx->num_of_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_indo is called\n";
    cout<<"In Hamiltonian_core_deriv_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  if(Norb!=dHao_dy->num_of_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_indo is called\n";
    cout<<"In Hamiltonian_core_deriv_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  if(Norb!=dHao_dz->num_of_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_indo is called\n";
    cout<<"In Hamiltonian_core_deriv_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }

  int sz = syst.Number_of_atoms; // number of atoms in this fragment

//  cout<<"sz= "<<sz<<endl;

  if(modprms.eri.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_deriv_indo: size of auxiliary eri array is not right\n"; exit(0);}
  if(modprms.V_AB.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_deriv_indo: size of auxiliary V_AB array is not right\n"; exit(0);}


  //----------- Compute core terms of the Hamiltonian ------------
  *Hao = 0.0;
  *dHao_dx = 0.0;
  *dHao_dy = 0.0;
  *dHao_dz = 0.0;

  *dSao_dx = 0.0;
  *dSao_dy = 0.0;
  *dSao_dz = 0.0;
  


  vector<int> sorb_indx;
  sorb_indx = compute_sorb_indices(sz, basis_ao, atom_to_ao_map, ao_to_atom_map);


  for(i=0;i<Norb;i++){  // global orbital indices
    // values of IP are different for cndo (just IPs) and cndo2 and indo ( 0.5*(IP + EA) )
    a = ao_to_atom_map[i];    
    Hao->M[i*Norb+i] += modprms.PT[basis_ao[i].element].IP[basis_ao[i].ao_shell];

    if(DF){
      cout<<"Setting diagonal element i = "<<i<<endl;
      cout<<"Contribution from IP = "<<Hao->M[i*Norb+i]<<endl;
    }

    /// The code below is the same for CNDO, CNDO2 and INDO - but the difference comes in use of different G1 and F2 parameters
    /// for CNDO and CNDO2 they are zero
    double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
    double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
    
    if(DF){
      cout<<"a= "<<a<<" G1= "<<G1<<" F2= "<<F2<<" Atom[a].Atom_Z= "<<syst.Atoms[a].Atom_Z
          <<" basis_ao[i].ao_shell_type= "<<basis_ao[i].ao_shell_type<<endl;
    }
   
    /// Eqs. 3.17 - 3.23 from Pople, Beveridge, Dobosh, JCP 47, 2026 (1967)
    ///
    int Z = syst.Atoms[a].Atom_Z;  // modprms.PT[elt].Z - e.g. 6 for C


    //================== The portion below does not contribute to the derivatives ====================

         if(Z==1){      Hao->M[i*Norb+i] -= 0.5*modprms.eri[a*sz+a];    }  // H
    else if(Z==3 || Z==11){ 

      if(basis_ao[i].ao_shell_type=="s"){   Hao->M[i*Norb+i] -= 0.5*modprms.eri[a*sz+a]; } // s
      else if(basis_ao[i].ao_shell_type=="p"){  Hao->M[i*Norb+i] -= (0.5*modprms.eri[a*sz+a] - G1/12.0); } // p

    }  // Li or Na

    else if(Z==4 || Z==12){ 

      if(basis_ao[i].ao_shell_type=="s"){  Hao->M[i*Norb+i] -= (1.5*modprms.eri[a*sz+a] - 0.5*G1); } // s
      else if(basis_ao[i].ao_shell_type=="p"){  Hao->M[i*Norb+i] -= (1.5*modprms.eri[a*sz+a] - 0.25*G1); } // p

    }  // Be or Mg

    else if( (Z>=5 && Z<=9) || (Z>=13 && Z<=18)){ // B - F or Al - Cl

      double Z_core = modprms.PT[basis_ao[i].element].Nval; // core charge of atom 

      if(basis_ao[i].ao_shell_type=="s"){      
        Hao->M[i*Norb+i] -= ( ( Z_core - 0.5 )*modprms.eri[a*sz+a] - (1.0/6.0)*( Z_core - 1.5 )*G1); // s
      }
      else if(basis_ao[i].ao_shell_type=="p"){
        Hao->M[i*Norb+i] -= ( ( Z_core - 0.5 )*modprms.eri[a*sz+a] - (1.0/3.0)*G1 - 0.08*( Z_core - 2.5 )*F2); // p
      }

    }
    else{  cout<<"Error: INDO is not implemented for elements beyond Cl\n"; exit(0);  }
    if(DF){ cout<<" + Contribution from Frank-Condon factors = "<<Hao->M[i*Norb+i]<<endl;   }

    
    //----------------- Coulombic terms: Contribution to diagonal elements --------------


    for(b=0;b<sz;b++){
//      cout<<"i= "<<i<<" a= "<<a<<" b= "<<b<<endl;

      if(b!=a){
        Hao->M[i*Norb+i] -= modprms.V_AB[a*sz+b];  //  = V_AB = Z_B * eri[A][B] -in INDO

        VECTOR deri, dV_AB;    
        compute_indo_core_parameters_derivs(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map,
                                            sorb_indx, modprms.indo_opt, a, b, c, deri, dV_AB);

//        cout<<"i= "<<i<<" a= "<<a<<" b= "<<b<<" c= "<<c<<" deri= "<<deri<<" dV_AB= "<<dV_AB<<endl;

        dHao_dx->M[i*Norb+i] -= dV_AB.x;
        dHao_dy->M[i*Norb+i] -= dV_AB.y;
        dHao_dz->M[i*Norb+i] -= dV_AB.z;


      }
    }// for b
    if(DF){  cout<<" + Contribution from Coulombic terms = "<<Hao->M[i*Norb+i]<<endl;   }


    //-------------- Off-diagonal terms of the core matrix ---------
    for(j=0;j<Norb;j++){

      if(j!=i){
        b = ao_to_atom_map[j];

        if(b==a){ ;; }  // different orbitals centered on the same atom - give zero (not true for hybrid orbitals)
        else{           // centered on different atoms - use overlap formula

          // Overlap is set to identity in INDO, so need to recompute it explicitly
          //double sao_ij = gaussian_overlap(basis_ao[i],basis_ao[j]); // 0, dIdA,dIdB), mem->aux, mem->n_aux);

          VECTOR dSda, dSdb, dSdc;
          double sao_ij = gaussian_overlap(&basis_ao[i],&basis_ao[j], 1, 1, dSda, dSdb);
          // i - on atom a
          // j - on atom b

          dSdc = 0.0;
          if(c==a){ dSdc += dSda; }
          if(c==b){ dSdc += dSdb; }

          double beta_ij = 0.5*(modprms.PT[basis_ao[i].element].beta0[basis_ao[i].ao_shell] + modprms.PT[basis_ao[j].element].beta0[basis_ao[j].ao_shell]);

          Hao->M[i*Norb+j] += beta_ij * sao_ij;

          dHao_dx->M[i*Norb+j] += beta_ij * dSdc.x;
          dHao_dy->M[i*Norb+j] += beta_ij * dSdc.y;
          dHao_dz->M[i*Norb+j] += beta_ij * dSdc.z;

          dSao_dx->M[i*Norb+j] += dSdc.x;
          dSao_dy->M[i*Norb+j] += dSdc.y;
          dSao_dz->M[i*Norb+j] += dSdc.z;


        }
      }// j!=i

    }// for j    

  }// for i

}

void Hamiltonian_core_deriv_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF,
  int c,
  MATRIX& dHao_dx, MATRIX& dHao_dy, MATRIX& dHao_dz, 
  MATRIX& dSao_dx, MATRIX& dSao_dy, MATRIX& dSao_dz
){

  Hamiltonian_core_deriv_indo
  ( syst, basis_ao, prms, modprms,  atom_to_ao_map, ao_to_atom_map,  &Hao, &Sao, DF, c,
    &dHao_dx, &dHao_dy, &dHao_dz,   &dSao_dx, &dSao_dy, &dSao_dz);

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
        (basis_ao[i].z_exp == basis_ao[j].z_exp)){
      ij_ij = ii_jj = eri_aa + 0.16*F2;      // xx_xx = yy_yy = zz_zz
    }
    else{
      ij_ij = 0.12*F2;                         // xy_xy = xz_xz = ...
      ii_jj = eri_aa - 0.08*F2;                // xx_yy = xx_zz = ...
    }
  }

  //=====================================================================================================

}


void Hamiltonian_Fock_indo(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
                           Control_Parameters& prms, Model_Parameters& modprms,
                           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                          ){

  
// This function constructs INDO (CNDO/2) Fock matrix
// Unrestricted formulation

  int i,j,k,n,I,J,K,a,b,A,B;

  int Norb = basis_ao.size(); // how many AOs are included in this fragment
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


  vector<double> Zeff(syst.Number_of_atoms, 0.0);
  vector<double> Mull_charges_gross(syst.Number_of_atoms, 0.0);
  vector<double> Mull_charges_net(syst.Number_of_atoms, 0.0);

  for(a=0;a<syst.Number_of_atoms;a++){ Zeff[a] = modprms.PT[syst.Atoms[a].Atom_element].Zeff; } // e.g. 4 for STO-3G C

  update_Mull_charges(ao_to_atom_map, Zeff, el->Mull_orb_pop_gross, el->Mull_orb_pop_net, Mull_charges_gross, Mull_charges_net);

  for(a=0;a<syst.Number_of_atoms;a++){ 
    syst.Atoms[a].Atom_mull_charge_gross = Mull_charges_gross[a]; 
    syst.Atoms[a].Atom_mull_charge_net = Mull_charges_net[a]; 
  }



    
  // Formation of the Fock matrix: add Coulomb and Exchange parts    
  for(i=0;i<Norb;i++){
    a = ao_to_atom_map[i];

    for(j=0;j<Norb;j++){
      b = ao_to_atom_map[j];

      if(i==j){  // Diagonal terms

//        cout<<"Diagonal terms, i= "<<i<<endl;
//        cout<<"el->Fao_alp->M[i*Norb+i] = "<<el->Fao_alp->M[i*Norb+i]<<endl;

        for(int kk=0;kk<atom_to_ao_map[a].size();kk++){    // for all orbitals on atom a
          k = atom_to_ao_map[a][kk];                       // global orbital index of AO kk on atom a


          double ii_kk, ik_ik; ii_kk = ik_ik = 0.0;
          double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
          double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
          double eri_aa = modprms.eri[a*syst.Number_of_atoms + a];

          get_integrals(i,k,basis_ao,eri_aa,G1,F2,ii_kk,ik_ik);

/*
          cout<<"kk= "<<kk<<" k= "<<k
              <<" basis_ao[i].ao_shell_type= "<<basis_ao[i].ao_shell_type
              <<" basis_ao[k].ao_shell_type= "<<basis_ao[k].ao_shell_type
              <<"basis_ao[i].exponents = ("<<basis_ao[i].x_exp<<","<<basis_ao[i].y_exp<<","<<basis_ao[i].z_exp
              <<"basis_ao[k].exponents = ("<<basis_ao[k].x_exp<<","<<basis_ao[k].y_exp<<","<<basis_ao[k].z_exp
              <<" el->P->M[k*Norb+k]= "<<el->P->M[k*Norb+k]
              <<" el->P_alp->M[k*Norb+k]= "<<el->P_alp->M[k*Norb+k]
              <<" ii_kk= "<<ii_kk<<" ik_ik= "<<ik_ik
              <<endl;
*/

          if(prms.use_rosh){ // Restricted open-shell
            el->Fao_alp->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - 0.5*el->P->M[k*Norb+k]*ik_ik);
            el->Fao_bet->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - 0.5*el->P->M[k*Norb+k]*ik_ik);

          }
          else{ // unrestricted
            el->Fao_alp->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - el->P_alp->M[k*Norb+k]*ik_ik);
            el->Fao_bet->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - el->P_bet->M[k*Norb+k]*ik_ik);

          }

//          cout<<"at this point el->Fao_bet->M[i*Norb+i]= "<<el->Fao_bet->M[i*Norb+i]<<endl;
  
        }// for kk - all orbitals on atom A
  

        // Contributions from all other atoms to the diagonal terms
        // don't worry that b determined above will be rewritten - this is ok
        for(b=0;b<syst.Number_of_atoms;b++){

          if(b!=a){

            // Compute density matrix due to all orbitals on atom b 
            el->Fao_alp->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*modprms.eri[a*syst.Number_of_atoms+b]; 
          }
        }// for b
  
  
      }// i==j
      else{      // Off-diagonal terms
  
        if(a==b){ // different orbitals are on the same atom
          double ij_ij,ii_jj; ij_ij = ii_jj = 0.0;

          double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
          double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
          double eri_aa = modprms.eri[a*syst.Number_of_atoms+a];
          get_integrals(i,j,basis_ao,eri_aa,G1,F2,ii_jj,ij_ij);

          //cout<<"i= "<<i<<" j= "<<j<<" ii_jj= "<<ii_jj<<", ij_ij= "<<ij_ij<<endl;

          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - 0.5*el->P->M[i*Norb+j])*ij_ij - 0.5*el->P->M[i*Norb+j]*ii_jj );
            el->Fao_bet->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - 0.5*el->P->M[i*Norb+j])*ij_ij - 0.5*el->P->M[i*Norb+j]*ii_jj );
          }
          else{  
            el->Fao_alp->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - el->P_alp->M[i*Norb+j])*ij_ij - el->P_alp->M[i*Norb+j]*ii_jj );
            el->Fao_bet->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - el->P_bet->M[i*Norb+j])*ij_ij - el->P_bet->M[i*Norb+j]*ii_jj );
          }
  
        }// a==b - different orbitals are on the same atom
  
        else{ // different orbitals are on different atoms

          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
          }
          else{
            el->Fao_alp->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
          }

        }
  
  
      }// i!=j
  
  
    }// for j
  }// for i

}


void Hamiltonian_Fock_indo(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
                           Control_Parameters& prms, Model_Parameters& modprms,
                           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                          ){

  Hamiltonian_Fock_indo(&el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

}




void Hamiltonian_Fock_derivs_indo
( Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int c, 
  MATRIX* dHao_dx, MATRIX* dHao_dy, MATRIX* dHao_dz,
  MATRIX* dFao_alp_dx, MATRIX* dFao_alp_dy, MATRIX* dFao_alp_dz,
  MATRIX* dFao_bet_dx, MATRIX* dFao_bet_dy, MATRIX* dFao_bet_dz
){
// This functions constructs INDO (CNDO/2) Fock matrix
// with gradients, also for unrestricted formulation

  int i,j,k,n,I,J,K,a,b,A,B;



  int Norb = basis_ao.size(); // how many AOs are included in this fragment
  if(Norb!=el->Hao->num_of_cols){  
    cout<<"In Hamiltonian_Fock_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }


  // Formation of the Fock matrix: Core part
  *el->Fao_alp = *el->Hao;
  *el->Fao_bet = *el->Hao;

  *dFao_alp_dx = *dHao_dx;
  *dFao_alp_dy = *dHao_dy;
  *dFao_alp_dz = *dHao_dz;

  *dFao_bet_dx = *dHao_dx;
  *dFao_bet_dy = *dHao_dy;
  *dFao_bet_dz = *dHao_dz;


  vector<int> sorb_indx;
  sorb_indx = compute_sorb_indices(syst.Number_of_atoms, basis_ao, atom_to_ao_map, ao_to_atom_map);


  // Update charges
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



    
  // Formation of the Fock matrix: add Coulomb and Exchange parts    
  for(i=0;i<Norb;i++){
    a = ao_to_atom_map[i];

    for(j=0;j<Norb;j++){
      b = ao_to_atom_map[j];

      if(i==j){  // Diagonal terms

        for(int kk=0;kk<atom_to_ao_map[a].size();kk++){    // for all orbitals on atom a
          k = atom_to_ao_map[a][kk];                       // global orbital index of AO kk on atom a


          double ii_kk, ik_ik; ii_kk = ik_ik = 0.0;
          double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
          double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
          double eri_aa = modprms.eri[a*syst.Number_of_atoms + a];

          get_integrals(i,k,basis_ao,eri_aa,G1,F2,ii_kk,ik_ik);


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
        for(b=0;b<syst.Number_of_atoms;b++){

          if(b!=a){

            // Compute density matrix due to all orbitals on atom b 
            el->Fao_alp->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*modprms.eri[a*syst.Number_of_atoms+b]; 


            // Derivatives

            VECTOR deri, dV_AB;    
            compute_indo_core_parameters_derivs(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map,
                                                sorb_indx, modprms.indo_opt, a, b, c, deri, dV_AB);

            dFao_alp_dx->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.x;
            dFao_alp_dy->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.y;
            dFao_alp_dz->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.z;

            dFao_bet_dx->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.x;
            dFao_bet_dy->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.y;
            dFao_bet_dz->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.z;


          }
        }// for b
  
  
      }// i==j
      else{      // Off-diagonal terms
  
        if(a==b){ // different orbitals are on the same atom
          double ij_ij,ii_jj; ij_ij = ii_jj = 0.0;

          double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
          double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
          double eri_aa = modprms.eri[a*syst.Number_of_atoms+a];

          get_integrals(i,j,basis_ao,eri_aa,G1,F2,ii_jj,ij_ij);

          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - 0.5*el->P->M[i*Norb+j])*ij_ij - 0.5*el->P->M[i*Norb+j]*ii_jj );
            el->Fao_bet->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - 0.5*el->P->M[i*Norb+j])*ij_ij - 0.5*el->P->M[i*Norb+j]*ii_jj );
          }
          else{  
            el->Fao_alp->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - el->P_alp->M[i*Norb+j])*ij_ij - el->P_alp->M[i*Norb+j]*ii_jj );
            el->Fao_bet->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - el->P_bet->M[i*Norb+j])*ij_ij - el->P_bet->M[i*Norb+j]*ii_jj );
          }
  
        }// a==b - different orbitals are on the same atom
  
        else{ // different orbitals are on different atoms

          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 

            VECTOR deri, dV_AB;    
            compute_indo_core_parameters_derivs(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map,
                                                sorb_indx, modprms.indo_opt, a, b, c, deri, dV_AB);

            dFao_alp_dx->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.x; 
            dFao_alp_dy->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.y; 
            dFao_alp_dz->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.z; 

            dFao_bet_dx->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.x; 
            dFao_bet_dy->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.y; 
            dFao_bet_dz->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.z; 

          }
          else{
            el->Fao_alp->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 


            VECTOR deri, dV_AB;    
            compute_indo_core_parameters_derivs(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map,
                                                sorb_indx, modprms.indo_opt, a, b, c, deri, dV_AB);

            dFao_alp_dx->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*deri.x; 
            dFao_alp_dy->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*deri.y; 
            dFao_alp_dz->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*deri.z; 

            dFao_bet_dx->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*deri.x; 
            dFao_bet_dy->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*deri.y; 
            dFao_bet_dz->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*deri.z; 

          }

        }
    
      }// i!=j
  
  
    }// for j
  }// for i

}

void Hamiltonian_Fock_derivs_indo
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int c, 
  MATRIX& dHao_dx,     MATRIX& dHao_dy,     MATRIX& dHao_dz,
  MATRIX& dFao_alp_dx, MATRIX& dFao_alp_dy, MATRIX& dFao_alp_dz,
  MATRIX& dFao_bet_dx, MATRIX& dFao_bet_dy, MATRIX& dFao_bet_dz
){

  Hamiltonian_Fock_derivs_indo( &el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map,
  c, &dHao_dx, &dHao_dy, &dHao_dz,  &dFao_alp_dx, &dFao_alp_dy, &dFao_alp_dz,  &dFao_bet_dx, &dFao_bet_dy, &dFao_bet_dz);

}





}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian

