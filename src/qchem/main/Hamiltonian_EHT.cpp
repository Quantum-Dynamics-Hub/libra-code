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

#include "Hamiltonian.h"

/****************************************************************************

  This file contains following functions:

  void Hamiltonian_core_eht(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                            vector<int>& fragment,vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                            MATRIX* Hao,MATRIX* Sao,Memory* mem)

  void Hamiltonian_Fock_eht(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                            vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                            Electronic* el,Electronic* el0, Memory* mem)




****************************************************************************/

void Hamiltonian_core_eht(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                          vector<int>& fragment,vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                          MATRIX* Hao,MATRIX* Sao,Memory* mem){
// prms - control parameters
// modprms - model parameters
// mol - nuclear geometry and charges
// basis_fo - basis of fragment obitals - basis_fo[i] - is index of basis_ao, which is i-th in the given group
// basis_ao - the full/complete pool of AOs (entire system)
// Hao, Sao - Hamiltonian and overlap matrices

// prms.eht_formula == 0 - unweighted 
// prms.eht_formula == 1 - weighted
// prms.eht_formula == 2 - Calzaferi
// prms.eht_formula == 3 - my developments (exhange parameters + long-range)

// prms.eht_sce_formula == 0 - EHT (no self-consistent electrostatics)
// prms.eht_sce_formula == 1 - SC-EHT (with self-consistent electrostatics based on atomic gross populations!)
// prms.eht_sce_formula == 2 - SC-EHT (with self-consistent electrostatics based on net orbital populations wrt reference value)

           
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


  //========================= Charge-corrected diagonal elements ==========================  
  // Depending of SCE, do different type of correction to diagonal elements
  if(prms.eht_sce_formula==0){ ; ; }
  else if(prms.eht_sce_formula==1){

    // Now modify IPs - diagonal elements
    for(int a=0;a<Norb;a++){
      int A = basis_fo[a];
      double Ai = 0.0;  // typically a positive number
      double Q = mol.Mull_charges_gross[basis_ao[A].at_indx];
      if(Q>0){ Ai = modprms.PT[basis_ao[A].element].J_param1[basis_ao[A].ao_shell]; }
      else{    Ai = modprms.PT[basis_ao[A].element].J_param2[basis_ao[A].ao_shell]; }

      Hao->M[a*Norb+a] -= (Ai * Q);

    }// for a
  }// eht_sce_formula == 1

  else if(prms.eht_sce_formula==2){

    // Now modify IPs - diagonal elements
    for(int a=0;a<Norb;a++){
      // Larsson & Pykko
/*
      double dn = (Mull_orb_pop_net[a] - Mull_orb_pop_net0[a]);
      double Ai = 0.0;

      if(dn>0.0){ Ai = PT[basis_ao[a].element].J_param1[basis_ao[a].ao_shell]; }
      else{ Ai = PT[basis_ao[a].element].J_param2[basis_ao[a].ao_shell]; }
      

      Hao->M[a*Norb+a] -= Ai * dn;
*/
    }// for a
  }// eht_sce_formula == 2

  else{
    cout<<"Warning (in System::build_core_eht):  no eht_sce_formula="<<prms.eht_sce_formula<<" is known. Skipping...\n";
  }


  //========================= Off-diagonal elements ==========================    
  double K_const = 0.0;// = 1.75;

  // Off-diagonal elements
  for(i=0;i<Norb;i++){
    I = basis_fo[i];
    for(j=i+1;j<Norb;j++){
        J = basis_fo[j];

        K_const = modprms.eht_k.get_K_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);

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

        }

        else if(prms.eht_formula==3){       // My developments

          double K1_const = modprms.eht_k.get_K1_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);
          double K2_const = modprms.eht_k.get_K2_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);
          double K3_const = modprms.eht_k.get_K3_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);
          double K4_const = modprms.eht_k.get_K4_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);

          double r2ab = (mol.R[basis_ao[I].at_indx] - mol.R[basis_ao[J].at_indx]).length2();

          Hao->M[i*Norb+j]  = 0.5*K_const*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j])*Sao->M[i*Norb+j];        
//          Hao->M[i*Norb+j] += K1_const;
          Hao->M[i*Norb+j] += K1_const*Sao->M[i*Norb+j]; // /sqrt(1.0 + r2ab);

          VECTOR dIdA,dIdB;
          Hao->M[i*Norb+j] += K2_const*KINETIC_INTEGRAL(basis_ao[I],basis_ao[J],dIdA,dIdB,mem->aux, mem->n_aux);


//          Hao->M[i*Norb+j] += K2_const*Sao->M[i*Norb+j]*Sao->M[i*Norb+j]; // /sqrt(1.0 + r2ab);

//          Hao->M[i*Norb+j] += K3_const*(Sao->M[i*Norb+j])/sqrt((K4_const*K4_const) + r2ab);


//          Hao->M[i*Norb+j] += K3_const*(Sao->M[i*Norb+j])/sqrt((K4_const*K4_const) + r2ab);
//          Hao->M[i*Norb+j] += K3_const*(log(1.1+Sao->M[i*Norb+j]))/sqrt((K4_const*K4_const) + r2ab);
          


          Hao->M[j*Norb+i] = Hao->M[i*Norb+j];

        }// == 3
     
    }// for bj
  }// for i



  // Correct diagonal terms after they have been used to construct off-diagonal terms:
  for(i=0;i<Norb;i++){
    I = basis_fo[i];

    if(prms.eht_formula==3){       // My developments

          double K1_const = modprms.eht_k.get_K1_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[I].element,basis_ao[I].ao_shell);
          double K2_const = modprms.eht_k.get_K2_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[I].element,basis_ao[I].ao_shell);
          double K3_const = modprms.eht_k.get_K3_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[I].element,basis_ao[I].ao_shell);
          double K4_const = modprms.eht_k.get_K4_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[I].element,basis_ao[I].ao_shell);

          double r2ab = (mol.R[basis_ao[I].at_indx] - mol.R[basis_ao[I].at_indx]).length2();

//          VECTOR dIdA,dIdB;
//          Hao->M[i*Norb+i] += K2_const*KINETIC_INTEGRAL(basis_ao[I],basis_ao[J],dIdA,dIdB,mem->aux, mem->n_aux);

//          Hao->M[i*Norb+i] += K1_const;
//          Hao->M[i*Norb+i] += K2_const*Sao->M[i*Norb+i]; 
//          Hao->M[i*Norb+i] += K3_const*Sao->M[i*Norb+i]*Sao->M[i*Norb+i]; 
//          Hao->M[i*Norb+i] += K3_const*(log(1.0+Sao->M[i*Norb+i]))/sqrt((K4_const*K4_const) + r2ab);


    }// eht_formula == 3

  }//



}

void Hamiltonian_Fock_eht(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                          vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                          Electronic* el,Electronic* el0, Memory* mem){

  // Update charges
  //*el->P = *el->P_alp + *el->P_bet;

  update_Mull_orb_pop(el->P, el->Sao, el->Mull_orb_pop_gross, el->Mull_orb_pop_net);

  // Here we need a local version of at_orbitals - only those atoms that carry orbitals spanning basis_fo
  // we assume that there is no overlaps of atoms between two subsystems

  update_Mull_charges(fragment, basis_fo, at_orbitals, mol.Zeff, el->Mull_orb_pop_gross, el->Mull_orb_pop_net,mol.Mull_charges_gross, mol.Mull_charges_net);


  // Recompute (core) EHT Hamiltonian - it may be charge-dependent (well, this is not strictly speaking a "core")
  Hamiltonian_core_eht(prms,modprms,mol,fragment,basis_fo,basis_ao, at_orbitals, el->Hao,el->Sao,mem);

  // Compute Fock matrices
  if(prms.eht_fock_opt==0){  // no self-consistency correction is needed, Hamiltonian itself is a Fock matrix

    *el->Fao_alp = *el->Hao;
    *el->Fao_bet = *el->Hao;

  }

  else if(prms.eht_fock_opt==1){  // with self-consistency correction

    *el->Fao_alp = 2.0 * (*el->Hao) - (*el0->Hao); 
    *el->Fao_bet = 2.0 * (*el->Hao) - (*el0->Hao); 
  }


  // On top - add possible electrostatics: field due to charges
  if(prms.eht_electrostatics==1){ // external field for effective 1-electron problem (Fock)

    vector<double> Wii(el->Norb,0.0);

    // Precompute electrostatic correction factor for all orbitals
    for(int i=0;i<el->Norb;i++){
      int I = basis_fo[i];

      double K3_const = modprms.eht_k.get_K3_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[I].element,basis_ao[I].ao_shell);
      double K4_const = modprms.eht_k.get_K4_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[I].element,basis_ao[I].ao_shell);
        
      Wii[i] = 0.0;

      for(int k=0;k<fragment.size();k++){  
        int n = fragment[k];      

        if(basis_ao[I].at_indx!=n){  // Contributions from all other atoms
          double dist = (mol.R[basis_ao[I].at_indx]-mol.R[n]).length();

          Wii[i] -= K3_const*mol.Mull_charges_gross[n]/sqrt((K4_const*K4_const) + dist*dist);

        }// if - different atoms

        
        // Contribution from periodic images - well, this is not really converging, but at least first environment will be
        // taken into account
        for(int nx=-prms.x_period;nx<=prms.x_period;nx++){
          for(int ny=-prms.y_period;ny<=prms.y_period;ny++){
            for(int nz=-prms.z_period;nz<=prms.z_period;nz++){
        
              if(nx==0 && ny==0 && nz==0){ }
              else{
                
                VECTOR TV = nx*prms.t1 + ny*prms.t2 + nz*prms.t3;

                double dist = (mol.R[basis_ao[I].at_indx]-mol.R[n]-TV).length();
                // This summation corresponds to k = 0 (Gamma-point)

                Wii[i] -= K3_const*mol.Mull_charges_gross[n]/sqrt((K4_const*K4_const) + dist*dist);

                //Sao->M[i*Norb+j] += OVERLAP_INTEGRAL(basis_ao[I],basis_ao[J],0,dIdA,dIdB,aux,n_aux,TV);
              }
        
            }// for nz
          }// for ny
        }// for nx 


      }// for k
    }// for i




    // Now use diagonal correction factors to update all elements of Fock matrix
    for(int i=0;i<el->Norb;i++){
      for(int j=0;j<el->Norb;j++){

        el->Fao_alp->M[i*el->Norb+j] += 0.5*(Wii[i]+Wii[j])*el->Sao->M[i*el->Norb+j];
        el->Fao_bet->M[i*el->Norb+j] += 0.5*(Wii[i]+Wii[j])*el->Sao->M[i*el->Norb+j];
 
      }// j
    }// i
    

  }// eht_electrostatics == 1



  else if(prms.eht_electrostatics==2){ // external field for effective 1-electron problem (Fock)

    VECTOR da,db,dc;

    for(int i=0;i<el->Norb;i++){
      int I = basis_fo[i];
//      AO ao_i(basis_ao[I]); ao_i.x_exp = 0; ao_i.y_exp = 0; ao_i.z_exp = 0;
//      for(int p=0;p<ao_i.primitives.size();p++){ ao_i.primitives[p].x_exp = ao_i.primitives[p].y_exp = ao_i.primitives[p].z_exp = 0; }


      for(int j=i+1;j<el->Norb;j++){
        int J = basis_fo[j];

          double Wij = 0.0;


//          AO ao_j(basis_ao[J]); ao_j.x_exp = 0; ao_j.y_exp = 0; ao_j.z_exp = 0;
//          for(int p1=0;p1<ao_j.primitives.size();p1++){ ao_j.primitives[p1].x_exp = ao_j.primitives[p1].y_exp = ao_j.primitives[p1].z_exp = 0; }


          double gamma = basis_ao[I].primitives[0].G_alpha + basis_ao[J].primitives[0].G_alpha;
          VECTOR R_AB; R_AB = *basis_ao[I].primitives[0].G_center - *basis_ao[J].primitives[0].G_center;

          VECTOR P_AB; P_AB = (basis_ao[I].primitives[0].G_alpha*(*basis_ao[I].primitives[0].G_center)
                             + basis_ao[J].primitives[0].G_alpha*(*basis_ao[J].primitives[0].G_center))/gamma;



          double K3_const = modprms.eht_k.get_K3_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);
          double K4_const = modprms.eht_k.get_K4_value(basis_ao[I].element,basis_ao[I].ao_shell,basis_ao[J].element,basis_ao[J].ao_shell);

//          cout<<"I= "<<I<<" J= "<<J<<" looking: "<<basis_ao[I].element<<" "<<basis_ao[I].ao_shell<<" | "<<basis_ao[J].element<<" "<<basis_ao[J].ao_shell<<endl;
//          cout<<" Found result is: K3 = "<<K3_const/eV<<" eV, K4 = "<<K4_const/Angst<<" Angst\n";


          for(int k=0;k<fragment.size();k++){  

              int n = fragment[k];      

//            if(basis_ao[I].at_indx!=n && basis_ao[J].at_indx!=n){ // whith this condition SiH4 is symmetric, but bodypy has unreasonable charges
                                                                  // without this - bodypy has reasonable charges, but SiH4 is broken symmetry (non-zero dipole moment)
              
//            Wij -= mol.Mull_charges_gross[n] * 
//                   NUCLEAR_ATTRACTION_INTEGRAL(basis_ao[I],basis_ao[J],mol.R[n],n,da,db,dc,mem->aux,mem->n_aux,mem->auxv,mem->n_auxv);


//            Wij -= mol.Mull_charges_gross[n] * 
//                   NUCLEAR_ATTRACTION_INTEGRAL(ao_i,ao_j,mol.R[n],n,da,db,dc,mem->aux,mem->n_aux,mem->auxv,mem->n_auxv);


              double dist = (P_AB-mol.R[n]).length();

//              if(dist>0.0){
//                Wij -= mol.Mull_charges_gross[n]/dist;
//              }
//              else{
//                cout<<"i= "<<i<<" j= "<<j<<" k= "<<k<<" dist= "<<dist<<endl;
//              }

          Wij -= K3_const*mol.Mull_charges_gross[n]/sqrt((K4_const*K4_const) + dist*dist);


//            }//if


          }// for k

//          Wij *= (el->Sao->M[i*el->Norb+j]); 



          double scl = 1.0;
          if(i==j){ scl = 0.0; }

          el->Fao_alp->M[i*el->Norb+j] += scl*Wij;
          el->Fao_bet->M[i*el->Norb+j] += scl*Wij;
                                         
          el->Fao_alp->M[j*el->Norb+i] += scl*Wij;
          el->Fao_bet->M[j*el->Norb+i] += scl*Wij;

          cout<<"correction to F_["<<I<<"]["<<J<<"] = "<<scl*Wij<<"  K3= "<<K3_const<<" K4= "<<K4_const<<endl;


      }// j
    }// i
  }// external_filed
    
}

