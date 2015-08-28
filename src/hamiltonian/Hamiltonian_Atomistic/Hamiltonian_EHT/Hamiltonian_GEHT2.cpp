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
#include "MOAO.h"

/****************************************************************************

  This file contains following functions:

  void Hamiltonian_core_geht2(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                              vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao, vector<vector<int> >& at_orbitals, 
                              MATRIX* Hao,MATRIX* Sao,Memory* mem, vector<double>& eri, vector<double>& V_AB)


  void Hamiltonian_Fock_geht2(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                              vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                              Electronic* el, Memory* mem,vector<double>& eri, vector<double>& V_AB)



****************************************************************************/


void Hamiltonian_core_geht2(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                            vector<int>& fragment,vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                            MATRIX* Hao,MATRIX* Sao,Memory* mem){

// Generalized EHT

// prms - control parameters
// modprms - model parameters
// mol - nuclear geometry and charges
// basis_fo - basis of fragment obitals - basis_fo[i] - is index of basis_ao, which is i-th in the given group
// basis_ao - the full/complete pool of AOs (entire system)
// Hao, Sao - Hamiltonian and overlap matrices


           
  int i,j,n, I,J, a,b;
  double delt, delt2, delt4;
  VECTOR dIdA, dIdB;
  

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
  for(int i=0;i<Norb;i++){
    int I = basis_fo[i];          // global index of orbital i
    int a = basis_ao[I].at_indx;  // global index of the atom on which orbital i is centered


    for(int j=i+1;j<Norb;j++){
      int J = basis_fo[j];          // global index of orbital j
      int b = basis_ao[J].at_indx;  // global index of the atom on which orbital j is centered

      K_const = modprms.meht_k.get_K_value(I,J);

      double dist = (mol.R[a]-mol.R[b]).length();

      if(dist<200.0){  // cutoff in a.u.

        //===============  Core EHT (tight-binding) contribution  =======================

        if(prms.eht_formula==0){  // Unweighted formula
  
          Hao->M[i*Norb+j] = 0.5*K_const*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j])*Sao->M[i*Norb+j]; 
  
        }
  
        else if(prms.eht_formula==1){  // Weighted formula:        
  
          delt = (Hao->M[i*Norb+i]-Hao->M[j*Norb+j])/(Hao->M[i*Norb+i]+Hao->M[j*Norb+j]);
          delt2 = delt*delt;
          delt4 = delt2*delt2;
        
          Hao->M[i*Norb+j] = 0.5*(K_const + delt2 + (1.0 - K_const)*delt4)*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j])*Sao->M[i*Norb+j];
  
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
  
        }// eht_formula==2 (Calzaferi)
  
        else if(prms.eht_formula==3){       // This is slightly generalized eht_formula == 0
  
          double K1_const = modprms.meht_k.get_K1_value(I,J); //!!!! Careful when using with self-consistent electrostatics !!!!
          Hao->M[i*Norb+j]  = (0.5*K_const*(Hao->M[i*Norb+i]+Hao->M[j*Norb+j]) )*Sao->M[i*Norb+j] + K1_const;        
  
        }// == 3
  
  
        //===============  Pseudopotential contribution  =======================
//        double C0 = modprms.meht_k.get_C0_value(I,J);
//        double C2 = modprms.meht_k.get_C2_value(I,J);
//        double C4 = modprms.meht_k.get_C4_value(I,J); // in a.u. 
  
        for(int n=0;n<mol.Nnucl;n++){

          for(int t=0;t<modprms.meht_k.eht_PP0[n].size();t++){

            double Ca = modprms.meht_k.eht_PPa[n][t];
            double C0 = modprms.meht_k.eht_PP0[n][t];
            double C2 = modprms.meht_k.eht_PP2[n][t];

//            cout<<"n= "<<n<<" t= "<<t<<" Ca= "<<Ca<<" C0= "<<C0<<" C2= "<<C2<<endl;

//            exit(0);


//          if((mol.R[n]-mol.R[a]).length2()<10.0 || (mol.R[n]-mol.R[b]).length2()<10.0){

            VECTOR tv(0.0,0.0,0.0);
            double pd = PSEUDO_02_INTEGRAL(C0, C2, Ca, mol.R[n],basis_ao[I],basis_ao[J],0, dIdA, dIdB,mem->aux,mem->n_aux, tv );
  
            Hao->M[i*Norb+j] += pd; //(C0 + C1*dist + C2*dist*dist)*exp(-C4*C4*dist) * Sao->M[i*Norb+j];

//          }// if close

          }// for t
        }// for all atoms in the system
        //=========================================================================
  
  
      }// dist < 200.0
      else{

        Hao->M[i*Norb+j] = 0.0;

      }

      // Use of symmetry
      Hao->M[j*Norb+i] = Hao->M[i*Norb+j];

     
    }// for bj



//---------------------------------------------------------------
    for(int n=0;n<mol.Nnucl;n++){

      for(int t=0;t<modprms.meht_k.eht_PP0[n].size();t++){

        double Ca = modprms.meht_k.eht_PPa[n][t];
        double C0 = modprms.meht_k.eht_PP0[n][t];
        double C2 = modprms.meht_k.eht_PP2[n][t];


//      if((mol.R[n]-mol.R[a]).length2()<10.0 || (mol.R[n]-mol.R[b]).length2()<10.0){

        VECTOR tv(0.0,0.0,0.0);
        double pd = PSEUDO_02_INTEGRAL(C0, C2, Ca, mol.R[n],basis_ao[I],basis_ao[I],0, dIdA, dIdB,mem->aux,mem->n_aux, tv );
  
        Hao->M[i*Norb+i] += pd; //(C0 + C1*dist + C2*dist*dist)*exp(-C4*C4*dist) * Sao->M[i*Norb+j];

//      }// if close

      }// for t
    }// for all atoms in the system



  }// for i


}


// Here we actually start with energy (taken from charge equilibration) and then derive the correct Fock matrix
// no additional Fock matrix corrections are needed
// 
void Hamiltonian_Fock_geht2(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
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

  
/*
  for(int i=0;i<modprms.eht_k.psps_data.size();i++){

    // P^n1 * S^n2 * P^n3 * S^n4

    int n1 = modprms.eht_k.psps_data[i].n1;
    int n2 = modprms.eht_k.psps_data[i].n2;
    int n3 = modprms.eht_k.psps_data[i].n3;
    int n4 = modprms.eht_k.psps_data[i].n4;
    double K = modprms.eht_k.psps_data[i].K;

    // Inappropriate use of the arrays, but save on time
    int i1; 

    el->dFao_alp_dP_alp->Init_Unit_Matrix(1.0);
    el->dFao_alp_dP_bet->Init_Unit_Matrix(1.0);

    for(i1=0;i1<n1;i1++){
      *el->dFao_alp_dP_alp =  *el->dFao_alp_dP_alp * *el->P_alp; 
      *el->dFao_alp_dP_bet =  *el->dFao_alp_dP_bet * *el->P_bet; 
    }

    for(i1=0;i1<n2;i1++){
      *el->dFao_alp_dP_alp =  *el->dFao_alp_dP_alp * *el->Sao; 
      *el->dFao_alp_dP_bet =  *el->dFao_alp_dP_bet * *el->Sao; 
    }

    for(i1=0;i1<n3;i1++){
      *el->dFao_alp_dP_alp =  *el->dFao_alp_dP_alp * *el->P_alp; 
      *el->dFao_alp_dP_bet =  *el->dFao_alp_dP_bet * *el->P_bet; 
    }

    for(i1=0;i1<n4;i1++){
      *el->dFao_alp_dP_alp =  *el->dFao_alp_dP_alp * *el->Sao; 
      *el->dFao_alp_dP_bet =  *el->dFao_alp_dP_bet * *el->Sao; 
    }


    *el->Fao_alp += K * *el->dFao_alp_dP_alp;
    *el->Fao_bet += K * *el->dFao_alp_dP_bet;



  }// for i

  *el->dFao_alp_dP_alp = 0.0;
  *el->dFao_alp_dP_bet = 0.0;

*/

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



/*
      double dist = (mol.R[a]-mol.R[b]).length();

      double K1_const = modprms.meht_k.get_K1_value(I,J); // a.u. of energy - exchange integral
      double K2_const = modprms.meht_k.get_K2_value(I,J); // a.u. of energy - exchange integral
      double K3_const = modprms.meht_k.get_K3_value(I,J) * (Angst / eV); // convert from a.u. of energy to a.u. of length


      double K4_const = modprms.meht_k.get_K4_value(I,J); // a.u. of length    
      if(K4_const<0.0){  K4_const = 0.0; }      
      double f = ERF(K4_const * dist);


      double J_AB = (K1_const*(1.0-f)  - f/sqrt(dist*dist + K3_const*K3_const) );


*/


      //***************************************************
      if(a==b){  // orbitals are on the same atom - add fixed exchange to all orbitals

        // Exchange effects - only local, to preserve rotational invariance
        double K1_const = modprms.meht_k.get_K1_value(I,J); // a.u. of energy - exchange integral

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

/* 
        if(prms.eht_electrostatics>=2){

          for(int B=0;B<mol.Nnucl;B++){  // over all atoms - effects of Coulomb terms

            if(B!=A){

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
            
            
              double dist = (mol.R[A]-mol.R[B]).length();

              double K3_const = modprms.meht_k.get_K3_value(orb_A,orb_B) * (Angst / eV); // convert from a.u. of energy to a.u. of length
              double J_AB = 1.0/sqrt(dist*dist + K3_const*K3_const);
            
            
              // These are additional degrees of freedom to use in more elaborate 
              // functional J_AB
       
              double K2_const = modprms.meht_k.get_K2_value(orb_A,orb_B); // a.u. of energy
              double K4_const = modprms.meht_k.get_K4_value(orb_A,orb_B); // a.u. of length
       
              if(K4_const<0.0){  K4_const = 0.0; }
       
              double f = ERF(K4_const * dist);
              J_AB = J_AB * f +  (1.0 - f) * K2_const;  // add long-range electrostatics smoothly
       
            
              //***************** This corresponds to contribution to energy:  J_AB * Q_A * Q_B ****************
              el->Fao_alp->M[i*Norb+j] += 0.5*J_AB * (dQA_dPij_alp * QB + QA * dQB_dPij_alp);
              el->Fao_bet->M[i*Norb+j] += 0.5*J_AB * (dQA_dPij_bet * QB + QA * dQB_dPij_bet);
              //************************************************************************************************
          
            }// B!=A
          
          }// for B

        }// electrostatics >=2
*/
     

      }// for A


    }// for j
  }// i

  }// if eht_electrostatics == 1



    
}


